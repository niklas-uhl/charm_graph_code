#include "NDMeshStreamer.h"
#include "TopoManager.h"
#include "htram_group.h"
#include "weighted_htram_smp.decl.h"
//#define PAPI
#ifdef PAPI
#include <papi.h>
#endif
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <map>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <limits>
#include <queue>
#include <algorithm>
#include <random>

#define INFO_PRINTS
//#define PRINT_HISTO //print histograms to file
#define LOCAL_TO_TRAM //add all outgoing updates (even local) to tram
//#define PQ_HOLD_ONLY
//#define PQ_EDGE_DIST //add cost of smallest edge when finding bucket
//#define VCOUNT
//#define ALL_TO_TRAM_HOLD //place all updates in the tram hold at first
//#define NODE_LOAD_BALANCE

// set data type for messages
using tram_proxy_t = CProxy_HTram;
using tram_t = HTram;

/* readonly */
//tram_proxy_t tram_proxy;
CProxy_HTramRecv nodeGrpProxy;
CProxy_HTramNodeGrp srcNodeGrpProxy;
CProxy_Main mainProxy;
CProxy_SsspChares arr;
CProxy_SharedInfo shared;
CProxy_processHeapShared process_heap_shared;
int N;	  // number of processors
long V;	  // number of vertices
int M = 1024; // divisor for dest_table (must be power of 2)
long num_global_edges; //number of global edges in the graph (used when graph is generated)
long average_degree; //average degree of graph
int generate_mode; //0 = read from file, 1 = generate automatically
int S; //seed for randomization
cost lmax; // long maximum
int histo_bucket_count = 2048; // number of buckets for the histogram needed for message prioritization
#define HISTO_BUCKET_COUNT 2048 //needed macro for array init
int histo_reduction_width = histo_bucket_count/8;
double reduction_delay = 0.1; // each histogram reduction happens at this interval
int initial_threshold = 3; //initial histo threshold
// tram constants
int buffer_size = 1024; // meaningless for smp; size changed in htram_group.h
double flush_timer = 0.01; // milliseconds
bool enable_buffer_flushing = false; //true = buffer flushes at interval specified by flush_timer
tram_proxy_t tram_proxy;


void fast_exit(void *obj, double time);

void start_reductions(void *obj, double time)
{
	arr.contribute_histogram(0);
}

struct ComparePairs
{
	bool operator()(const Update &lhs, const Update &rhs) const
	{
		// Compare the second integers of the pairs
		return lhs.distance > rhs.distance; // '>' for min heap, '<' for max heap
	}
};

struct histoInstance
{
public:
  int fnz;
  int width;
  int *reducedValues;
};
  
class histogramSequence
{
private:
  int maxBuckets;
  std::vector<histoInstance> histos;
public: 

  histogramSequence(int _maxBuckets)
  {
    maxBuckets = _maxBuckets;
    histoInstance h;
    h.fnz = 0;
    h.width = 10;
    h.reducedValues = new int[10];
    histos.push_back(h);
  }

  void insert(int fnz, int width, long* histo)
  {
    histoInstance h;
    h.fnz = fnz;
    h.width = width;
    h.reducedValues = new int[width];
    for (int i = 0; i< width; i++) h.reducedValues[i] = histo[i];
    histos.push_back(h);
  }

  void putout()
  {
    // using cout instead of ckout to avoid buffer overflow. (should be output to a file)
	std::ofstream out_file;
	out_file.open("histos.txt");
    for (int i=0; i<histos.size(); i++)
      {
        for (int j=0; j<histos[i].fnz; j++) out_file << "0 ";
		for (int j=0; j<histos[i].width; j++) out_file << histos[i].reducedValues[j] << " ";
		for (int j=0; j < ( maxBuckets - histos[i].fnz - histos[i].width); j++) out_file << "0 "; 
		out_file << endl;
      }
  }
};


class Main : public CBase_Main
{
private:
	long start_vertex;
	long *partition_index;
	double start_time;
	double read_time;
	double total_time;
	long max_index;
	int threshold_change_counter;
	int previous_threshold;
	int reduction_counts = 0;
	int no_incoming = 0;
	std::vector<double> reduction_times;
	bool first_qd_done = false;
	bool second_qd_done = false;
	int activeBucketMax = 10;
	int current_phase = 0; //0=initial, 1=bfs, 2=converged_bfs
	int last_first_nonzero = 0;
	long previous_updates_created = 0;
	long previous_updates_processed = 0;
	long previous_distance_changes = 0;
	double tram_percentile = 0.01;
	double heap_percentile = 0.01;
	#ifdef PRINT_HISTO
	histogramSequence *histoSeq;
	#endif

public:

	double compute_begin;
	double compute_time;

	/**
	 * Read in graph from csv (currently sequential)
	 */
	Main(CkArgMsg *m)
	{
		N = CkNumPes();
		if (!m->argv[1])
		{
			ckout << "Missing vertex count" << endl;
			CkExit(0);
		}
		V = atol(m->argv[1]); // read in number of vertices
		if (!m->argv[2])
		{
			ckout << "Missing file name/edge count" << endl;
			CkExit(0);
		}
		std::string file_name = m->argv[2]; // read file name or edge count
		if (!m->argv[3])
		{
			ckout << "Missing random seed" << endl;
			CkExit(0);
		}
		S = atoi(m->argv[3]); // randomize seed
		if (!m->argv[4])
		{
			ckout << "Missing start vertex" << endl;
			CkExit(0);
		}
		start_vertex = atol(m->argv[4]); // number of beginning vertex
		if (!m->argv[5])
		{
			ckout << "Missing generate mode" << endl;
			CkExit(0);
		}
		generate_mode = atoi(m->argv[5]); //0 for read from csv, 1 for generation
		if (!m->argv[6])
		{
			ckout << "Missing tram percentile" << endl;
			CkExit(0);
		}
		tram_percentile = std::stod(m->argv[6]);
		if (!m->argv[7])
		{
			ckout << "Missing heap percentile" << endl;
			CkExit(0);
		}
		heap_percentile = std::stod(m->argv[7]);
		#ifdef PRINT_HISTO
		histoSeq = new histogramSequence(histo_bucket_count);
		#endif
		// create TRAM proxy
		nodeGrpProxy = CProxy_HTramRecv::ckNew();
		srcNodeGrpProxy = CProxy_HTramNodeGrp::ckNew();
		CkCallback ignore_cb(CkCallback::ignore);
		process_heap_shared = CProxy_processHeapShared::ckNew();
		tram_proxy = tram_proxy_t::ckNew(nodeGrpProxy.ckGetGroupID(), 
		srcNodeGrpProxy.ckGetGroupID(), buffer_size, enable_buffer_flushing, flush_timer, false, true, ignore_cb);
		shared = CProxy_SharedInfo::ckNew();
		arr = CProxy_SsspChares::ckNew(tram_proxy, N);
		mainProxy = thisProxy;
		arr.initiate_pointers();
		partition_index = new long[N + 1]; // last index=maximum index
		lmax = std::numeric_limits<cost>::max();
		start_time = CkWallTimer();
		if(generate_mode==2)
		{
			long side_length = (int) std::sqrt((double) V);
			#ifdef INFO_PRINTS
			ckout << "Side length: " << side_length << endl;
			#endif
			num_global_edges = 4 * side_length * (side_length - 1); //will ignore command line input
			#ifdef INFO_PRINTS
			ckout << "2-D Graph will be automatically generated with " << V << " vertices and " << num_global_edges << " edges" << endl;
			#endif
			for(int i=0; i<N+1; i++)
			{
				partition_index[i] = i * (V / N);
				if(i==N) partition_index[i] = V;
			}
			arr.generate_2d_graph(partition_index, N+1);
		}
		else if(generate_mode==1)
		{
			num_global_edges = std::stol(file_name);
			#ifdef INFO_PRINTS
			ckout << "Graph will be automatically generated with " << V << " vertices" << endl;
			#endif
			average_degree = num_global_edges / V;
			//for each pe, generate a random vertex and edge count, and send to pes
			long remaining_vertices = V;
			long current_start_index = 0; //tracks start vertex for indices 
			std::mt19937 generator(S);
			std::uniform_int_distribution<long> edge_count_distribution(((num_global_edges * 4) / (N * 5)), (num_global_edges * 6 / (N * 5))); //average += 20%
			std::uniform_int_distribution<long> vertex_count_distribution(((V * 4) / (N * 5)), ((V * 6) / (N * 5))); //average += 20%
			long *vertex_counts = new long[N];
			long *edge_counts = new long[N];
			for(int i=0; i<N; i++)
			{
				partition_index[i] = current_start_index;
				long vertex_count = vertex_count_distribution(generator);
				long edge_count = edge_count_distribution(generator);
				if((i == N - 1) && (vertex_count > remaining_vertices)) vertex_count = remaining_vertices; //make sure num_vertices = V
				remaining_vertices -= vertex_count;
				vertex_counts[i] = vertex_count;
				edge_counts[i] = edge_count;
				current_start_index += vertex_count;
			}
			partition_index[N] = V;
			#ifdef INFO_PRINTS
			ckout << "Partition index: ["; 
			for(int i=0; i<N+1; i++)
			{
				ckout << partition_index[i] << ", ";
			}
			ckout << "]" << endl;
			#endif
			for(int i=0; i<N; i++)
			{
				arr[i].generate_local_graph(vertex_counts[i], edge_counts[i], partition_index, N+1);
			}
		}
		else
		{
			#ifdef INFO_PRINTS
			ckout << "Graph will be read from file" << endl;
			#endif
			unsigned int seed = (unsigned int)S;
			srand(seed);
			// read file
			std::ifstream file(file_name);
			std::string readbuf;
			std::string delim = ",";
			// iterate through edge list
			CkVec<LongEdge> edges;
			int *incoming_count = new int[V]; //how many times does each vertex appear in the edge list
			for(int i=0; i<V; i++) incoming_count[i] = 0;
			max_index = 0; //maximum vertex index
			long edges_read = 0; //number of edges read
			while (getline(file, readbuf))
			{
				// get nodes on each edge
				std::string token = readbuf.substr(0, readbuf.find(delim));
				std::string token2 = readbuf.substr(readbuf.find(delim) + 1, readbuf.length());
				// make random distance
				cost edge_distance = (cost) rand() % 1000 + 1;
				// string to int
				long node_num = std::stol(token); //v
				long node_num_2 = std::stol(token2); //w
				incoming_count[node_num_2]++;
				// find the maximum vertex index
				if (node_num > max_index)
					max_index = node_num;
				if (node_num_2 > max_index)
					max_index = node_num_2;
				LongEdge new_edge;
				new_edge.begin = node_num;
				new_edge.end = node_num_2;
				new_edge.distance = edge_distance;
				edges.insertAtEnd(new_edge);
				edges_read++;
				//if(edges_read%10000000==0) ckout << "Read " << edges_read << " edges" << endl;
				// ckout << "One loop iteration complete" << endl;
			}
			average_degree = edges_read / V;
			for(int i=0; i<max_index; i++)
			{
				if(incoming_count[i]==0) no_incoming++;
			}
			#ifdef INFO_PRINTS
			ckout << "Vertices with no incoming edges: " << no_incoming << endl;
			ckout << "Max index: " << max_index << endl;
			#endif
			//ckout << "Loop complete" << endl;
			file.close();
			read_time = CkWallTimer() - start_time;
			// assign nodes to location
			std::vector<LongEdge> edge_lists[N];
			long average = edges.size() / N;
			for (int i = 0; i < edges.size(); i++)
			{
				int dest_proc = i / average;
				if (dest_proc >= N)
					dest_proc = N - 1;
				else if (i % average == 0)
					partition_index[dest_proc] = edges[i].begin;
				edge_lists[dest_proc].insert(edge_lists[dest_proc].end(), edges[i]);
			}
			partition_index[N] = max_index + 1;
			// reassign edges to move to correct pe
			for (int i = 0; i < N - 1; i++)
			{
				for (int j = edge_lists[i].size() - 1; j >= 0; --j)
				{
					// TODO
					if (edge_lists[i][j].begin >= partition_index[i + 1])
					{
						edge_lists[i + 1].insert(edge_lists[i + 1].begin(), edge_lists[i][j]);
						edge_lists[i].erase(edge_lists[i].begin() + j);
					}
				}
			}
			// add nodes to node lists
			// send subgraphs to nodes
			for (int i = 0; i < N; i++)
			{
				arr[i].get_graph(edge_lists[i].data(), edge_lists[i].size(), partition_index, N + 1);
			}
		}
	}

	/**
	 * Start algorithm from source vertex
	 */
	void begin(cost max_sum)
	{
		// ready to begin algorithm
		shared.max_path_value(max_sum);
		if(generate_mode==1||generate_mode==2) read_time = CkWallTimer() - start_time;
		#ifdef INFO_PRINTS
		ckout << "The sum of the maximum out-edges is " << max_sum << endl;
		#endif
		Update new_edge;
		new_edge.dest_vertex = start_vertex;
		new_edge.distance = 0;
		int dest_proc = 0;
		for (int i = 0; i < N; i++)
		{
			if (start_vertex >= partition_index[i] && start_vertex < partition_index[i + 1])
			{
				dest_proc = i;
				break;
			}
			if (i == N - 1)
				dest_proc = N - 1;
		}
		// quiescence detection
		//CkCallback cb(CkIndex_Main::quiescence_detected(), mainProxy);
		//CkStartQD(cb);
		// temp callback to test flushing
		threshold_change_counter = 0;
		previous_threshold = initial_threshold;
		CcdCallFnAfter(start_reductions, (void *) this, reduction_delay);
		CcdCallFnAfter(fast_exit, (void *) this, 30000.0); //end after 5 s
		compute_begin = CkWallTimer();
		#ifdef INFO_PRINTS
		ckout << "Beginning at time: " << compute_begin << endl;
		#endif
		arr.start_papi();
		arr[dest_proc].start_algo(new_edge);
	}

	/**
	 * Before printing distances, check if all the buffers are empty
	 * If not, flush the buffer (allowing the execution to continue)
	 * also restart qd
	 * If empty, end execution by printing the distances
	 */
	void quiescence_detected()
	{
		/*
		if (first_qd_done == false) 
		{
			first_qd_done = true;
			#ifdef INFO_PRINTS
			ckout << "First Quiescence detected at time: " << CkWallTimer() << endl;
			#endif
			CkCallback cb(CkIndex_Main::quiescence_detected(), mainProxy);
			CkStartQD(cb);
			// Ask everyone to call flush
			arr.get_bucket_limit(initial_threshold, initial_threshold+2, 0); // dummy values.. hopefully no damage. Just get tflush called
 	  	}
 	  	else 
		{
			second_qd_done = true;
			#ifdef INFO_PRINTS
			ckout << "Second Quiescence detected at time: " << CkWallTimer() << "  starting reductions. "<< endl;
			#endif
			arr.contribute_histogram(0);
 	  	}
		*/
	}


	/**
	 * Receive histo values from pes
	 * The idea is to get the distribution of update values, then 
	 * do local processing to select thresholds
	*/
	void reduce_histogram(long *histo_values, int histo_length)
	{
		reduction_times.push_back(CkWallTimer());
		reduction_counts++;
		long histogram_sum = 0;
		int first_nonzero = -1;
		long updates_processed = histo_values[histo_reduction_width + 1];
		long updates_created = histo_values[histo_reduction_width];
		long bfs_processed = histo_values[histo_reduction_width+2];
		long done_vertex_count = histo_values[histo_reduction_width+3];
		long updates_noted = histo_values[histo_reduction_width+4];
		long bfs_noted = histo_values[histo_reduction_width+5];
		long distance_changes = histo_values[histo_reduction_width+6];
		int heap_threshold = 0;
		int tram_threshold = 0;
		int bfs_threshold = heap_threshold;
		long active_counter = 0;
		//calculate the total histogram sum
		for(int i=0; i<histo_reduction_width; i++)
		{
			histogram_sum += histo_values[i];
			if((histo_values[i]>0)&&(first_nonzero==-1))
			{
				first_nonzero = i + last_first_nonzero;
			}
		}
		#ifdef PRINT_HISTO
		histoSeq->insert(last_first_nonzero, histo_reduction_width, histo_values);
		#endif
		#ifdef INFO_PRINTS
		ckout << "Updates: created: " << updates_created << ", noted: " << updates_noted << ", processed: " << updates_processed << 
		", distance changes: " << distance_changes << ", Done vertices: " << done_vertex_count << ", BFS noted: " << bfs_noted;
		#endif
		/*
		if ((bfs_processed == done_vertex_count) && (bfs_processed > 1000)) // not quite correct.. this should be affter we are sure bfs_processed has converged.. maybe via qd
		{
		    ckout << "all reachable vertices are done. " << bfs_processed << ":" << done_vertex_count << " at time: " << CkWallTimer() << endl; 
		    compute_time = CkWallTimer() - compute_begin;
		    arr.print_distances();
		    return;
		}
		*/
		if((updates_processed-updates_created==1)&&(updates_created>1000)&&(updates_created==previous_updates_created)&&(updates_processed==previous_updates_processed))
		{
			ckout << endl << "updates_processed and updates_created match" << endl;
			#ifdef INFO_PRINTS
			ckout << "Threshold: " << previous_threshold << endl;
			#endif
			compute_time = CkWallTimer() - compute_begin;
			arr.print_distances();
			return;
		}
		previous_updates_created = updates_created;
		previous_updates_processed = updates_processed;
		//calculate target percentile
		double heap_percent; //heap percentage
		double tram_percent; //tram percentage
		if(histogram_sum <= N * 100) 
		{
			heap_percent = 0.9999;
			tram_percent = 0.9999;
		}
		else
		{
			heap_percent = heap_percentile;
			tram_percent = tram_percentile;
		}
		previous_distance_changes = distance_changes;
		previous_updates_processed = updates_processed;
		//select bucket limit
		for(int i=0; i<histo_reduction_width; i++)
		{
			active_counter += histo_values[i];
			if((double) active_counter >= histogram_sum * heap_percent) 
			{
				heap_threshold = i + last_first_nonzero;
				break;
			}
			
		}
		bfs_threshold = heap_threshold;
		active_counter = 0;
		for(int i=0; i<histo_reduction_width; i++)
		{
			active_counter += histo_values[i];
			if((double) active_counter >= histogram_sum * tram_percent)
			{
				tram_threshold = i + last_first_nonzero;
				break;
			}
		}
		//in case of floating point weirdness
		if(heap_threshold >= histo_bucket_count) heap_threshold = histo_bucket_count - 1;
		if(tram_threshold >= histo_bucket_count) tram_threshold = histo_bucket_count - 1;
		if(histogram_sum==0)
		{
			heap_threshold = histo_bucket_count - 1;
			tram_threshold = histo_bucket_count - 1;
			bfs_threshold = histo_bucket_count - 1;
		}
		if(heap_threshold != previous_threshold)
		{
			previous_threshold = heap_threshold;
			threshold_change_counter++;
		}
		#ifdef INFO_PRINTS
		ckout << ", Heap threshold: " << heap_threshold << ", Tram: " << 
		tram_threshold << ", BFS threshold: " << bfs_threshold << ", first nonzero: " << first_nonzero << ", t= " << CkWallTimer() << endl;
		#endif
		//arr.contribute_histogram(first_nonzero-1);
		last_first_nonzero = first_nonzero;
		arr.current_thresholds(heap_threshold, tram_threshold, bfs_threshold, first_nonzero - 1, current_phase);
		
		//start next reduction round
		//CcdCallFnAfter(start_reductions, (void *) this, reduction_delay);
		
	}

	/**
	 * returns when all buffers are checked
	 */
	void check_buffer_done(long *msg_stats, int N)
	{
		/*
		int net_messages = msg_stats[1] - msg_stats[0]; // updates_processed - updates_created
		if (net_messages == 1)							// difference of 1 because of initial send
		{
			//ckout << "Real quiescence, terminate" << endl;
			compute_time = CkWallTimer() - compute_begin;
			//arr.stop_periodic_flush();
			arr.print_distances();
		}
		else
		{
			//ckout << "False quiescence, continue execution" << endl;
			CkCallback cb(CkIndex_Main::quiescence_detected(), mainProxy);
			CkStartQD(cb);
			arr.keep_going();
		}
		*/
	}

	void done(long *msg_stats, int N)
	{
		// ends program, prints that program is ended
		// ckout << "Completed" << endl;
		// CkPrintf("Memory usage at end: %f\n", CmiMemoryUsage()/(1024.0*1024.0));
		total_time = CkWallTimer() - start_time;
		ckout << "Actual edges: " << msg_stats[4+histo_bucket_count] << endl;
		ckout << "Read time: " << read_time << endl;
		ckout << "Compute time: " << compute_time << endl;
		ckout << "Total time: " << total_time << endl;
		ckout << "Wasted updates: " << msg_stats[0] - V << endl;
		ckout << "Wasted updates normalized to |E|: " << (double) (msg_stats[0] - V) / msg_stats[4+histo_bucket_count] << endl;
		ckout << "Rejected updates: " << msg_stats[1] << endl;
		ckout << "Rejected updates normalized to |E|: " << (double) msg_stats[1] / msg_stats[4+histo_bucket_count] << endl;
		ckout << "Number of threshold changes: " << threshold_change_counter << endl;
		ckout << "Number of reductions: " << reduction_counts << endl;
		ckout << "Updates noted: " << msg_stats[3+histo_bucket_count] << endl;
		ckout << "Distance changes: " << msg_stats[5+histo_bucket_count] << ", per vertex: " << msg_stats[5+histo_bucket_count] * 1.0 / V << endl;
		#ifdef VCOUNT
		long vcount_sum = 0;
		ckout << "Vcount: [ ";
		for(int i=0; i<histo_bucket_count+1; i++)
		{
			ckout << msg_stats[i+2] << ", ";
			vcount_sum += msg_stats[i+2];
		}
		ckout << endl;
		ckout << "Vcount sum: " << vcount_sum << endl;
		#endif
		#ifdef PAPI
		ckout << "Total insts: " << msg_stats[6+histo_bucket_count] << endl;
		ckout << "Insts per edge: " << msg_stats[6+histo_bucket_count] * 1.0 / msg_stats[4+histo_bucket_count] << endl;
		#endif
		arr.get_max_cost();
	}

	void done_max_cost(cost max_cost)
	{
		ckout << "Maximum vertex cost, not counting unreachable: " << max_cost << endl;
		#ifdef PRINT_HISTO
		histoSeq->putout();
		#endif
		CkExit(0);
	}
};

void fast_exit(void *obj, double time)
{
	ckout << "Ending program now at time " << CkWallTimer() << endl;
	((Main*) obj)->compute_time = CkWallTimer() - ((Main*) obj)->compute_begin;
	arr.print_distances();
	CkExit(0);
}

/**
 * This holds information that needs to be broadcasted
 * but that is calculated after the Main method
*/
class SharedInfo : public CBase_SharedInfo
{
	public:
	cost max_path;
	int event_id;
	int bracketed_id;
	int othercaller_id;
	std::vector<std::atomic_int *> chunks_remaining;

	SharedInfo()
	{
		event_id = traceRegisterUserEvent("Contrib reduction");
		bracketed_id = traceRegisterUserEvent("Vector inserts");
		othercaller_id = traceRegisterUserEvent("Other caller");
		#ifdef PAPI
		if(CkNodeFirst(CkMyNode())==CkMyPe())
		{
			int retval = PAPI_library_init(PAPI_VER_CURRENT);
			if(retval != PAPI_VER_CURRENT)
			{
				fprintf(stderr,"PAPI library init error!\n");
				CkExit(1);
			}
		}
		#endif
		for(int i=0; i<N; i++)
		{
			chunks_remaining.push_back(new std::atomic_int[HISTO_BUCKET_COUNT]);
			for (int j=0; j<histo_bucket_count; j++)
			{
				chunks_remaining[i][j] = 0;
			}
		}
	}

	void max_path_value(cost max_path_val)
	{
		max_path = max_path_val;
	}

};

/**
 * Array of chares for Dijkstra
 */
class SsspChares : public CBase_SsspChares
{
private:
	Node *local_graph; //structure to hold vertices assigned to this pe
	long start_vertex; //global index of lowest vertex assigned to this pe
	long num_vertices=0; //number of vertices assigned to this pe
	long updates_created_locally=0; //number of update messages sent
	long updates_processed_locally=0; //number of update messages received
	long *partition_index; //defines boundaries of indices for each pe
	long wasted_updates=0; //number of updates that don't have the final answer
	long rejected_updates=0; //number of updates that don't decrease a distance value/create more messages
	tram_proxy_t tram_proxy;
	tram_t *tram; //tram library
	SharedInfo *shared_local;
	std::priority_queue<Update, std::vector<Update>, ComparePairs> pq; //heap of messages
	long *histogram; //local histogram of data, from 0 to max_size, divided into histo_bucket_count buckets
	long *vcount; //array of vertex distances, calculated with same formula as histogram
	int heap_threshold; //highest bucket where messages can be pushed to heap
	int tram_threshold; //highest bucket where messages can be pushed to tram
	int bfs_threshold;
	double bucket_multiplier; //constant to calculate bucket
	std::vector<Update> *tram_hold; //hold buffer for messages not in tram limit
	std::vector<Update> *pq_hold; //hold for heap messages
	long bfs_created=0; //bfs created messages
	long bfs_processed=0; //bfs processed messages
	int updates_noted=0; //updates that have either updated a vertex value, or are confirmed to not be an improvement
	int *dest_table; //destination table for faster pe calculation
	std::vector<Update> *bfs_hold; // bfs hold (control distance explosion due to bfs)
	int current_phase=0;
	long actual_edges=0; //when graph is generated, here's how many edges actually got generated
	long bfs_noted=0;
	std::vector<Update> local_updates;
	long *info_array;
    long distance_changes=0;
	long updates_in_tram = 0;
	#ifdef PAPI
	int eventset;
	#endif
	int chunk_size = 100000;
	std::vector<Update> *hold_to_process;

public:

	/**
	 * Gets the destination processor for a given vertex 
	*/
	int get_dest_proc(long vertex)
	{
		int dest_proc = 0;
		for (int j = 0; j < N; j++)
		{
			// find first partition that begins at a higher edge count;
			if (vertex >= partition_index[j] && vertex < partition_index[j + 1])
			{
				dest_proc = j;
				break;
			}
			if (j == N - 1)
			{
				dest_proc = N - 1;
			}
		}
		return dest_proc;
	}

	int get_dest_proc_fast(long vertex)
	{
		//look up x/M and 1+x/M
		int xm_pe, xm_plus_one_pe;
		long dest_table_index = vertex / M;
		// if this points to the end of dest_table
		if(dest_table_index >= (V / M) - 1)
		{
			xm_pe = dest_table[(V / M) - 1];
			int dest_proc = xm_pe;
			for (int j = xm_pe; j < N; j++)
			{
				// find first partition that begins at a higher edge count;
				if (vertex >= partition_index[j] && vertex < partition_index[j + 1])
				{
					dest_proc = j;
					break;
				}
				if (j == N - 1)
				{
					dest_proc = N - 1;
				}
			}
			return dest_proc;
		}
		xm_pe = dest_table[dest_table_index];
		xm_plus_one_pe = dest_table[dest_table_index+1];
		int dest_proc = xm_pe;
		for (int j = xm_pe; j <= xm_plus_one_pe; j++)
		{
			// find first partition that begins at a higher edge count;
			if (vertex >= partition_index[j] && vertex < partition_index[j + 1])
			{
				dest_proc = j;
				break;
			}
			if (j == N - 1)
			{
				dest_proc = N - 1;
			}
		}
		return dest_proc;
	}

	int get_dest_proc_local(Update upd) {
		int dest_proc = get_dest_proc_fast(upd.dest_vertex);
/*
		if(dest_proc==CkMyPe()) {
//		local_updates.push_back(upd);
		return -1;
		}
*/
		return dest_proc;
	}

	SsspChares(CProxy_HTram htram)
	{
		tram_proxy = htram;
	}

	void initiate_pointers()
	{
		tram = tram_proxy.ckLocalBranch();
		tram->set_func_ptr_retarr(SsspChares::process_update_caller, get_dest_proc_local_caller, done_caller, this);
		shared_local = shared.ckLocalBranch();
		#ifdef PAPI
		eventset = PAPI_NULL;
		int result = PAPI_create_eventset(&eventset);
		if (result!=PAPI_OK)
		{
			printf("Error PAPI create eventset: %s\n",PAPI_strerror(result));
		}
		result=PAPI_add_event(eventset,PAPI_TOT_INS);
    	if (result!=PAPI_OK) {
      		printf("Error PAPI add_event %s\n",
       		PAPI_strerror(result));
    	}
		#endif
	}

	bool idle_triggered(double)
	{
		process_heap();
		return true;
	}

	void initialize_data(long *partition, int dividers)
	{
		histogram = new long[histo_bucket_count];
		vcount = new long[histo_bucket_count+1]; //histo buckets plus infty
		for(int i=0; i<histo_bucket_count; i++)
		{
			histogram[i] = 0;
			vcount[i] = 0;
		}
		vcount[histo_bucket_count] = 0;
		partition_index = new long[dividers];
		for (int i = 0; i < dividers; i++)
		{
			partition_index[i] = partition[i];
		}
		start_vertex = partition_index[thisIndex];
		num_vertices = partition_index[CkMyPe()+1] - partition_index[CkMyPe()];
		dest_table = new int[V / M];
		for(int i=0, j=0; i < V; j++, i=j*M)
		{
			dest_table[j] = get_dest_proc(i);
		}
		local_graph = new Node[num_vertices];
		heap_threshold = initial_threshold;
		tram_threshold = initial_threshold + 2;
		bfs_threshold = heap_threshold;
		tram_hold = new std::vector<Update>[histo_bucket_count];
		pq_hold = new std::vector<Update>[histo_bucket_count];
		hold_to_process = new std::vector<Update>[histo_bucket_count];
		for(int i=0; i<histo_bucket_count; i++)
		{
			tram_hold[i].reserve(4096);
			pq_hold[i].reserve(4096);
			hold_to_process[i].reserve(4096);
		}
		bfs_hold = new std::vector<Update>[histo_bucket_count];
		info_array = new long[histo_reduction_width+7];
		bucket_multiplier = histo_bucket_count / (histo_bucket_count * log(V));
		CkCallWhenIdle(CkIndex_SsspChares::idle_triggered(0), this);
	}

	void generate_2d_graph(long *partition, int dividers)
	{
		initialize_data(partition, dividers);
		bucket_multiplier = histo_bucket_count / (histo_bucket_count * sqrt(V));
		#ifdef INFO_PRINTS
		ckout << "Generating local graph on PE " << CkMyPe() << " with " << num_vertices << " vertices" << endl;
		#endif
		cost *largest_outedges = new cost[num_vertices];
		long side_length = (int) std::sqrt((double) V);
		bucket_multiplier = histo_bucket_count / (histo_bucket_count * sqrt(V));
		for(int i=0; i<num_vertices; i++)
		{
			Node new_node;
			new_node.home_process = thisIndex;
			new_node.distance = lmax;
			std::vector<Edge> adj;
			new_node.adjacent = adj;
			vcount[histo_bucket_count]++;
			long largest_outedge = 0;
			long this_vertex = (long) i + start_vertex;
			std::mt19937 generator(this_vertex + S);
			std::uniform_int_distribution<cost> edge_weight_distribution(1,1000);
			long x_index = this_vertex / side_length;
			long y_index = this_vertex % side_length;
			for(int j=-1; j<=1; j+=2)
			{
				long neighbor_x = x_index + j;
				long neighbor_y = y_index;
				if((neighbor_x >= 0) && (neighbor_y >= 0) && (neighbor_x < side_length) && (neighbor_y < side_length))
				{
					actual_edges++;
					Edge new_edge;
					new_edge.end = neighbor_x * side_length + neighbor_y;
					new_edge.distance = edge_weight_distribution(generator);
					if(new_edge.distance > largest_outedge) largest_outedge = new_edge.distance;
					new_node.adjacent.push_back(new_edge);
				}
			}
			for(int j=-1; j<=1; j+=2)
			{
				long neighbor_x = x_index;
				long neighbor_y = y_index + j;
				if((neighbor_x >= 0) && (neighbor_y >= 0) && (neighbor_x < side_length) && (neighbor_y < side_length))
				{
					actual_edges++;
					Edge new_edge;
					new_edge.end = neighbor_x * side_length + neighbor_y;
					new_edge.distance = edge_weight_distribution(generator);
					if(new_edge.distance > largest_outedge) largest_outedge = new_edge.distance;
					new_node.adjacent.push_back(new_edge);
				}
			}
			if((x_index==0 && y_index == 0) || (x_index==side_length-1 && y_index == 0) || (x_index==0 && y_index == side_length-1) || (x_index==side_length-1 && y_index == side_length - 1))
			{
				if(new_node.adjacent.size() != 2) ckout << "Edge count wrong for vertex " << this_vertex << " should be 2 not " << new_node.adjacent.size() << endl;
			}
			else if(x_index==0 || y_index == 0 || x_index == side_length - 1 || y_index == side_length - 1)
			{
				if(new_node.adjacent.size() != 3) ckout << "Edge count wrong for vertex " << this_vertex << " should be 3 not " << new_node.adjacent.size() << endl;
			}
			else
			{
				if(new_node.adjacent.size() != 4) ckout << "Edge count wrong for vertex " << this_vertex << " should be 4 not " << new_node.adjacent.size() << endl;
			}
			std::sort(new_node.adjacent.begin(), new_node.adjacent.end(), [](Edge a, Edge b){return a.distance < b.distance;});
			local_graph[i] = new_node;
			largest_outedges[i] = largest_outedge;
		}
		cost max_edges_sum = 0;
		for(int i=0; i<num_vertices; i++)
		{
			max_edges_sum += largest_outedges[i];
		}
		#ifdef INFO_PRINTS
		ckout << "PE " << CkMyPe() << " generated " << actual_edges << " edges" << endl;
		#endif
		CkCallback cb(CkReductionTarget(Main, begin), mainProxy);
		contribute(sizeof(cost), &max_edges_sum, CkReduction::sum_long, cb);
	}

	void generate_local_graph(long _num_vertices, long _num_edges, long *partition, int dividers)
	{
		#ifdef INFO_PRINTS
		ckout << "Generating local graph on PE " << CkMyPe() << " with " << _num_vertices << " vertices and " << _num_edges << " edges" << endl;
		#endif
		initialize_data(partition, dividers);
		cost *largest_outedges = new cost[num_vertices];
		for(int i=0; i<num_vertices; i++)
		{
			Node new_node;
			new_node.home_process = thisIndex;
			new_node.distance = lmax;
			std::vector<Edge> adj;
			new_node.adjacent = adj;
			vcount[histo_bucket_count]++;
			long largest_outedge = 0;
			std::mt19937 generator( (long) i + start_vertex);
			std::uniform_int_distribution<long> edge_count_distribution(0,2*average_degree);
			std::uniform_int_distribution<long> edge_dest_distribution(0,V-1);
			std::uniform_int_distribution<cost> edge_weight_distribution(1,1000);
			long num_edges = edge_count_distribution(generator);
			long *edge_destinations = new long[num_edges];
			for(int j=0; j<num_edges; j++)
			{
				edge_destinations[j] = -1;
			}
			if((CkMyPe()==N-1) && (i >= _num_vertices)) continue; 
			for(int j=0; j<num_edges; j++)
			{
				actual_edges++;
				Edge new_edge;
				bool repeated = true;
				long candidate_end = edge_dest_distribution(generator);
				//logic to keep destinations different
				while(repeated)
				{
					bool different = true;
					for(int k=0; k<j; k++)
					{
						if(edge_destinations[k] == candidate_end)
						{
							different = false;
							break;
						}
					}
					if(different)
					{
						new_edge.end = candidate_end;
						repeated = false;
						edge_destinations[j] = candidate_end;
					}
					else candidate_end = edge_dest_distribution(generator);
				}
				new_edge.distance = edge_weight_distribution(generator);
				if(new_edge.distance > largest_outedge) largest_outedge = new_edge.distance;
				new_node.adjacent.push_back(new_edge);
			}
			std::sort(new_node.adjacent.begin(), new_node.adjacent.end(), [](Edge a, Edge b){return a.distance < b.distance;});
			local_graph[i] = new_node;
			largest_outedges[i] = largest_outedge;
		}
		cost max_edges_sum = 0;
		for(int i=0; i<num_vertices; i++)
		{
			max_edges_sum += largest_outedges[i];
		}
		CkCallback cb(CkReductionTarget(Main, begin), mainProxy);
		contribute(sizeof(cost), &max_edges_sum, CkReduction::sum_long, cb);
	}

	void get_graph(LongEdge *edges, long E, long *partition, int dividers)
	{
		actual_edges = E;
		initialize_data(partition, dividers);
		cost *largest_outedges = new cost[num_vertices];
		if (num_vertices != 0)
		{
			for (int i = 0; i < num_vertices; i++)
			{
				Node new_node;
				new_node.home_process = thisIndex;
				new_node.distance = lmax;
				std::vector<Edge> adj;
				new_node.adjacent = adj;
				local_graph[i] = new_node;
				vcount[histo_bucket_count]++;
				largest_outedges[i] = 0;
			}
			for (int i = 0; i < E; i++)
			{
				Edge new_edge;
				new_edge.end = edges[i].end;
				new_edge.distance = edges[i].distance;
				int new_edge_origin = edges[i].begin - start_vertex;
				local_graph[new_edge_origin].adjacent.push_back(new_edge);
				if(edges[i].distance > largest_outedges[new_edge_origin]) largest_outedges[new_edge_origin] = edges[i].distance;
			}
			for(int i=0; i < num_vertices; i++)
			{
				std::sort(local_graph[i].adjacent.begin(), local_graph[i].adjacent.end(), [](Edge a, Edge b){return a.distance < b.distance;});
			}
		}
		//reduce largest edge
		cost max_edges_sum = 0;
		if (num_vertices != 0)
		{
			for(int i=0; i<num_vertices; i++)
			{
				max_edges_sum += largest_outedges[i];
			}
		}
		CkCallback cb(CkReductionTarget(Main, begin), mainProxy);
		contribute(sizeof(cost), &max_edges_sum, CkReduction::sum_long, cb);
	}

	void start_papi()
	{
		#ifdef PAPI
		int result = PAPI_start(eventset);
		if (result != PAPI_OK)
		{
			printf("Error PAPI start: %s\n",PAPI_strerror(result));
		}
		#endif
		traceBegin();
	}

	/**
	 * Method that accepts initial update to source vertex
	*/
	void start_algo(Update new_vertex_and_distance)
	{
		process_update(new_vertex_and_distance);
		//arr[thisIndex].process_heap();
	}

	static void process_update_caller(void *p, Update *new_vertex_and_distances, int count)
	{
		//ckout << "PE " << CkMyPe() << " receiving " << count << " updates" << endl;
		for(int i=0; i<count; i++)
		{
			((SsspChares *)p)->process_update(new_vertex_and_distances[i]);
		}
		//((SsspChares *)p)->process_heap();

	}

	static int get_dest_proc_local_caller(void *p, Update new_upd)
	{
		return ((SsspChares *)p)->get_dest_proc_local(new_upd);
	}

	static void done_caller(void *p) 
	{
		((SsspChares *)p)->process_local_updates();
	}

	/**
	 * Gets the histogram bucket for any given distance
	*/
	int get_histo_bucket(cost distance)
	{
		double bucket = distance * bucket_multiplier;
		int result = (int) bucket;
		if(result >= histo_bucket_count) return histo_bucket_count - 1;
		else return result;
	}

	void generate_updates(long local_index, bool bfs)
	{
		for (int i = 0; i < local_graph[local_index].adjacent.size(); i++)
		{ 
			// calculate distance pair for neighbor
			Update new_update;
			new_update.dest_vertex = local_graph[local_index].adjacent[i].end;
			new_update.distance = local_graph[local_index].distance + local_graph[local_index].adjacent[i].distance;
			//we are going to send this, so add to the histogram and the send update count
			int neighbor_bucket = get_histo_bucket(new_update.distance);
			histogram[neighbor_bucket]++;
			updates_created_locally++;
			//if exceeds limit, put in hold
			#ifndef ALL_TO_TRAM_HOLD
			if((neighbor_bucket > tram_threshold) && !bfs)
			{
				tram->sendItemPrioDeferredDest(new_update, neighbor_bucket);//tram_hold[neighbor_bucket].push_back(new_update);
			}
			else
			{
				//calculated dest proc
				int dest_proc = get_dest_proc_fast(new_update.dest_vertex);
				#ifndef LOCAL_TO_TRAM
				if(dest_proc==CkMyPe())
				{
					process_update(new_update);
				}
				else
				{
					tram->sendItemPrioDeferredDest(new_update, 0);//tram->insertValue(new_update, dest_proc);
				} 
				#else
				tram->sendItemPrioDeferredDest(new_update, 0);//tram->insertValue(new_update, dest_proc);//this gets called
				#endif
			}
			#else
      		tram->sendItemPrioDeferredDest(new_update, neighbor_bucket);
//			tram_hold[neighbor_bucket].push_back(new_update);
			if(neighbor_bucket <= tram_threshold) updates_in_tram++;
#if 0
			if(updates_in_tram == 8192)
			{
      			tram->insertBuckets(tram_threshold);
      			updates_in_tram = 0;
    		}
#endif
			#endif
		}
	}

	void generateUpdatesOtherPe(long local_index, Node *local_graph_other)
	{
		for (int i = 0; i < local_graph_other[local_index].adjacent.size(); i++)
		{ 
			// calculate distance pair for neighbor
			Update new_update;
			new_update.dest_vertex = local_graph_other[local_index].adjacent[i].end;
			new_update.distance = local_graph_other[local_index].distance + local_graph_other[local_index].adjacent[i].distance;
			//we are going to send this, so add to the histogram and the send update count
			int neighbor_bucket = get_histo_bucket(new_update.distance);
			histogram[neighbor_bucket]++;
			updates_created_locally++;
			//if exceeds limit, put in hold
			#ifndef ALL_TO_TRAM_HOLD
			if(neighbor_bucket > tram_threshold)
			{
				tram->sendItemPrioDeferredDest(new_update, neighbor_bucket);//tram_hold[neighbor_bucket].push_back(new_update);
			}
			else
			{
				//calculated dest proc
				int dest_proc = get_dest_proc_fast(new_update.dest_vertex);
				#ifndef LOCAL_TO_TRAM
				if(dest_proc==CkMyPe())
				{
					process_update(new_update);
				}
				else
				{
					tram->sendItemPrioDeferredDest(new_update, 0);//tram->insertValue(new_update, dest_proc);
				} 
				#else
				tram->sendItemPrioDeferredDest(new_update, 0);//tram->insertValue(new_update, dest_proc);//this gets called
				#endif
			}
			#else
      		tram->sendItemPrioDeferredDest(new_update, neighbor_bucket);
			if(neighbor_bucket <= tram_threshold) updates_in_tram++;
			#endif
		}

	}

	void processHeapChunk(int first, int last, int send_chare, int bucket, std::vector<Update> *whole_bucket, long data_count, Node *local_graph_other, long vertex_count, long start_vertex)
	{
		Update *data = whole_bucket -> data();
		for(int i=first; i<last; i++)
		{
			Update new_vertex_and_distance = data[i];
			long dest_vertex = new_vertex_and_distance.dest_vertex;
			cost new_distance = new_vertex_and_distance.distance;
			int this_histo_bucket = get_histo_bucket(new_distance);
			long local_index = dest_vertex - start_vertex;
			if (new_distance == local_graph_other[local_index].distance)
			{
					// for all neighbors
					generateUpdatesOtherPe(local_index, local_graph_other);
			}
			else
			{
				rejected_updates++;
			}
			wasted_updates++;
			histogram[this_histo_bucket]--;
			updates_processed_locally++;
		}
		shared_local -> chunks_remaining[send_chare][bucket]--;
		if(shared_local -> chunks_remaining[send_chare][bucket]==0)
		{
			whole_bucket -> clear();
		}
	}

	/**
	 * Takes a distance update and immediately adds it to the local heap/pq
	 */
	inline void process_update(Update new_vertex_and_distance)
	{
		long dest_vertex = new_vertex_and_distance.dest_vertex;
		long local_index = dest_vertex - start_vertex;
		cost this_cost = new_vertex_and_distance.distance;
		int this_bucket = get_histo_bucket(this_cost);
		if (this_cost < local_graph[local_index].distance)
		{
			/*
			vcount[this_bucket]++;
			updates_noted++;
			if(local_graph[local_index].distance == lmax)
			{ 
				bfs_noted++;
				local_graph[local_index].distance = this_cost;
				vcount[histo_bucket_count]--;
				if (this_bucket > bfs_threshold)
				{
					local_graph[local_index].send_updates = true;
			      	bfs_hold[this_bucket].push_back(new_vertex_and_distance);
				}
				else
				{
					bfs_processed++;
					generate_updates(local_index, true);
					wasted_updates++;
					histogram[this_bucket]--;
					updates_processed_locally++;
				}
			}
			else
			{ 
				vcount[get_histo_bucket(local_graph[local_index].distance)]--;
				local_graph[local_index].distance = this_cost;
				local_graph[local_index].send_updates = true;
				if(this_bucket > heap_threshold)
				{
					pq_hold[this_bucket].push_back(new_vertex_and_distance);
				}
				else pq.push(new_vertex_and_distance);
			}
			*/
			#ifdef VCOUNT
			vcount[this_bucket]++;
			if(local_graph[local_index].distance == lmax)
			{
				vcount[histo_bucket_count]--;
			}
			else vcount[get_histo_bucket(local_graph[local_index].distance)]--;
			#endif
			local_graph[local_index].distance = this_cost;
			distance_changes++;
			updates_noted++;
			int pq_bucket;
			if(local_graph[local_index].adjacent.size()>0)
			{
				#ifdef PQ_EDGE_DIST
				pq_bucket = get_histo_bucket(local_graph[local_index].adjacent[0].distance + this_cost);
				#ifndef PQ_HOLD_ONLY
				if(pq_bucket > heap_threshold)
				{
					pq_hold[pq_bucket].push_back(new_vertex_and_distance);
				}
				else pq.push(new_vertex_and_distance);
				#else
				pq_hold[pq_bucket].push_back(new_vertex_and_distance);
				#endif
				#else
				#ifndef PQ_HOLD_ONLY
				if(this_bucket > heap_threshold)
				{
					pq_hold[this_bucket].push_back(new_vertex_and_distance);
				}
				else pq.push(new_vertex_and_distance);
				#else
				pq_hold[this_bucket].push_back(new_vertex_and_distance);
				#endif
				#endif
			}
			else
			{
				wasted_updates++;
				histogram[this_bucket]--;
				updates_processed_locally++;
			}
		}
		else
		{
			wasted_updates++;
			rejected_updates++;
			histogram[this_bucket]--;
			updates_noted++;
			updates_processed_locally++;
		}
		
	}

	/**
	 * Update distances, but locally (the incoming pair comes from this PE)
	 * this is not an entry method
	 * returns true (runs when pe is idle)
	 */
	void process_heap()
	{
		#ifdef PQ_HOLD_ONLY
		for(int i=0; i<=heap_threshold; i++)//iterate to heap threshold
		{
			long items_processed = 0;
			#ifdef NODE_LOAD_BALANCE
			if(pq_hold[i].size()>chunk_size)
			{
				//load all updates into hold_to_process
				traceBeginUserBracketEvent(shared_local -> bracketed_id);
				for(int j=0; j<pq_hold[i].size(); j++)
				{
					hold_to_process[i].push_back(pq_hold[i][j]);
				}
				traceEndUserBracketEvent(shared_local -> bracketed_id);
				//make separate vectors
				int chunk_count = (pq_hold[i].size() / chunk_size) + 1;
				if(shared_local -> chunks_remaining[thisIndex][i] != 0) continue;
				else shared_local -> chunks_remaining[thisIndex][i] = chunk_count;
				for(int j=0; j<pq_hold[i].size(); j+=chunk_size)
				{
					//fill array
					int low = j;
					int high = j + chunk_size;
					if(high>pq_hold[i].size()) high = pq_hold[i].size();
					traceBeginUserBracketEvent(shared_local -> othercaller_id);
					process_heap_shared[CkMyNode()].processHeapOtherCaller(low, high, thisIndex, i, (intptr_t) &hold_to_process[i], high-low, (intptr_t) local_graph, num_vertices, start_vertex);
					traceEndUserBracketEvent(shared_local -> othercaller_id);
				}
				pq_hold[i].clear();
				//arr[thisIndex].process_heap();
				//break;
			}
			else
			{
				items_processed = pq_hold[i].size();
				for(int j=0; j<pq_hold[i].size(); j++)//iterate pq bucket in reverse
				{
					Update new_vertex_and_distance = pq_hold[i][j];
					long dest_vertex = new_vertex_and_distance.dest_vertex;
					cost new_distance = new_vertex_and_distance.distance;
					int this_histo_bucket = get_histo_bucket(new_distance);
					long local_index = dest_vertex - start_vertex;
					if (new_distance == local_graph[local_index].distance)
					{
							// for all neighbors
							generate_updates(local_index, false);
					}
					else
					{
						rejected_updates++;
					}
					wasted_updates++;
					histogram[this_histo_bucket]--;
					updates_processed_locally++;
				}
				if(items_processed>0)
				{
					pq_hold[i].clear();
					arr[thisIndex].process_heap();
					break;
				}
			}
			#else
			items_processed = pq_hold[i].size();
			for(int j=0; j<pq_hold[i].size(); j++)//iterate pq bucket in reverse
			{
				Update new_vertex_and_distance = pq_hold[i][j];
				long dest_vertex = new_vertex_and_distance.dest_vertex;
				cost new_distance = new_vertex_and_distance.distance;
				int this_histo_bucket = get_histo_bucket(new_distance);
				long local_index = dest_vertex - start_vertex;
				if (new_distance == local_graph[local_index].distance)
				{
						// for all neighbors
						generate_updates(local_index, false);
				}
				else
				{
					rejected_updates++;
				}
				wasted_updates++;
				histogram[this_histo_bucket]--;
				updates_processed_locally++;
			}
			if(items_processed>0)
			{
				pq_hold[i].clear();
				arr[thisIndex].process_heap();
				break;
			}
			#endif
		}
		#else
		int heap_count = 0;
		while (pq.size() > 0)
		{
			if (++heap_count > 100) { thisProxy[thisIndex].process_heap(); break; } // give other eps a chance to run  
			Update new_vertex_and_distance = pq.top();
			long dest_vertex = new_vertex_and_distance.dest_vertex;
			cost new_distance = new_vertex_and_distance.distance;
			int this_histo_bucket = get_histo_bucket(new_distance);
			if(this_histo_bucket > heap_threshold)
			{
				break;
			}
			pq.pop();
			if (dest_vertex >= partition_index[thisIndex] && dest_vertex < partition_index[thisIndex + 1])
			{
				long local_index = dest_vertex - start_vertex;
				//  if the incoming distance is actually smaller
				if (new_distance == local_graph[local_index].distance)
				{
						generate_updates(local_index, false);
				}
				else
				{
					rejected_updates++;
				}
			}
			wasted_updates++;
			histogram[this_histo_bucket]--;
			updates_processed_locally++;
		}
		#endif
	}


	/**
	 * Checks if anything is in the buffer (false quiescence)
	 */
	void check_buffer()
	{
		// ckout << "Checking message stats" << endl;
		/*
		int msg_stats[2];
		msg_stats[0] = updates_created_locally;
		msg_stats[1] = updates_processed_locally;
		CkCallback cb(CkReductionTarget(Main, check_buffer_done), mainProxy);
		contribute(2 * sizeof(int), msg_stats, CkReduction::sum_int, cb);
		*/
	}

	/**
	 * Contribute to a reduction to get the overall histogram to pe 0/main chare
	*/
	void contribute_histogram(int behind_first_nonzero)
	{
		long donecount = 0;
		traceUserEvent(shared_local -> event_id);
		CkCallback cb(CkReductionTarget(Main, reduce_histogram), mainProxy);
		int first_nonzero = behind_first_nonzero + 1;
		for(int i = first_nonzero; i < (first_nonzero + histo_reduction_width); i++)
		{
			if(i >= histo_bucket_count) info_array[i-first_nonzero] = 0;
			else info_array[i-first_nonzero] = histogram[i];
		}
		info_array[histo_reduction_width] = updates_created_locally;
		info_array[histo_reduction_width+1] = updates_processed_locally;
		info_array[histo_reduction_width+2] = bfs_processed;
		for (int i=0; i<=behind_first_nonzero;i++) donecount += vcount[i];
		info_array[histo_reduction_width+3] = donecount;
		info_array[histo_reduction_width+4] = updates_noted;
		info_array[histo_reduction_width+5] = bfs_noted;
		info_array[histo_reduction_width+6] = distance_changes;
		contribute((histo_reduction_width+7) * sizeof(long), info_array, CkReduction::sum_long, cb);
	}

	/**
	 * Called when some of the buffers aren't full, meaning we need to keep the algorithm going
	 */
	void keep_going()
	{
		/*
		// everything in the tram hold gets added to tram
		for(int i=0; i<histo_bucket_count; i++)
		{
			for(int j=0; j<tram_hold[i].size(); j++)
			{
				int dest_proc = get_dest_proc_fast(tram_hold[i][j].dest_vertex);
				if(dest_proc==CkMyPe())
				{
					pq.push(tram_hold[i][j]);
				}
				else tram->insertValue(tram_hold[i][j], dest_proc);
			}
			tram_hold[i].clear();
		}
		for(int i=0; i<histo_bucket_count; i++)
		{
			for(int j=0; j<pq_hold[i].size(); j++)
			{
				pq.push(pq_hold[i][j]);
			}
			pq_hold[i].clear();
		}
		tram->tflush();
		//arr[thisIndex].process_heap();
		*/
	}

	void clear_pq_hold()
	{
		//we should maintain lower bound
		for(int i=0; i<=heap_threshold; i++)
		{
			for(int j=0; j<pq_hold[i].size(); j++)
			{
				pq.push(pq_hold[i][j]);
			}
			pq_hold[i].clear();
		}
	}

	void process_local_updates()
	{
		for(int i=0; i<local_updates.size(); i++)
		{
			process_update(local_updates[i]);
		}
		local_updates.clear();
	}

	/**
	 * Broadcasts bucket limit
	*/
	void current_thresholds(int _heap_threshold, int _tram_threshold, int _bfs_threshold, int behind_first_nonzero, int phase)
	{
		heap_threshold = _heap_threshold;
		tram_threshold = _tram_threshold;
		bfs_threshold = _bfs_threshold;
		current_phase = phase;
		//after every reduction, push out messages in hold that are in limit
		//replace this loop with call to tram->changethreshold(tram_threshold)
		tram->shareArrayOfBuckets(tram_hold, histo_bucket_count);
    int direct_threshold = behind_first_nonzero + 8 ;
	//int direct_threshold = tram_threshold;
    if(direct_threshold > tram_threshold-1) direct_threshold = tram_threshold-1;
    float selectivity = 1.0;
    //if(behind_first_nonzero > 68) selectivity = 1.0;
    tram->changeThreshold(direct_threshold, tram_threshold, selectivity);
		#if 0
		for(int i=0; i<=tram_threshold; i++)
		{
			for(int j=0; j<tram_hold[i].size(); j++)
			{
				int dest_proc = get_dest_proc_fast(tram_hold[i][j].dest_vertex);
				if(dest_proc==CkMyPe())
				{
					//process_update(tram_hold[i][j]); //todo: put it in a vector, then loop over it and call process_update
					local_updates.push_back(tram_hold[i][j]);
				}
				else tram->insertValue(tram_hold[i][j], dest_proc); 
			}
			tram_hold[i].clear();
		}
		#endif
		#ifndef PQ_HOLD_ONLY
		arr[thisIndex].clear_pq_hold();
		//add user event
		#endif
		process_local_updates();
		/*
		for(int i=0; i<=bfs_threshold; i++)
		{
			for(int j=0; j<bfs_hold[i].size(); j++)
			{ 
				long local_index = bfs_hold[i][j].dest_vertex - start_vertex;
				cost this_cost = bfs_hold[i][j].distance;
				int this_bucket = get_histo_bucket(this_cost);
				if (this_cost == local_graph[local_index].distance)
				{
					generate_updates(local_index, true);
				}
				else if (this_cost > local_graph[local_index].distance)
				{
					rejected_updates++;
				}
				bfs_processed++;
				wasted_updates++;
				histogram[this_bucket]--;
				updates_processed_locally++;
			}
			bfs_hold[i].clear();
		}
		*/
		if(rand()%5==0) tram->tflush();
//    tram->sanityCheck();
//    tram->flush_everything();
		arr[thisIndex].process_heap();
		contribute_histogram(behind_first_nonzero);
	}

	/**
	 * Print out the final distances calculated by the algorithm
	 */
	void print_distances()
	{
		/*
		 //enable only for smaller graphs
		for (int i = 0; i < num_vertices; i++)
		{
			ckout << "Partition " << thisIndex << " vertex num " << local_graph[i] << " distance " << local_graph[i].distance << endl;
		}
		*/
		traceEnd();
		#ifdef PAPI
		long long values[1] = {(long long) 0};
		int result = PAPI_stop(eventset,values);
		if(result != PAPI_OK)
		{
			printf("Error PAPI stop %s\n", PAPI_strerror(result));
		}
		//ckout << "PE " << CkMyPe() << " total instructions: " << values[0] << endl;
		long msg_stats[7+histo_bucket_count];
		#else
		long msg_stats[6+histo_bucket_count];
		#endif
		msg_stats[0] = wasted_updates;
		msg_stats[1] = rejected_updates;
		for(int i=0; i<histo_bucket_count + 1; i++)
		{
			msg_stats[i+2] = vcount[i];
		}
		msg_stats[3+histo_bucket_count] = updates_noted;
		msg_stats[4+histo_bucket_count] = actual_edges;
		msg_stats[5+histo_bucket_count] = distance_changes;
		#ifdef PAPI
		msg_stats[6+histo_bucket_count] = values[0];
		#endif

		CkCallback cb(CkReductionTarget(Main, done), mainProxy);

		#ifdef PAPI
		contribute((7+histo_bucket_count) * sizeof(long), msg_stats, CkReduction::sum_long, cb);
		#else
		contribute((6+histo_bucket_count) * sizeof(long), msg_stats, CkReduction::sum_long, cb);
		#endif
		// mainProxy.done();
	}

	void get_max_cost()
	{
		cost max_cost = 0;
		for (int i = 0; i < num_vertices; i++)
		{
			cost vertex_cost = local_graph[i].distance;
			if(vertex_cost != lmax)
			{
				if(vertex_cost > max_cost)
				{
					max_cost = vertex_cost;
				}
			}
		}
		CkCallback cb(CkReductionTarget(Main, done_max_cost), mainProxy);
		contribute(sizeof(cost), &max_cost, CkReduction::max_long, cb);
	}

	void stop_periodic_flush()
	{
		tram->stop_periodic_flush();
	}
};

class processHeapShared : public CBase_processHeapShared
{
	public:
	processHeapShared(){}

	void processHeapOtherCaller(int low, int high, int send_chare, int bucket, intptr_t whole_bucket_pointer, long data_count, intptr_t local_graph_pointer, long vertex_count, long startVertex)
	{
		SsspChares *c = arr[CkMyPe()].ckLocal();
		Node *graph = (Node *) local_graph_pointer;
		std::vector<Update> *whole_bucket = (std::vector<Update> *) whole_bucket_pointer;
		c -> processHeapChunk(low, high, send_chare, bucket, whole_bucket, data_count, graph, vertex_count, startVertex);
	}
};

#include "weighted_htram_smp.def.h"
