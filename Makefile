# Change Charm build locations based on where the charm build is
#CHARMC =/ccs/home/rrao/charm/ofi-linux-x86_64-cxi-slurmpmi2cray-gcc/bin/charmc $(OPTS)
CHARMC=charmc $(OPTS)
CHARMC_SMP=charmc $(OPTS)
#CHARMC_SMP =/work/hdd/mzu/rao1/charm/ofi-linux-x86_64-cxi-slurmpmi2cray-smp-gcc/bin/charmc $(OPTS)
# delta: /work/hdd/rao1/charm/ofi-linux-x86_64-cxi-slurmpmi2cray-smp-gcc/bin/charmc
# frontier: /ccs/home/rrao/charm/ofi-linux-x86_64-cxi-slurmpmi2cray-smp-gcc/bin/charmc
# local: /Users/ritvik/charm/netlrts-darwin-x86_64-smp/bin/charmc

CHARMCFLAGS = $(OPTS) -g -O3

BINARY= graph_ckio graph_parallel graph_parallel_ckio weighted_nonsmp weighted_smp weighted_htram_nonsmp weighted_htram_smp weighted_htram_nonsmp_projections weighted_htram_smp_projections weighted_htram_smp_papi
all: $(BINARY)

.PHONY = clean

run_sssp_smp: weighted_htram_smp
	./weighted_htram_smp 10000 160000 100 1 1 +p8 +ppn 1 +setcpuaffinity ++local

graph_serial: graph_serial.cpp
	$(CHARMC) graph_serial.ci
	$(CHARMC) graph_serial.cpp -o $@

graph_ckio: graph_ckio.cpp
	$(CHARMC) graph_ckio.ci
	$(CHARMC) graph_ckio.cpp -o $@ -module CkIO

graph_parallel: graph_parallel.cpp
	$(CHARMC) graph_parallel.ci
	$(CHARMC) graph_parallel.cpp -o $@

graph_parallel_ckio: graph_parallel_ckio.cpp
	$(CHARMC) graph_parallel_ckio.ci
	$(CHARMC) graph_parallel_ckio.cpp -o $@ -module CkIO

weighted_nonsmp: weighted.cpp weighted.ci weighted_node_struct.h
	$(CHARMC) weighted.ci
	$(CHARMC) $(CHARMCFLAGS) $< -o $@ -module NDMeshStreamer

weighted_smp: weighted.cpp weighted.ci weighted_node_struct.h
	$(CHARMC_SMP) weighted.ci
	$(CHARMC_SMP) $(CHARMCFLAGS) $< -o $@ -module NDMeshStreamer

weighted_htram_nonsmp: weighted_htram_nonsmp.cpp weighted_htram_nonsmp.ci weighted_node_struct.h libtramnonsmp.a
	$(CHARMC) weighted_htram_nonsmp.ci -DTRAM_NON_SMP
	$(CHARMC) $(CHARMCFLAGS) libtramnonsmp.a -language charm++ -o $@ $< -std=c++1z -DTRAM_NON_SMP

weighted_htram_nonsmp_projections: weighted_htram_nonsmp.cpp weighted_htram_nonsmp.ci weighted_node_struct.h libtramnonsmp.a
	$(CHARMC) weighted_htram_nonsmp.ci -DTRAM_NON_SMP
	$(CHARMC) $(CHARMCFLAGS) libtramnonsmp.a -language charm++ -o $@ $< -std=c++1z -DTRAM_NON_SMP -tracemode projections

weighted_htram_smp: weighted_htram_smp.cpp weighted_htram_smp.ci weighted_node_struct.h libhtram_group.a
	$(CHARMC_SMP) weighted_htram_smp.ci -DTRAM_SMP -DGROUPBY
	$(CHARMC_SMP) $(CHARMCFLAGS) libhtram_group.a -language charm++ -o $@ $< -std=c++1z -DTRAM_SMP -DGROUPBY

weighted_htram_smp_papi: weighted_htram_smp.cpp weighted_htram_smp.ci weighted_node_struct.h libhtram_group.a
	$(CHARMC_SMP) weighted_htram_smp.ci -DTRAM_SMP -DGROUPBY
	$(CHARMC_SMP) $(CHARMCFLAGS) libhtram_group.a -language charm++ -o $@ $< -std=c++1z -DTRAM_SMP -DGROUPBY -I/opt/cray/pe/papi/7.0.1.2/include -L/opt/cray/pe/papi/7.0.1.2/lib -lpapi

weighted_htram_smp_projections: weighted_htram_smp.cpp weighted_htram_smp.ci weighted_node_struct.h libhtram_group.a
	$(CHARMC_SMP) weighted_htram_smp.ci -DTRAM_SMP -DGROUPBY
	$(CHARMC_SMP) $(CHARMCFLAGS) libhtram_group.a -language charm++ -o $@ $< -std=c++1z -DTRAM_SMP -DGROUPBY -tracemode projections

libtramnonsmp.a : tramNonSmp.o
	$(CHARMC) tramNonSmp.o -o libtramnonsmp.a -language charm++

tramNonSmp.o : tramNonSmp.C tramNonSmp.def.h tramNonSmp.decl.h
	$(CHARMC) -c tramNonSmp.C -o tramNonSmp.o

tramNonSmp.def.h tramNonSmp.decl.h : tramNonSmp.ci
	$(CHARMC) tramNonSmp.ci

libhtram_group.a: htram_group.o
	$(CHARMC_SMP) $< -o $@ -language charm++ -tracemode projections

htram_group.o: htram_group.C htram_group.def.h htram_group.decl.h
	$(CHARMC_SMP) -O3 -c htram_group.C

htram_group.def.h htram_group.decl.h: htram_group.ci
	$(CHARMC_SMP) htram_group.ci

clean:
	$(MAKE) clean-libs
	rm -f *.o *.decl.h *.def.h $(BINARY) charmrun* *.stamp
	rm -rf weighted_htram_Smp_projections

remove-out:
	rm -f *.out

clean-libs:
	rm -f *.def.h *.decl.h
	rm -f *.log.gz *.projrc *.topo *.sts *.sum
	rm -f *.a
