CXX = g++ 
CXXFLAGS = -O2 -Wall -std=c++11

all: clusterMC_ising_thermo_prop 

clusterMC_ising_thermo_prop: clusterMC_ising_thermo_prop.o nodes.o
	$(CXX) -o clusterMC_ising_thermo_prop clusterMC_ising_thermo_prop.o \
	    nodes.o

clean:
	rm -f clusterMC_ising_thermo_prop *.o
