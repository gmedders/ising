CXX = g++ 
CXXFLAGS = -O2 -Wall

all: clusterMC_ising_thermo_prop 

clusterMC_ising_thermo_prop: clusterMC_ising_thermo_prop.o nodes.o
	$(CXX) -o clusterMC_ising_thermo_prop clusterMC_ising_thermo_prop.o \
	    nodes.o

nodes.o: nodes.cpp nodes.h
	        $(CXX) -c nodes.cpp $(CXXFLAGS)

clean:
	rm -f clusterMC_ising_thermo_prop *.o
