CXX = g++ 
CXXFLAGS = -O2 -Wall -std=c++11

all: clusterMC_ising_thermo_prop 

clusterMC_ising_thermo_prop: clusterMC_ising_thermo_prop.o nodes.o
	$(CXX) -o clusterMC_ising_thermo_prop clusterMC_ising_thermo_prop.o \
	    nodes.o

check: test-random-selector test-random-selector2

test-random-selector: test-random-selector.o
	$(CXX) -o test-random-selector test-random-selector.o

test-random-selector2: test-random-selector2.o
	$(CXX) -o test-random-selector2 test-random-selector2.o

clean:
	rm -f clusterMC_ising_thermo_prop test-random-selector \
	   test-random-selector2  *.o
