#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <iterator>
#include <random>

#include "nodes.h"
#include "random_selector.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

size_t nsample(10000);

////-----------  Templates to select random elements from a container  ---------//
//// Taken from http://stackoverflow.com/a/16421677
//
//template<typename Iter, typename RandomGenerator>
//Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
//    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
//    std::advance(start, dis(g));
//    return start;
//}
//
//template<typename Iter>
//Iter select_randomly(Iter start, Iter end) {
//    static std::random_device rd;
//    static std::mt19937 gen(rd());
//    return select_randomly(start, end, gen);
//}

//----------------------------------------------------------------------------//

void do_ising(ising::nodes& lattice, double* spin, const double T)
{
    double av_m(0);
    double E_av(0);
    double E2_av(0);
    double beta = 1.0/T;

    // Cluster expansion probability for a 2D ising model with all sites filled
    const double p = 1.0 - std::exp(-2.0 * beta);


    double Cv = beta*beta/lattice.nsites*(E2_av - E_av*E_av);
    std::cout << std::setw(15) << T
	      << std::setw(15) << av_m
	      << std::setw(15) << Cv
	      << std::setw(15) << nsample << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    if (argc < 1) {
        std::cerr << "usage: clusterMC_ising_thermo_prop > a.dat" << std::endl;
        return EXIT_FAILURE;
    }

    // Define lattice size and initialize it
    const int nx = 12;
    const int ny = 12;
    const int nz = 1;

    ising::nodes lattice;
    lattice.init(nx, ny, nz);

    // Set up the initial spins
    double initial_spin[lattice.nsites];

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,1);
    for(size_t i = 0; i < lattice.nsites; ++i){
	if(distribution(generator) == 0)
	    initial_spin[i] = -1;
	else
	    initial_spin[i] = 1;
    }

    // Define the initial temperature and increments
    double T0 = 0.25;
    double dT = 0.25;

    for(size_t i = 0; i < 12; ++i){
	double T = T0 + dT*i;

	double spin[lattice.nsites];
	std::copy(initial_spin, initial_spin + lattice.nsites, spin);

	do_ising(lattice, spin, T);
    }
}
