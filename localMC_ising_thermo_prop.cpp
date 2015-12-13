#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <random>

#include "nodes.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

size_t nsample(10000);

//----------------------------------------------------------------------------//

void do_ising(ising::nodes& lattice, int* spin, const double T,
	      std::default_random_engine& generator)
{
    double M_av(0);
    double beta = 1.0/T;

    // Set up random numbers
    std::uniform_int_distribution<int> rand_lattice_site(0,lattice.nsites - 1);
    std::uniform_real_distribution<double> rand_01(0.0, 1.0);

    // The sampling is determined by the total number of spins that have been
    // flipped. Required number of spins grows with number of lattice sites
    int n_spins_to_flip(lattice.nsites * nsample);

    for(int n = 0; n < n_spins_to_flip; ++n){

	int active_site = rand_lattice_site(generator);

        // Now loop over the neighbors (nbr) of the active_site
        int E0 = (0);
        for(size_t ibr = 0; ibr < lattice.neighbors[active_site].size(); ++ibr){
            int nbr = lattice.neighbors[active_site][ibr];

            E0 -= spin[active_site] * spin[nbr];
        }

        // Change in energy for ising model is (E' - E0) == -2.0*E0
        double dE = -2.0*E0;

        // Monte-Carlo acceptance criteria
        double p = std::min(1.0, std::exp(-beta*dE));
        if( rand_01(generator) < p){
            spin[active_site] *= -1;
        }

	int M(0);
	for(int i = 0; i < lattice.nsites; ++i)
	    M += spin[i];

	M_av += ((double)std::abs(M));

#if 0
	std::cout << " <<<<<<< step: " << n << ", "
                  << M/((double)lattice.nsites)
                  << " = "
                  << ((double)M_av)/((double)n + 1)/((double)lattice.nsites)
                  << " >>>>>>>" << std::endl;
	for(size_t i = 0; i < lattice.nx; ++i){
	    for(size_t j = 0; j < lattice.ny; ++j){
		if(spin[i*lattice.ny + j] < 0)
		    std::cout << std::setw(2) << '-';
		else
		    std::cout << std::setw(2) << '+';
	    }
	    std::cout << std::endl;
	}
#endif

    }

    std::cout << std::setw(12) << T
	      << std::setw(12) << ((double)M_av)/((double)n_spins_to_flip)
						/((double)lattice.nsites)
	      << std::setw(12) << n_spins_to_flip << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: clusterMC_ising_thermo_prop nx ny nz > a.dat" << std::endl;
        return EXIT_FAILURE;
    }

    // Define lattice size and initialize it
    int nx;
    {
        std::istringstream iss(argv[1]);
        iss >> nx; 
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[1]
                << "' to int" << std::endl;
            return EXIT_FAILURE;
        }
    }
    int ny;
    {
        std::istringstream iss(argv[2]);
        iss >> ny; 
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[2]
                << "' to int" << std::endl;
            return EXIT_FAILURE;
        }
    }
    int nz;
    {
        std::istringstream iss(argv[3]);
        iss >> nz; 
        if (!iss || !iss.eof()) {
            std::cerr << "could not convert '" << argv[3]
                << "' to int" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cerr << "Initializing a " << nx << " x " << ny << " x "
                      << nz << " grid." << std::endl;

    ising::nodes lattice;
    lattice.init(nx, ny, nz);

    // Set up the initial spins
    int initial_spin[lattice.nsites];

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,1);
    for(int i = 0; i < lattice.nsites; ++i){
	if(distribution(generator) == 0)
	    initial_spin[i] = -1;
	else
	    initial_spin[i] = 1;
    }
    int mi = 0;
    for(int i = 0; i < lattice.nsites; ++i)
	mi += initial_spin[i];

    // Define the initial temperature and increments
    double T0 = 1.00;
    double dT = 0.2;
    const size_t nT = 16;

    int spin[lattice.nsites];
    std::copy(initial_spin, initial_spin + lattice.nsites, spin);

    for(size_t i = 0; i < nT; ++i){
	double T = T0 + dT*i;

	do_ising(lattice, spin, T, generator);
    }
}
