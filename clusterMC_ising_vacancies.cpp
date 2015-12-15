#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <algorithm>
#include <random>

#include "nodes.h"

#define ANYWHERE doit

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

size_t nmove(10000);
//size_t nmove(1000000);

//----------------------------------------------------------------------------//

void swap_spins(int* spin, int& orig, int& dest)
{
    int tmp = spin[dest];
    spin[dest] = spin[orig];
    spin[orig] = tmp;
}

//----------------------------------------------------------------------------//

void print_cell(ising::nodes& lattice, const int* spin, std::string comment)
{    
    std::cerr << comment << std::endl;
    for(int i = 0; i < lattice.nx; ++i){
        for(int j = 0; j < lattice.ny; ++j){
            if(spin[i*lattice.ny + j] < 0)
                std::cerr << std::setw(2) << '-';
            else if(spin[i*lattice.ny + j] == 0)
                std::cerr << std::setw(2) << 'O';
            else
                std::cerr << std::setw(2) << '+';
        }
        std::cerr << std::endl;
    }
}

//----------------------------------------------------------------------------//

int calcE_for_two_connected_sites(ising::nodes& lattice, int* spin,
                                  int& orig_site, int& dest_site)
{
    int E0_orig(0);
    for(size_t ibr = 0; ibr < lattice.neighbors[orig_site].size();++ibr){
        int nbr = lattice.neighbors[orig_site][ibr];
        if(nbr != dest_site){
            E0_orig -= spin[orig_site] * spin[nbr];
        }
    }

    int E0_dest(0);
    for(size_t ibr = 0; ibr < lattice.neighbors[dest_site].size();++ibr){
        int nbr = lattice.neighbors[dest_site][ibr];
        if(nbr != orig_site){
            E0_dest -= spin[dest_site] * spin[nbr];
        }
    }

    return E0_orig + E0_dest;
}

//----------------------------------------------------------------------------//

void do_ising(ising::nodes& lattice, int* spin, const double T,
	      std::default_random_engine& generator)
{
    double M_av(0);
    double beta = 1.0/T;

    // Cluster expansion probability for a 2D ising model with all sites filled
    const double pCluster = 1.0 - std::exp(-2.0 * beta);

    // Set up random numbers
    std::uniform_int_distribution<int> rand_lattice_site(0,lattice.nsites - 1);
    std::uniform_real_distribution<double> rand_01(0.0, 1.0);

    // Randomly zero out the spins for half the sites
    int noccupied(0);
    for(int i = 0; i < lattice.nsites; ++i)
        if(spin[i] != 0)
            ++noccupied;
    if(noccupied > (lattice.nsites / 2)){
        do{

            int delete_this_site = rand_lattice_site(generator);
            spin[delete_this_site] = 0;

            noccupied = 0;
            for(int i = 0; i < lattice.nsites; ++i)
                if(spin[i] != 0)
                    ++noccupied;

        }while(noccupied > (lattice.nsites / 2));
    }

    if(noccupied != lattice.nsites/2){
        std::cerr << "Function to create vacancies not working as expected.\n"
                  << " expected " << lattice.nsites/2 << " occupied sites but"
                  << " found " << noccupied<< std::endl;
        exit(1);
    }

    size_t n_av(0);
    {
        std::string comment("Initial configuration");
        print_cell(lattice, spin, comment);
    }

    // Required number of attempted moves to grow with number of lattice sites
    for(size_t imove = 0; imove < nmove*lattice.nsites; ++imove){

        // For move, select a random site (orig)
        //           and a random neighboring site, including diagonals (dest)
        int orig_site = rand_lattice_site(generator);
#ifdef ANYWHERE
        int dest_site = rand_lattice_site(generator);
        if(orig_site == dest_site)
            continue;
#else
        std::uniform_int_distribution<int>
                   rand_neighbor(0,lattice.diagNeighbors[orig_site].size() - 1);
        int dest_site =
                     lattice.diagNeighbors[orig_site][rand_neighbor(generator)];
#endif

        // Ensure that the move is not trivial
        if(spin[orig_site] != spin[dest_site]){

            // Calculate the energy before the move
            const int E0 = calcE_for_two_connected_sites(lattice, spin,
                                                         orig_site, dest_site);

            // Create trial move by swapping the spins. Recalc Energy
            swap_spins(spin, orig_site, dest_site);
            const int E1 = calcE_for_two_connected_sites(lattice, spin,
                                                         orig_site, dest_site);

            const int dE = E1 - E0;

            // Monte-Carlo acceptance criteria
            double p = std::min(1.0, std::exp(-beta*dE));
            double r = rand_01(generator);
            bool accepted(false);
            // Test if the move was REJECTED. If so, re-swap spins to undo move
            //if( rand_01(generator) > p){
            if( r > p){
                std::cerr << "Rejected!" << std::endl;
                swap_spins(spin, orig_site, dest_site);
            }else{
                accepted = true;
                std::cerr << "dE = " << E1 << " - " << E0 << " ||| MC:"
                          << r << " < " << p << ". ";
                std::cerr << "Accepted!" << std::endl;
            }

            int M(0);
            for(int i = 0; i < lattice.nsites; ++i)
                M += spin[i];

            M_av += ((double)std::abs(M));
            ++n_av;

#if 1
            if(imove % 1 == 0 && accepted == true){
                int orig_ix = 1 + (orig_site/lattice.nz/lattice.ny)%lattice.nx;
                int orig_iy = 1 + (orig_site/lattice.nz)%lattice.ny;
                int dest_ix = 1 + (dest_site/lattice.nz/lattice.ny)%lattice.nx;
                int dest_iy = 1 + (dest_site/lattice.nz)%lattice.ny;
                std::stringstream ss;
	        ss << " <<<<<<< step: " << imove<< ", swapped [(" 
                   << orig_ix << ',' << orig_iy << ") <=> ("
                   << dest_ix << ',' << dest_iy << ")]"
                   << ". M_inst = "
                   << M/((double)noccupied)
                   << ". <M> = "
                   << ((double)M_av)/((double)n_av)/((double)noccupied)
                   << " >>>>>>>" << std::endl;

                print_cell(lattice, spin, ss.str());
            }
#endif
        }

        const int n_spins_to_flip(1*noccupied);
        int n_flipped_spins(0);
        do{ // Sample until we have flipped the required number of spins

            int initial_site(-1);
            do{ // Make sure that the site is occupied
                int tmp = rand_lattice_site(generator);
                if(spin[tmp] != 0)
                    initial_site = tmp;
            }while(initial_site < 0);

            // The "pocket" will contain the list of candidate sites whose spins
            // we want to attempt to fliip
            std::vector<int> pocket;
            pocket.push_back(initial_site);
            // The "cluster" is the MC cluster that is grown with acceptance
            // probability P to ensure detailed balance
            std::vector<int> cluster;
            cluster.push_back(initial_site);

            do{ // Grow cluster until it can no longer expand

                // Randomize the ordering of sites within the pocket
                std::random_shuffle(pocket.begin(), pocket.end());

                // Select the last element of the pocket and delete it
                int active_site = pocket.back();
                pocket.pop_back();

                // Now loop over the neighbors of the active_site
                for(size_t ibr = 0; ibr < lattice.neighbors[active_site].size();
                    ++ibr){
                    int nbr = lattice.neighbors[active_site][ibr];

                    // If the active site and neighbor have the same spin,
                    // && current neighbor not already in cluster
                    // && detailed balance is satisified ( rand < pCluster ) 
                    if(( spin[active_site] == spin[nbr] )
                       && ( std::find(cluster.begin(), cluster.end(), nbr)
                                                               == cluster.end())
                       && ( rand_01(generator) < pCluster))
                    {// Then add the neighbor to the cluster and pocket
                        cluster.push_back(nbr);
                        pocket.push_back(nbr);
                    }

                }

            }while(pocket.size() != 0);

            // Now, flip the spins of the entire cluster
            for(size_t i = 0; i < cluster.size(); ++i)
                spin[cluster[i]] *= -1;

            n_flipped_spins += cluster.size();
        }while(n_flipped_spins < n_spins_to_flip);

    }

    std::cout << std::setw(12) << T
	      << std::setw(12) << ((double)M_av)/((double)n_av)
						/((double)noccupied)
	      << std::setw(12) << n_av << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: localMC_ising_vacancies nx ny nz > a.dat" << std::endl;
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

    std::default_random_engine generator(19103);
    std::uniform_int_distribution<int> distribution(0,1);
#if 0
    for(int i = 0; i < lattice.nsites; ++i){
        if(i < lattice.nsites/2){
            initial_spin[i] = -1;
        }else{
            initial_spin[i] = 0;
        }
    }
#else
    for(int i = 0; i < lattice.nsites; ++i){
	if(distribution(generator) == 0)
	    initial_spin[i] = -1;
	else
	    initial_spin[i] = 1;
    }
#endif
    int mi = 0;
    for(int i = 0; i < lattice.nsites; ++i)
	mi += initial_spin[i];

    // Define the initial temperature and increments
    double T0 = 1.0e-16;
    double dT = 0.2;
    const size_t nT = 1;
    //const size_t nT = 16;

    for(size_t i = 0; i < nT; ++i){
	double T = T0 + dT*i;

        int spin[lattice.nsites];
        std::copy(initial_spin, initial_spin + lattice.nsites, spin);

	do_ising(lattice, spin, T, generator);
    }
}
