#ifdef ENABLE_MPI
#include <mpi.h>
#endif

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
#define FIXED_NUMBER yes

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

size_t nmove(10);
//size_t nmove(10000);
int ndesiredOccupied;

//----------------------------------------------------------------------------//

int read_command_line_arg(char* argv)
{
    int val;
    std::istringstream iss(argv);
    iss >> val; 
    if (!iss || !iss.eof()) {
        std::cerr << "could not convert '" << argv
            << "' to int" << std::endl;
        return EXIT_FAILURE;
    }
    return val;
}

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

int do_ising(ising::nodes& lattice, int* spin, const double T,
             std::default_random_engine& generator,
             double& M_av, double& numNeighbor_av, double& numVertNeighbor_av)
{
    int n_av(0);

    // Cluster expansion probability for a 2D ising model with all sites filled
    double beta = 1.0/T;
    const double pCluster = 1.0 - std::exp(-2.0 * beta);

    // Set up random numbers
    std::uniform_int_distribution<int> rand_lattice_site(0,lattice.nsites - 1);
    std::uniform_real_distribution<double> rand_01(0.0, 1.0);

    // Randomly zero out the spins for half the sites
    int noccupied(0);
    for(int i = 0; i < lattice.nsites; ++i)
        if(spin[i] != 0)
            ++noccupied;
    if(noccupied > (ndesiredOccupied)){
        do{

            int delete_this_site = rand_lattice_site(generator);
            spin[delete_this_site] = 0;

            noccupied = 0;
            for(int i = 0; i < lattice.nsites; ++i)
                if(spin[i] != 0)
                    ++noccupied;

        }while(noccupied > (ndesiredOccupied));
    }

    if(noccupied != ndesiredOccupied){
        std::cerr << "Function to create vacancies not working as expected.\n"
                  << " expected " << ndesiredOccupied << " occupied sites but"
                  << " found " << noccupied<< std::endl;
        exit(1);
    }

#if 0
    {
        std::string comment("Initial configuration");
        print_cell(lattice, spin, comment);
    }
#endif

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

        bool rejected(true);
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
            // Test if the move was REJECTED. If so, re-swap spins to undo move
            if( r > p){
                swap_spins(spin, orig_site, dest_site);
                rejected = false;
            }
        }

        // Calculate magnetization
        int M(0);
        for(int i = 0; i < lattice.nsites; ++i)
            M += spin[i];

        M_av += ((double)std::abs(M));
        ++n_av;

        // Calculate number density
        int numNeighbor(0);
        int numVertNeighbor(0);
        for(int i = 0; i < lattice.nsites; ++i){
            if(spin[i] != 0){
                for(size_t ibr = 0; ibr < lattice.neighbors[i].size(); ++ibr){
                    int nbr = lattice.neighbors[i][ibr];
                    if(spin[nbr] != 0){
                        ++numNeighbor;
                        if(lattice.neighborVertical[i][ibr])
                            ++numVertNeighbor;
                    }
                }
            }
        }

        numNeighbor_av += ((double)numNeighbor);
        numVertNeighbor_av += ((double)numVertNeighbor);
        // Incremented n_av above in calculation of magnetization;

#if 0
        if(imove % 1 == 0 && rejected == false){
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
               << ((double)M_av)/((double)n_av)/((double)noccupied) << ' '
               << ((double)numNeighbor)/((double)noccupied) << ' '
               << ((double)numNeighbor - numVertNeighbor)/((double)noccupied) << ' '
               << ((double)numVertNeighbor)/((double)noccupied) << ' '
               << " >>>>>>>" << std::endl;

            print_cell(lattice, spin, ss.str());
        }
#endif

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

    return n_av;

//    std::cout << std::setw(12) << T
//	      << std::setw(12) << ((double)M_av)/((double)n_av)
//						/((double)noccupied)
//	      << std::setw(12) << ((double)numNeighbor_av)/((double)n_av)
//						         /((double)noccupied)
//	      << std::setw(12) << ((double)numNeighbor_av - numVertNeighbor_av)
//                                                       /((double)numNeighbor_av)
//	      << std::setw(12) << ((double)numVertNeighbor_av)
//                                                       /((double)numNeighbor_av)
//	      << std::setw(12) << n_av << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr << "usage: clusterMC_ising_vacancies "
#ifdef FIXED_NUMBER
                  << "nx ny nz #_occupied > a.dat"
#else
                  << "nx ny nz \%_occupied > a.dat"
#endif
                  << std::endl;
        return EXIT_FAILURE;
    }

    // Define lattice size and initialize it
    int nx = read_command_line_arg(argv[1]);
    int ny = read_command_line_arg(argv[2]);
    int nz = read_command_line_arg(argv[3]);

    int my_rank(0), my_size(1);
#ifdef ENABLE_MPI    
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &my_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //num_angular_samples /= my_size;
#endif

#ifdef ENABLE_MPI
    if(my_rank == 0)
#endif
    // Now have all information needed to set up the lattice sites
    std::cout << "# clusterMC_ising_vacancies.cpp\n"
              << "# Initializing a " << nx << " x " << ny << " x "
                      << nz << " grid." << std::endl;

    ising::nodes lattice;
    lattice.init(nx, ny, nz);

#ifdef ENABLE_MPI
    if(my_rank == 0)
#endif
    lattice.report();

    // Determine how many sites are to be occupied
    {
        int tmp = read_command_line_arg(argv[4]);
        // If passed a negative number, make a fully occupied lattice
        if(tmp < 0){
            ndesiredOccupied = lattice.nsites;
        }else{
#ifdef FIXED_NUMBER
            ndesiredOccupied = tmp;
#else
            ndesiredOccupied = (int)(lattice.nsites/tmp);
#endif
        }
        if(ndesiredOccupied > lattice.nsites){
            std::cerr << "DYING!!! Asked for too many sites for this lattice\n"
                      << "# Requested " << ndesiredOccupied << " sites to be"
                      << " occupied out of " << lattice.nsites << " total sites"
                      << std::endl;
            return EXIT_FAILURE;
        }
#ifdef ENABLE_MPI
    if(my_rank == 0)
#endif
        std::cout << "# Requesting " << ndesiredOccupied << " sites to be"
                  << " occupied out of " << lattice.nsites << " total sites"
                  << std::endl;
    }

    // Set up the initial spins
    int initial_spin[lattice.nsites];

#ifdef ENABLE_MPI
    std::default_random_engine generator(19103 + 11*my_rank);
#else
    std::default_random_engine generator(19103);
#endif

    std::uniform_int_distribution<int> distribution(0,1);
    for(int i = 0; i < lattice.nsites; ++i)
	if(distribution(generator) == 0)
	    initial_spin[i] = -1;
	else
	    initial_spin[i] = 1;

    // Define the initial temperature and increments
    double T0 = 1.0e-16;
    double dT = 0.2;
    const size_t nT = 1;
    //const size_t nT = 32;

    for(size_t i = 0; i < nT; ++i){
	double T = T0 + dT*i;

        int spin[lattice.nsites];
        std::copy(initial_spin, initial_spin + lattice.nsites, spin);

        double my_M_av(0);
        double my_numNeighbor_av(0);
        double my_numVertNeighbor_av(0);

	int my_n_av = do_ising(lattice, spin, T, generator, my_M_av,
                               my_numNeighbor_av, my_numVertNeighbor_av);

        // This is checked in do_ising()
        const int noccupied = ndesiredOccupied;

#ifdef ENABLE_MPI
        double M_av;
        double numNeighbor_av;
        double numVertNeighbor_av;
        int n_av;

        std::cout << "On rank " << my_rank << ": " << my_M_av << " " << my_n_av
                  << std::endl;

        MPI_Reduce(&my_M_av, &M_av, 1,
                   MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&my_numNeighbor_av, &numNeighbor_av, 1,
                   MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&my_numVertNeighbor_av, &numVertNeighbor_av, 1,
                   MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&my_n_av, &n_av, 1,
                   MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(my_rank == 0)
            std::cout << "Total: " << M_av << " " << n_av
                << std::endl;


        // Not all processes are expected to have made the same number of
        // observations. Therefore, can't divide by my_size. Instead, divide
        // total quantity by total number of observations when printing

        if(my_rank == 0)
#else
        double M_av = my_M_av;
        double numNeighbor_av = my_numNeighbor_av;
        double numVertNeighbor_av = my_numVertNeighbor_av;
        int n_av = my_n_av;
#endif
        std::cout << std::setw(12) << T
    	          << std::setw(12) << ((double)M_av)/((double)n_av)
    	            				/((double)noccupied)
    	          << std::setw(12) << ((double)numNeighbor_av)/((double)n_av)
    	                                                    /((double)noccupied)
    	          << std::setw(12) <<((double)numNeighbor_av-numVertNeighbor_av)
                                                       /((double)numNeighbor_av)
    	          << std::setw(12) << ((double)numVertNeighbor_av)
                                                       /((double)numNeighbor_av)
    	          << std::setw(12) << n_av << std::endl;
    }
#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
}
