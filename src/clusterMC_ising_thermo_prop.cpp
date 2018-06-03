#include <cassert>
#include <cmath>
#include <cstdlib>

#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <algorithm>
#include <random>

#include "helpers.h"
#include "nodes.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

size_t nsample(10000);

//----------------------------------------------------------------------------//

void do_ising(ising::nodes &lattice, const double T,
              std::default_random_engine &generator) {
  int nsteps(0);
  int M_av(0);
  double E_av(0.0);
  double E2_av(0.0);
  double beta = 1.0 / T;

  // Cluster expansion probability for a 2D ising model with all sites filled
  const double p = 1.0 - std::exp(-2.0 * beta);

  // Set up random numbers
  std::uniform_int_distribution<int> rand_lattice_site(0, lattice.nsites - 1);
  std::uniform_real_distribution<double> rand_01(0.0, 1.0);

  // The sampling is determined by the total number of spins that have been
  // flipped. Required number of spins grows with number of lattice sites
  int n_spins_to_flip(lattice.nsites * nsample);
  int n_flipped_spins(0);

  do { // Sample until we have flipped the required number of spins

    int initial_site = rand_lattice_site(generator);

    // The "pocket" will contain the list of candidate sites whose spins
    // we want to attempt to fliip
    std::vector<int> pocket;
    pocket.push_back(initial_site);
    // The "cluster" is the MC cluster that is grown with acceptance
    // probability P to ensure detailed balance
    std::vector<int> cluster;
    cluster.push_back(initial_site);

    do { // Grow cluster until it can no longer expand

      // Randomize the ordering of sites within the pocket
      std::random_shuffle(pocket.begin(), pocket.end());

      // Select the last element of the pocket and delete it
      int active_site = pocket.back();
      pocket.pop_back();

      // Now loop over the neighbors of the active_site
      for (size_t ibr = 0; ibr < lattice.neighbors[active_site].size(); ++ibr) {
        int nbr = lattice.neighbors[active_site][ibr];

        // If the active site and neighbor have the same spin,
        // && current neighbor has not already been placed in cluster
        // && detailed balance is satisified ( rand < p )
        if ((lattice.spin[active_site] == lattice.spin[nbr]) &&
            (std::find(cluster.begin(), cluster.end(), nbr) == cluster.end()) &&
            (rand_01(generator) <
             p)) { // Then add the neighbor to the cluster and pocket
          cluster.push_back(nbr);
          pocket.push_back(nbr);
        }
      }

    } while (pocket.size() != 0);

    // Now, flip the spins of the entire cluster
    for (size_t i = 0; i < cluster.size(); ++i)
      lattice.spin[cluster[i]] *= -1;

    n_flipped_spins += cluster.size();

    // Initial magnetization
    int M(0);
    for (int i = 0; i < lattice.nsites; ++i)
      M += lattice.spin[i];

    M_av += std::abs(M);
    ++nsteps;

  } while (n_flipped_spins < n_spins_to_flip);

  double Cv = beta * beta / lattice.nsites * (E2_av - E_av * E_av);
  std::cout << std::setw(12) << T << std::setw(12)
            << ((double)M_av) / ((double)nsteps) / ((double)lattice.nsites)
            << std::setw(12) << n_flipped_spins / nsteps << std::setw(12)
            << nsteps << std::setw(12) << Cv << std::setw(12) << nsample
            << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "usage: clusterMC_ising_thermo_prop nx ny nz > a.dat"
              << std::endl;
    return EXIT_FAILURE;
  }

  // Define lattice size and initialize it
  int nx(0), ny(0), nz(0);
  try {
    nx = ising::read_command_line_int(argv[1]);
    ny = ising::read_command_line_int(argv[2]);
    nz = ising::read_command_line_int(argv[3]);
  } catch (const std::exception& e){
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "# clusterMC_ising_thermo_prop.cpp\n"
            << "# Initializing a " << nx << " x " << ny << " x " << nz
            << " grid." << std::endl;

  ising::nodes lattice;
  lattice.init(nx, ny, nz);

  std::default_random_engine generator;
  lattice.generate_random_spins(generator);

  // Save the initial spins
  int initial_spin[lattice.nsites];
  std::copy(lattice.spin, lattice.spin + lattice.nsites, initial_spin);

  // Define the initial temperature and increments
  double T0 = 1.00;
  double dT = 0.2;
  const size_t nT = 32;

  for (size_t i = 0; i < nT; ++i) {
    double T = T0 + dT * i;

    std::copy(initial_spin, initial_spin + lattice.nsites, lattice.spin);
    do_ising(lattice, T, generator);
  }
}
