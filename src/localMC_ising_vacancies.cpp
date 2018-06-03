#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <algorithm>
#include <random>

#include "helpers.h"
#include "nodes.h"

#define FIXED_NUMBER yes
#define CONDENSED yes
#define FROZEN yes

#define VERBOSE yes
////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

// size_t nmove(10);
size_t nmove(40000);
int ndesiredOccupied;
int my_rank(0), my_size(1);

//----------------------------------------------------------------------------//

int do_ising(ising::nodes &lattice, const double T,
             std::default_random_engine &generator, double &M_av,
             double &numNeighbor_av, double &numVertNeighbor_av) {
  int n_av(0);

  double beta = 1.0 / T;

  // Set up random numbers
  std::uniform_int_distribution<int> rand_lattice_site(0, lattice.nsites - 1);
  std::uniform_real_distribution<double> rand_01(0.0, 1.0);

  // Randomly zero out the spins for half the sites
  int noccupied = ising::generated_desired_occupancy(lattice, ndesiredOccupied);

#ifdef VERBOSE
#ifdef ENABLE_MPI
  if (my_rank == 0)
#endif
  {
    std::string comment("Initial configuration, T = ");
    comment += std::to_string(T);
    ising::print_cell(lattice, comment);
  }
#endif

  // Required number of attempted moves to grow with number of lattice sites
  for (size_t imove = 0; imove < nmove * lattice.nsites; ++imove) {

    bool do_move = true;
    // For move, select two random sites (orig, dest) of different spin
    int orig_site = rand_lattice_site(generator);
    int dest_site = rand_lattice_site(generator);
    if (orig_site == dest_site ||
        (lattice.spin[orig_site] == lattice.spin[dest_site])) {
      do_move = false;
    }

    // Do not allow moves between sites when one is frozen
    if (lattice.frozen[orig_site] || lattice.frozen[dest_site])
      do_move = false;

    // If move is allowed
    if (do_move) {
      // Calculate the energy before the move
      const double E0 = lattice.calcE_for_two_sites(orig_site, dest_site);

      // Create trial move by swapping the spins. Recalc Energy
      ising::swap_spins(lattice.spin, orig_site, dest_site);
      const double E1 = lattice.calcE_for_two_sites(orig_site, dest_site);

      const double dE = E1 - E0;

      // Monte-Carlo acceptance criteria
      double p = std::exp(-beta * dE);
      double r = rand_01(generator);
      // Test if the move was REJECTED. If so, re-swap spins to undo move
      if (r > p) {
        ising::swap_spins(lattice.spin, orig_site, dest_site);
      }
    }

    const int n_spins_to_flip(1 * noccupied);
    for (int n = 0; n < n_spins_to_flip; ++n) {

      int active_site = rand_lattice_site(generator);

      // Calculate energy before the spin flip
      double E0 = lattice.calcE_for_one_site(active_site);

      // Create trial move by swapping the spins. Recalc Energy
      lattice.spin[active_site] *= -1;
      double E1 = lattice.calcE_for_one_site(active_site);

      const double dE = E1 - E0;

      // Monte-Carlo acceptance criteria
      double p = std::exp(-beta * dE);
      double r = rand_01(generator);
      // Test if the move was REJECTED. If so, re-flip spin to undo move
      if (r > p) {
        lattice.spin[active_site] *= -1;
      }
    }

    ising::collect_stats(lattice, n_av, M_av, numNeighbor_av,
                         numVertNeighbor_av);
  }

#ifdef VERBOSE
#ifdef ENABLE_MPI
  if (my_rank == 0)
#endif
  {
    std::string comment("Final configuration, T = ");
    comment += std::to_string(T);
    ising::print_cell(lattice, comment);
  }
#endif

  return n_av;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cerr << "usage: localMC_ising_vacancies "
#ifdef FIXED_NUMBER
              << "nx ny nz kInteraction #_occupied > a.dat"
#else
              << "nx ny nz kInteraction \%_occupied > a.dat"
#endif
              << std::endl;
    return EXIT_FAILURE;
  }

  // Define lattice size and initialize it
  int nx = ising::read_command_line_int(argv[1]);
  int ny = ising::read_command_line_int(argv[2]);
  int nz = ising::read_command_line_int(argv[3]);

#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &my_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  nmove /= my_size;
#endif

  double kInteraction = ising::read_command_line_double(argv[4]);

#ifdef ENABLE_MPI
  if (my_rank == 0)
#endif
    // Now have all information needed to set up the lattice sites
    std::cout << "# localMC_ising_vacancies.cpp\n"
              << "# Initializing a " << nx << " x " << ny << " x " << nz
              << " grid." << std::endl;

  ising::nodes lattice;
  lattice.init(nx, ny, nz, kInteraction);

  // Determine how many sites are to be occupied
  {
    int tmp = ising::read_command_line_int(argv[5]);
    // If passed a negative number, make a fully occupied lattice
    if (tmp < 0) {
      ndesiredOccupied = lattice.nsites;
    } else {
#ifdef FIXED_NUMBER
      ndesiredOccupied = tmp;
#else
      ndesiredOccupied = (int)(lattice.nsites / tmp);
#endif
    }
    if (ndesiredOccupied > lattice.nsites) {
      std::cerr << "DYING!!! Asked for too many sites for this lattice\n"
                << "# Requested " << ndesiredOccupied << " sites to be"
                << " occupied out of " << lattice.nsites << " total sites"
                << std::endl;
      return EXIT_FAILURE;
    }
#ifdef ENABLE_MPI
    if (my_rank == 0)
#endif
      std::cout << "# Requesting " << ndesiredOccupied << " sites to be"
                << " occupied out of " << lattice.nsites << " total sites"
                << std::endl;
  }

#ifdef ENABLE_MPI
  if (my_rank == 0)
#endif
    lattice.report();

  // Set up the initial spins
  int initial_spin[lattice.nsites];
  std::fill(initial_spin, initial_spin + lattice.nsites, 0);

#ifdef ENABLE_MPI
  std::default_random_engine generator(19103 + 11 * my_rank);
#else
  std::default_random_engine generator(19103);
#endif

#ifdef FROZEN
#ifdef ENABLE_MPI
  if (my_rank == 0)
#endif
    std::cout << "# Freezing occupancy of bottom layer" << std::endl;
  // Freeze occupancy of bottom layer
  for (int iz = 0; iz < 1; ++iz) {
    for (int iy = 0; iy < lattice.ny; ++iy) {
      for (int ix = 0; ix < lattice.nx; ++ix) {
        lattice.frozen[lattice.find_site_index(ix, iy, iz)] = true;
      }
    }
  }
#endif

#ifdef CONDENSED
  int num_inited_sites(0);
  for (int iz = 0; iz < lattice.nz; ++iz) {
    for (int iy = 0; iy < lattice.ny; ++iy) {
      for (int ix = 0; ix < lattice.nx; ++ix) {

        if (num_inited_sites < ndesiredOccupied) {
          initial_spin[lattice.find_site_index(ix, iy, iz)] = 1;
          ++num_inited_sites;
        }
      }
    }
  }
#else
  std::uniform_int_distribution<int> distribution(0, 1);
  for (int i = 0; i < lattice.nsites; ++i)
    if (distribution(generator) == 0)
      initial_spin[i] = -1;
    else
      initial_spin[i] = 1;
#endif

  // Define the initial temperature and increments
  double T0 = 1.0e-16;
  double dT = 0.2;
  const size_t nT = 32;

  for (size_t i = 0; i < nT; ++i) {
    double T = T0 + dT * i;

    std::copy(initial_spin, initial_spin + lattice.nsites, lattice.spin);

    double my_M_av(0);
    double my_numNeighbor_av(0);
    double my_numVertNeighbor_av(0);

    int my_n_av = do_ising(lattice, T, generator, my_M_av, my_numNeighbor_av,
                           my_numVertNeighbor_av);

    // This is checked in do_ising()
    const int noccupied = ndesiredOccupied;

#ifdef ENABLE_MPI
    double M_av;
    double numNeighbor_av;
    double numVertNeighbor_av;
    int n_av;

    MPI_Reduce(&my_M_av, &M_av, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&my_numNe ighbor_av, &numNeighbor_av, 1, MPI_DOUBLE_PRECISION,
               MPI_SUM, 0, MPI_COMM_WOR
LD);
    MPI_Reduce(&my_numVertNeighbor_av, &numVertNeighbor_av, 1,
               MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_n_av, &n_av, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Not all processes are expected to have made the same number of
    // observations. Therefore, can't divide by my_size. Instead, divide
    // total quantity by total number of observations when printing

    if (my_rank == 0)
#else
    double M_av = my_M_av;
    double numNeighbor_av = my_numNeighbor_av;
    double numVertNeighbor_av = my_numVertNeighbor_av;
    int n_av = my_n_av;
#endif
      std::cout << std::setw(12) << T << std::setw(12)
                << ((double)M_av) / ((double)n_av) / ((double)noccupied)
                << std::setw(12)
                << ((double)numNeighbor_av) / ((double)n_av) /
                       ((double)noccupied)
                << std::setw(12)
                << ((double)numNeighbor_av - numVertNeighbor_av) /
                       ((double)numNeighbor_av)
                << std::setw(12)
                << ((double)numVertNeighbor_av) / ((double)numNeighbor_av)
                << std::setw(12) << n_av << std::endl;
  }
#ifdef ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
}
