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

#include "tclap/CmdLine.h"

#include "helpers.h"
#include "monte_carlo.h"
#include "nodes.h"

#define FROZEN yes
#define VERBOSE yes

////////////////////////////////////////////////////////////////////////////////

namespace {

int my_rank(0), my_size(1);

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // if (argc != 6) {
  //   std::cerr << "usage: clusterMC_ising_vacancies "
  //             << "nx ny nz k_interaction occupancy_ratio > a.dat" <<
  //             std::endl;
  //   return EXIT_FAILURE;
  // }

  // Define lattice size and initialize it
  int nx, ny, nz, nsteps;
  double k_interaction, occupancy;
  bool localMC;
  try {
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
    TCLAP::ValueArg<int> x_arg("x", "nx", "Number of sites in x-dimension",
                               true, 1, "ndim_x", cmd);
    TCLAP::ValueArg<int> y_arg("y", "ny", "Number of sites in y-dimension",
                               true, 1, "ndim_y", cmd);
    TCLAP::ValueArg<int> z_arg("z", "nz", "Number of sites in z-dimension",
                               false, 1, "ndim_z", cmd);
    TCLAP::ValueArg<int> n_arg("n", "nsteps", "Number of MC simulation steps",
                               true, 0, "nsteps", cmd);
    TCLAP::ValueArg<double> k_arg("k", "k_interaction",
                                  "Parameter that controls the surface "
                                  "tension. Only useful with vacancies.",
                                  false, 0.0, "Surface tension", cmd);
    TCLAP::ValueArg<double> o_arg("o", "occupancy",
                                  "Fractional occupancy of the lattice", false,
                                  -1, "occupancy", cmd);
    TCLAP::SwitchArg l_arg(
        "l", "localMC",
        "Do brute-force local Monte-Carlo for individual spins.", cmd, false);
    cmd.parse(argc, argv);
    nx = x_arg.getValue();
    ny = y_arg.getValue();
    nz = z_arg.getValue();
    nsteps = n_arg.getValue();
    k_interaction = k_arg.getValue();
    occupancy = o_arg.getValue();
    localMC = l_arg.getValue();

    if (occupancy > 1.0)
      throw std::invalid_argument("Occupancy ratio must be <= 1.0");

  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &my_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  nsteps /= my_size;
#endif

  auto lattice = std::make_shared<ising::nodes>();
  lattice->init(nx, ny, nz, k_interaction);

  // Determine how many sites are to be occupied
  // If passed a negative number, make a fully occupied lattice
  int ndesiredOccupied = lattice->nsites;
  if (occupancy >= 0)
    ndesiredOccupied = (int)(lattice->nsites * occupancy);

#ifdef FROZEN
  // Freeze occupancy of bottom layer
  for (int iz = 0; iz < 1; ++iz) {
    for (int iy = 0; iy < lattice->ny; ++iy) {
      for (int ix = 0; ix < lattice->nx; ++ix) {
        lattice->frozen[lattice->find_site_index(ix, iy, iz)] = true;
      }
    }
  }
#endif

  int seed = 19103 + 11 * my_rank;
  std::default_random_engine generator(seed);

  lattice->generate_random_spins(generator);
  const int noccupied = lattice->generate_desired_occupancy(ndesiredOccupied);

  // Save the initial spins
  int initial_spin[lattice->nsites];
  std::copy(lattice->spin, lattice->spin + lattice->nsites, initial_spin);

  // Define the initial temperature and increments
  double T0 = 1.0e-16;
  double dT = 0.2;
  const size_t nT = 32;

  for (size_t i = 0; i < nT; ++i) {
    double T = T0 + dT * i;

    std::copy(initial_spin, initial_spin + lattice->nsites, lattice->spin);

    seed += 11 * my_rank;
    ising::monte_carlo mc_simulation(seed, lattice, T);

    if (occupancy >= 0.0)
      mc_simulation.add_mc_move("particle_swap");

    if (localMC)
      mc_simulation.add_mc_move("local_mc_spin_flip");
    else
      mc_simulation.add_mc_move("cluster_mc_spin_flip");

    int my_n_av = mc_simulation.do_n_steps(nsteps);

    double my_M_av(mc_simulation.M_av);
    double my_numNeighbor_av(mc_simulation.numNeighbor_av);
    double my_numVertNeighbor_av(mc_simulation.numVertNeighbor_av);

#ifdef ENABLE_MPI
    double M_av;
    double numNeighbor_av;
    double numVertNeighbor_av;
    int n_av;

    MPI_Reduce(&my_M_av, &M_av, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&my_numNeighbor_av, &numNeighbor_av, 1, MPI_DOUBLE_PRECISION,
               MPI_SUM, 0, MPI_COMM_WORLD);
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
