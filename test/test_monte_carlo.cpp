#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "monte_carlo.h"
#include "gtest/gtest.h"

TEST(nodes, ParticleSwap) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  auto lattice = std::make_shared<ising::nodes>();
  int nx = 2, ny = 2, nz = 2;
  lattice->init(nx, ny, nz);
  // initial_spin
  // level 0:
  // -1  1
  // -1 -1
  // level 1:
  //  1   1
  //  1  -1
  std::vector<int> initial_spin = {-1, 1, -1, -1, 1, 1, 1, -1};
  std::copy(&initial_spin[0], &initial_spin[0] + lattice->nsites,
            lattice->spin);

  const double T = 10000.0;
  auto this_mc_move = std::make_unique<ising::particle_swap>(T);

  EXPECT_DOUBLE_EQ(this_mc_move->beta, 1.0 / T);

  std::default_random_engine generator(19104);
  std::uniform_int_distribution<int> rand_lattice_site(0, lattice->nsites - 1);
  this_mc_move->step(lattice, generator, rand_lattice_site);

  // Particles 1 and 7 should have swapped. Verify this
  EXPECT_EQ(initial_spin.at(1), lattice->spin[7]);
}

// Same test as ParticleSwap, but mediated through monte_carlo class
TEST(nodes, MonteCarloAddParticleSwapMove) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  auto lattice = std::make_shared<ising::nodes>();
  int nx = 2, ny = 2, nz = 2;
  lattice->init(nx, ny, nz);
  // initial_spin
  // level 0:
  // -1  1
  // -1 -1
  // level 1:
  //  1   1
  //  1  -1
  std::vector<int> initial_spin = {-1, 1, -1, -1, 1, 1, 1, -1};
  std::copy(&initial_spin[0], &initial_spin[0] + lattice->nsites,
            lattice->spin);

  int seed(19104);
  const double T = 10000.0;
  ising::monte_carlo mc_simulation(seed, lattice, T);

  mc_simulation.add_mc_move("particle_swap");

  mc_simulation.do_n_steps(1);

  // Particles 1 and 7 should have swapped. Verify this
  EXPECT_EQ(initial_spin.at(1), lattice->spin[7]);
}

TEST(nodes, LocalSpinFlip) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  auto lattice = std::make_shared<ising::nodes>();
  int nx = 2, ny = 2, nz = 2;
  lattice->init(nx, ny, nz);
  // initial_spin
  // level 0:
  // -1  1
  // -1 -1
  // level 1:
  //  1   1
  //  1  -1
  std::vector<int> initial_spin = {-1, 1, -1, -1, 1, 1, 1, -1};
  std::copy(&initial_spin[0], &initial_spin[0] + lattice->nsites,
            lattice->spin);

  const double T = 10000.0;
  auto this_mc_move = std::make_unique<ising::local_mc_spin_flip>(T);

  EXPECT_DOUBLE_EQ(this_mc_move->beta, 1.0 / T);

  std::default_random_engine generator(19104);
  std::uniform_int_distribution<int> rand_lattice_site(0, lattice->nsites - 1);
  this_mc_move->step(lattice, generator, rand_lattice_site);

  // Particle 1 spin should have flipped. Verify this
  EXPECT_EQ(initial_spin.at(1), -1 * lattice->spin[1]);
}

// Same test as ClusterSpinFlip, but mediated through monte_carlo class
TEST(nodes, MonteCarloAddLocalSpinFlipMove) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  auto lattice = std::make_shared<ising::nodes>();
  int nx = 2, ny = 2, nz = 2;
  lattice->init(nx, ny, nz);
  // initial_spin
  // level 0:
  // -1  1
  // -1 -1
  // level 1:
  //  1   1
  //  1  -1
  std::vector<int> initial_spin = {-1, 1, -1, -1, 1, 1, 1, -1};
  std::copy(&initial_spin[0], &initial_spin[0] + lattice->nsites,
            lattice->spin);

  int seed(19104);
  const double T = 10000.0;
  ising::monte_carlo mc_simulation(seed, lattice, T);

  mc_simulation.add_mc_move("local_mc_spin_flip");

  mc_simulation.do_n_steps(1);

  // Particle 1 spin should have flipped. Verify this
  EXPECT_EQ(initial_spin.at(1), -1 * lattice->spin[1]);
}

TEST(nodes, ClusterSpinFlip) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  auto lattice = std::make_shared<ising::nodes>();
  int nx = 2, ny = 2, nz = 2;
  lattice->init(nx, ny, nz);
  // initial_spin
  // level 0:
  // -1  1
  // -1 -1
  // level 1:
  //  1   1
  //  1  -1
  std::vector<int> initial_spin = {-1, 1, -1, -1, 1, 1, 1, -1};
  std::copy(&initial_spin[0], &initial_spin[0] + lattice->nsites,
            lattice->spin);

  const double T = 10000.0;
  auto this_mc_move = std::make_unique<ising::cluster_mc_spin_flip>(T);

  EXPECT_DOUBLE_EQ(this_mc_move->pCluster, 1.0 - std::exp(-2.0 / T));

  std::default_random_engine generator(19104);
  std::uniform_int_distribution<int> rand_lattice_site(0, lattice->nsites - 1);
  this_mc_move->step(lattice, generator, rand_lattice_site);

  // I've stored the state of the spins after a single pass through
  // cluster_mc_spin_flip.
  std::vector<int> final_spin = {-1, -1, -1, 1, -1, -1, -1, 1};

  // Assert that the lattice reflects these final_spins
  for (int i = 0; i < lattice->nsites; ++i) {
    EXPECT_EQ(final_spin.at(i), lattice->spin[i]);
  }
}

// Same test as ClusterSpinFlip, but mediated through monte_carlo class
TEST(nodes, MonteCarloAddClusterSpinFlipMove) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  auto lattice = std::make_shared<ising::nodes>();
  int nx = 2, ny = 2, nz = 2;
  lattice->init(nx, ny, nz);
  // initial_spin
  // level 0:
  // -1  1
  // -1 -1
  // level 1:
  //  1   1
  //  1  -1
  std::vector<int> initial_spin = {-1, 1, -1, -1, 1, 1, 1, -1};
  std::copy(&initial_spin[0], &initial_spin[0] + lattice->nsites,
            lattice->spin);

  int seed(19104);
  const double T = 10000.0;
  ising::monte_carlo mc_simulation(seed, lattice, T);

  mc_simulation.add_mc_move("cluster_mc_spin_flip");

  mc_simulation.do_n_steps(1);

  // I've stored the state of the spins after a single pass through
  // cluster_mc_spin_flip.
  std::vector<int> final_spin = {-1, -1, -1, 1, -1, -1, -1, 1};

  // Assert that the lattice reflects these final_spins
  for (int i = 0; i < lattice->nsites; ++i) {
    EXPECT_EQ(final_spin.at(i), lattice->spin[i]);
  }
}

TEST(nodes, MonteCarloAddTwoMoves) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  auto lattice = std::make_shared<ising::nodes>();
  int nx = 2, ny = 2, nz = 2;
  lattice->init(nx, ny, nz);
  // initial_spin
  // level 0:
  // -1  1
  // -1 -1
  // level 1:
  //  1   1
  //  1  -1
  std::vector<int> initial_spin = {-1, 1, -1, -1, 1, 1, 1, -1};
  std::copy(&initial_spin[0], &initial_spin[0] + lattice->nsites,
            lattice->spin);

  int seed(19104);
  const double T = 10000.0;
  ising::monte_carlo mc_simulation(seed, lattice, T);

  mc_simulation.add_mc_move("cluster_mc_spin_flip");
  mc_simulation.add_mc_move("particle_swap");

  int n_avg = mc_simulation.do_n_steps(1);

  // I've stored the state of the spins after a single pass through
  // cluster_mc_spin_flip.
  std::vector<int> final_spin = {-1, -1, -1, -1, -1, -1, 1, 1};

  // Assert that the lattice reflects these final_spins
  for (int i = 0; i < lattice->nsites; ++i) {
    EXPECT_EQ(final_spin.at(i), lattice->spin[i]);
  }
  EXPECT_EQ(n_avg, 1);
  EXPECT_EQ(mc_simulation.M_av, 4);
  EXPECT_EQ(mc_simulation.numNeighbor_av, 48);
  EXPECT_EQ(mc_simulation.numVertNeighbor_av, 16);

  n_avg = mc_simulation.do_n_steps(1);

  EXPECT_EQ(n_avg, 1);
  EXPECT_EQ(mc_simulation.M_av, 0);
  EXPECT_EQ(mc_simulation.numNeighbor_av, 48);
  EXPECT_EQ(mc_simulation.numVertNeighbor_av, 16);
}
