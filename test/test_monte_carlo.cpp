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
