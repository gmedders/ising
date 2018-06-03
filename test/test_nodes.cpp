#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "helpers.h"
#include "nodes.h"
#include "gtest/gtest.h"

TEST(nodes, PosMod) {
  int n_pos = 2;
  int pos = 0;
  EXPECT_EQ(ising::pos_mod(pos, n_pos), 0);

  pos = 1;
  EXPECT_EQ(ising::pos_mod(pos, n_pos), 1);

  pos = 2;
  EXPECT_EQ(ising::pos_mod(pos, n_pos), 0);

  pos = -1;
  EXPECT_EQ(ising::pos_mod(pos, n_pos), 1);
}

TEST(nodes, BoundsPBC) {
  ising::pbc imaging(2, 2, 2);
  EXPECT_EQ(imaging.find_site_index(0, 0, 0), 0);
  EXPECT_EQ(imaging.find_site_index(1, 0, 0), 1);
  EXPECT_EQ(imaging.find_site_index(0, 1, 0), 2);
  EXPECT_EQ(imaging.find_site_index(0, 0, 1), 4);
  EXPECT_EQ(imaging.find_site_index(1, 1, 1), 7);
}

TEST(nodes, FindSiteIndexFullPBC) {
  ising::nodes lattice;
  int nx = 2, ny = 2, nz = 2;
  auto imaging = std::make_unique<ising::pbc>(nx, ny, nz);
  lattice.init(nx, ny, nz, 0.0, std::move(imaging));
  EXPECT_EQ(lattice.find_site_index(0, 0, 0), 0);
  EXPECT_EQ(lattice.find_site_index(1, 0, 0), 1);
  EXPECT_EQ(lattice.find_site_index(0, 1, 0), 2);
  EXPECT_EQ(lattice.find_site_index(0, 0, 1), 4);
  EXPECT_EQ(lattice.find_site_index(1, 1, 1), 7);
}

TEST(nodes, EnergyForOneSite) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  lattice.init(2, 2, 2, 0.0);

  // initial_spin
  // level 0:
  //  1  0
  // -1 -1
  // level 1:
  //  1  0
  //  1  1
  std::vector<int> initial_spin = {1, 0, -1, -1, 1, 0, 1, 1};
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));
  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  int active_site = 2;
  double E0 = lattice.calcE_for_one_site(active_site);
  EXPECT_DOUBLE_EQ(E0, 2.0);
}

TEST(nodes, EnergyChangeOnSpinFlip) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  lattice.init(2, 2, 2, 0.0);

  // initial_spin
  // level 0:
  //  1  0
  // -1 -1
  // level 1:
  //  1  0
  //  1  1
  std::vector<int> initial_spin = {1, 0, -1, -1, 1, 0, 1, 1};
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));
  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  int active_site = 2;
  double E0 = lattice.calcE_for_one_site(active_site);
  EXPECT_DOUBLE_EQ(E0, 2.0);

  // Create trial move by swapping the spins. Recalc Energy
  lattice.spin[active_site] *= -1;
  double E1 = lattice.calcE_for_one_site(active_site);
  EXPECT_DOUBLE_EQ(E1, -2.0);
}

TEST(nodes, EnergyChangeOnMove) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  lattice.init(2, 2, 2, 0.0);

  // initial_spin
  // level 0:
  //  1  0
  // -1 -1
  // level 1:
  //  1  0
  //  1  1
  std::vector<int> initial_spin = {1, 0, -1, -1, 1, 0, 1, 1};
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));
  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  // Calculate the energy before the move
  int orig_site = 2;
  int dest_site = 4;
  double E0 = lattice.calcE_for_two_sites(orig_site, dest_site);
  EXPECT_DOUBLE_EQ(E0, -2.0);

  // Create trial move by swapping the spins. Recalc Energy
  ising::swap_spins(lattice.spin, orig_site, dest_site);
  double E1 = lattice.calcE_for_two_sites(orig_site, dest_site);
  EXPECT_DOUBLE_EQ(E1, 2.0);
}
