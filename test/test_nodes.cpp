#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "helpers.h"
#include "nodes.h"
#include "gtest/gtest.h"

TEST(nodes, AddConnectionOnlyIfSiteIsValid) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  int nx = 2, ny = 2, nz = 2;
  auto imaging = std::make_unique<ising::no_pbc>(nx, ny, nz);
  lattice.init(nx, ny, nz, 0.0, std::move(imaging));

  std::vector<int> connections;
  std::vector<bool> vertical;

  EXPECT_EQ(int(connections.size()), int(connections.size()));

  int prev_size = int(connections.size());
  lattice.attempt_to_add(connections, vertical, 0, 0, 0, false);

  EXPECT_EQ(int(connections.size()), int(connections.size()));
  EXPECT_EQ(int(connections.size()), prev_size + 1);

  // In no-PBC, the site (-1, -1, -1) is illegal, so there should be no increase
  // in the size of connections and vertical
  prev_size = int(connections.size());
  lattice.attempt_to_add(connections, vertical, -1, -1, -1, false);

  EXPECT_EQ(int(connections.size()), int(connections.size()));
  EXPECT_EQ(int(connections.size()), prev_size);
}

TEST(nodes, ConnectivityPBC) {
  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  int nx = 2, ny = 2, nz = 2;
  auto imaging = std::make_unique<ising::pbc>(nx, ny, nz);
  lattice.init(nx, ny, nz, 0.0, std::move(imaging));

  EXPECT_EQ(int(lattice.neighbors.size()), int(lattice.nsites));
  EXPECT_EQ(int(lattice.neighborVertical.size()), int(lattice.nsites));

  for (auto& neighbor: lattice.neighbors)
    EXPECT_EQ(int(neighbor.size()), 6);

  for (auto neighborVertical: lattice.neighborVertical){
    int n_vert(0);
    for (auto isVerticalNeighbor: neighborVertical) {
      if (isVerticalNeighbor)
        ++n_vert;
    }
    EXPECT_EQ(n_vert, 2);
  }

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
