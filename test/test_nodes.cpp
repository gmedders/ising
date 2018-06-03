#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

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

// TEST(helpers, collectStatsNoVacancy2D) {
//   // initial_spin:
//   //  1  1 -1
//   // -1 -1  1
//   //  1  1  1
//   std::vector<int> initial_spin = {1, 1, -1, -1, -1, 1, 1, 1, 1};
//
//   // nx = 3, ny = 3, nz = 1, kInteraction = 0
//   ising::nodes lattice;
//   lattice.init(3, 3, 1, 0.0);
//   ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));
//
//   std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);
//
//   int n_av(0);
//   double M_av(0.0);
//   double numNeighbor_av(0.0);
//   double numVertNeighbor_av(0.0);
//   collect_stats(lattice, n_av, M_av, numNeighbor_av, numVertNeighbor_av);
//
//   EXPECT_EQ(n_av, 1);
//   EXPECT_DOUBLE_EQ(numNeighbor_av, 36.0);
//   EXPECT_DOUBLE_EQ(M_av, 3.0);
//   EXPECT_DOUBLE_EQ(numVertNeighbor_av, 0.0);
// }
