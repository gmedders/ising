#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include "helpers.h"
#include "gtest/gtest.h"

TEST(helpers, ReadCommandLineIntValidInput) {
  EXPECT_EQ(ising::read_command_line_int("12"), 12);
}

TEST(helpers, ReadCommandLineIntInvalidInput) {
  EXPECT_THROW(ising::read_command_line_int("12.0"), std::invalid_argument);
  EXPECT_THROW(ising::read_command_line_int("abc"), std::invalid_argument);
  EXPECT_THROW(ising::read_command_line_int(""), std::invalid_argument);
}

TEST(helpers, ReadCommandLineDoubleValidInput) {
  EXPECT_EQ(ising::read_command_line_double("12"), 12.0);
  EXPECT_EQ(ising::read_command_line_double("12.0"), 12.0);
}

TEST(helpers, ReadCommandLineDoubleInvalidInput) {
  EXPECT_THROW(ising::read_command_line_double("abc"), std::invalid_argument);
  EXPECT_THROW(ising::read_command_line_double(""), std::invalid_argument);
}

TEST(helpers, SwapSpins) {
  int orig_site = 1;
  int dest_site = 2;
  std::vector<int> before = {-1, -1, 1, -1};
  std::vector<int> after = {-1, 1, -1, -1};
  ising::swap_spins(&before[0], orig_site, dest_site);
  EXPECT_EQ(before, after);
}

TEST(helpers, collectStatsNoVacancy2D) {
  // initial_spin:
  //  1  1 -1
  // -1 -1  1
  //  1  1  1
  std::vector<int> initial_spin = {1, 1, -1, -1, -1, 1, 1, 1, 1};

  // nx = 3, ny = 3, nz = 1, kInteraction = 0
  ising::nodes lattice;
  lattice.init(3, 3, 1, 0.0);
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));

  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  int n_av(0);
  double M_av(0.0);
  double numNeighbor_av(0.0);
  double numVertNeighbor_av(0.0);
  collect_stats(lattice, n_av, M_av, numNeighbor_av, numVertNeighbor_av);

  EXPECT_EQ(n_av, 1);
  EXPECT_DOUBLE_EQ(numNeighbor_av, 36.0);
  EXPECT_DOUBLE_EQ(M_av, 3.0);
  EXPECT_DOUBLE_EQ(numVertNeighbor_av, 0.0);
}

TEST(helpers, collectStatsTwiceNoVacancy2D) {
  // initial_spin:
  //  1  1 -1
  // -1 -1  1
  //  1  1  1
  std::vector<int> initial_spin = {1, 1, -1, -1, -1, 1, 1, 1, 1};

  // nx = 3, ny = 3, nz = 1, kInteraction = 0
  ising::nodes lattice;
  lattice.init(3, 3, 1, 0.0);
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));

  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  int n_av(0);
  double M_av(0.0);
  double numNeighbor_av(0.0);
  double numVertNeighbor_av(0.0);
  collect_stats(lattice, n_av, M_av, numNeighbor_av, numVertNeighbor_av);

  EXPECT_EQ(n_av, 1);
  EXPECT_DOUBLE_EQ(numNeighbor_av, 36.0);
  EXPECT_DOUBLE_EQ(M_av, 3.0);
  EXPECT_DOUBLE_EQ(numVertNeighbor_av, 0.0);

  // The stats are accumulated, so the values double when observed twice
  collect_stats(lattice, n_av, M_av, numNeighbor_av, numVertNeighbor_av);

  EXPECT_EQ(n_av, 2);
  EXPECT_DOUBLE_EQ(numNeighbor_av, 72.0);
  EXPECT_DOUBLE_EQ(M_av, 6.0);
  EXPECT_DOUBLE_EQ(numVertNeighbor_av, 0.0);
}

TEST(helpers, collectStatsNoVacancy3D) {
  // initial_spin
  // level 0:
  //  1  1
  // -1 -1
  // level 1:
  //  1  1
  //  1  1
  std::vector<int> initial_spin = {1, 1, -1, -1, 1, 1, 1, 1};

  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  lattice.init(2, 2, 2, 0.0);
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));

  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  int n_av(0);
  double M_av(0.0);
  double numNeighbor_av(0.0);
  double numVertNeighbor_av(0.0);
  collect_stats(lattice, n_av, M_av, numNeighbor_av, numVertNeighbor_av);

  EXPECT_EQ(n_av, 1);
  EXPECT_DOUBLE_EQ(numNeighbor_av, 48.0);
  EXPECT_DOUBLE_EQ(M_av, 4.0);
  EXPECT_DOUBLE_EQ(numVertNeighbor_av, 16.0);
}

TEST(helpers, collectStatsVacancy3D) {
  // initial_spin
  // level 0:
  //  1  0
  // -1 -1
  // level 1:
  //  1  0
  //  1  1
  std::vector<int> initial_spin = {1, 0, -1, -1, 1, 0, 1, 1};

  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  lattice.init(2, 2, 2, 0.0);
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));

  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  int n_av(0);
  double M_av(0.0);
  double numNeighbor_av(0.0);
  double numVertNeighbor_av(0.0);
  collect_stats(lattice, n_av, M_av, numNeighbor_av, numVertNeighbor_av);

  EXPECT_EQ(n_av, 1);
  EXPECT_DOUBLE_EQ(numNeighbor_av, 28.0);
  EXPECT_DOUBLE_EQ(M_av, 2.0);
  EXPECT_DOUBLE_EQ(numVertNeighbor_av, 12.0);
}

TEST(helpers, GenerateDesiredOccupancy) {
  // initial_spin
  // level 0:
  //  1  0
  // -1 -1
  // level 1:
  //  1  0
  //  1  1
  std::vector<int> initial_spin = {1, 0, -1, -1, 1, 0, 1, 1};

  // nx = 2, ny = 2, nz = 2, kInteraction = 0
  ising::nodes lattice;
  lattice.init(2, 2, 2, 0.0);
  ASSERT_EQ(int(initial_spin.size()), int(lattice.nsites));

  std::copy(&initial_spin[0], &initial_spin[0] + lattice.nsites, lattice.spin);

  // Set up random numbers
  std::default_random_engine generator(19103);
  std::uniform_int_distribution<int> rand_lattice_site(0, lattice.nsites - 1);

  int ndesiredOccupied(4);
  int noccupied = ising::generated_desired_occupancy(
      lattice, generator, rand_lattice_site, ndesiredOccupied);

  EXPECT_EQ(noccupied, ndesiredOccupied);
}
