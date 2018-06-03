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

  EXPECT_EQ(lattice.find_site_index(-1, 0, 0), 1);
  EXPECT_EQ(lattice.find_site_index(0, -1, 0), 2);
  EXPECT_EQ(lattice.find_site_index(0, 0, -1), 4);
  EXPECT_EQ(lattice.find_site_index(-1, -1, -1), 7);

  EXPECT_EQ(lattice.find_site_index(3, 0, 0), 1);
  EXPECT_EQ(lattice.find_site_index(0, 3, 0), 2);
  EXPECT_EQ(lattice.find_site_index(0, 0, 3), 4);
  EXPECT_EQ(lattice.find_site_index(3, 3, 3), 7);
}

TEST(nodes, FindSiteIndexNoPBC) {
  ising::nodes lattice;
  int nx = 2, ny = 2, nz = 2;
  auto imaging = std::make_unique<ising::no_pbc>(nx, ny, nz);
  lattice.init(nx, ny, nz, 0.0, std::move(imaging));
  EXPECT_EQ(lattice.find_site_index(0, 0, 0), 0);

  EXPECT_EQ(lattice.find_site_index(1, 0, 0), 1);
  EXPECT_EQ(lattice.find_site_index(0, 1, 0), 2);
  EXPECT_EQ(lattice.find_site_index(0, 0, 1), 4);
  EXPECT_EQ(lattice.find_site_index(1, 1, 1), 7);

  EXPECT_EQ(lattice.find_site_index(-1, 0, 0), -1);
  EXPECT_EQ(lattice.find_site_index(0, -1, 0), -1);
  EXPECT_EQ(lattice.find_site_index(0, 0, -1), -1);
  EXPECT_EQ(lattice.find_site_index(-1, -1, -1), -1);

  EXPECT_EQ(lattice.find_site_index(3, 0, 0), -1);
  EXPECT_EQ(lattice.find_site_index(0, 3, 0), -1);
  EXPECT_EQ(lattice.find_site_index(0, 0, 3), -1);
  EXPECT_EQ(lattice.find_site_index(3, 3, 3), -1);
}

TEST(nodes, FindSiteIndexSlab) {
  ising::nodes lattice;
  int nx = 2, ny = 2, nz = 2;
  auto imaging = std::make_unique<ising::slab>(nx, ny, nz);
  lattice.init(nx, ny, nz, 0.0, std::move(imaging));
  EXPECT_EQ(lattice.find_site_index(0, 0, 0), 0);

  EXPECT_EQ(lattice.find_site_index(1, 0, 0), 1);
  EXPECT_EQ(lattice.find_site_index(0, 1, 0), 2);
  EXPECT_EQ(lattice.find_site_index(0, 0, 1), 4);
  EXPECT_EQ(lattice.find_site_index(1, 1, 1), 7);

  EXPECT_EQ(lattice.find_site_index(-1, 0, 0), 1);
  EXPECT_EQ(lattice.find_site_index(0, -1, 0), 2);
  EXPECT_EQ(lattice.find_site_index(0, 0, -1), -1);
  EXPECT_EQ(lattice.find_site_index(-1, -1, -1), -1);

  EXPECT_EQ(lattice.find_site_index(3, 0, 0), 1);
  EXPECT_EQ(lattice.find_site_index(0, 3, 0), 2);
  EXPECT_EQ(lattice.find_site_index(0, 0, 3), -1);
  EXPECT_EQ(lattice.find_site_index(3, 3, 3), -1);
}
