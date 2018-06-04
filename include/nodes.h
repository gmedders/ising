#ifndef NODES_H
#define NODES_H

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "imaging.h"

namespace ising {

struct nodes {
  inline nodes();
  ~nodes();

  void init(int, int, int, double, std::unique_ptr<bounds>);
  void init(int, int, int, double);
  void init(int, int, int);
  void report();
  int calculate_noccupied();
  void determine_connectivity();
  int find_site_index(int, int, int);
  void attempt_to_add(std::vector<int> &, std::vector<bool> &, int, int, int,
                      bool);

  void generate_random_spins(std::default_random_engine &generator);

  double calcE_for_one_site(int &);
  double calcE_for_two_sites(int &, int &);

  int generate_desired_occupancy(int &);

  int nx;
  int ny;
  int nz;
  double kInteraction;
  int nsites;

  bool *frozen;
  int *spin;

  std::unique_ptr<bounds> Handler;

  std::vector<std::vector<int>> neighbors;
  std::vector<std::vector<bool>> neighborVertical;
};

inline nodes::nodes()
    : nx(0), ny(0), nz(0), kInteraction(0), nsites(0), frozen(0), spin(0),
      neighbors(0), neighborVertical(0) {}

} // namespace ising

#endif // NODES_H
