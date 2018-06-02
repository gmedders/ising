#ifndef ISING_H
#define ISING_H

#include <cstdlib>
#include <iostream>
#include <vector>

namespace ising {

struct nodes {
  inline nodes();
  ~nodes();

  void init(int, int, int, double);
  void init(int, int, int);
  void report();
  void determine_connectivity();
  int find_site_index(int, int, int);
  void attempt_to_add(std::vector<int> &, std::vector<bool> &, int, int, int,
                      bool);

  double calcE_for_one_site(ising::nodes &, int &);
  double calcE_for_two_sites(ising::nodes &, int &, int &);

  int nx;
  int ny;
  int nz;
  double kInteraction;
  int nsites;

  bool *frozen;
  int *spin;

  std::vector<std::vector<int>> neighbors;
  std::vector<std::vector<bool>> neighborVertical;
  //    std::vector< std::vector<int> > diagNeighbors;
};

inline nodes::nodes()
    : nx(0), ny(0), nz(0), kInteraction(0), nsites(0), frozen(0), spin(0),
      neighbors(0), neighborVertical(0) {}

} // namespace ising

#endif // ISING_H
