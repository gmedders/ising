#ifndef NODES_H
#define NODES_H

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

namespace ising {

int pos_mod(int &a, int &b);

struct bounds {

  int nx;
  int ny;
  int nz;
  bounds(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz){};

  virtual int find_site_index(int, int, int) = 0;
  virtual void report() = 0;
};

struct pbc : public bounds {
  pbc(int nx, int ny, int nz) : bounds(nx, ny, nz){};
  int find_site_index(int ix, int iy, int iz) {
    return pos_mod(iz, nz) * ny * nx + pos_mod(iy, ny) * nx + pos_mod(ix, nx);
  }
  void report() {
    std::cerr << "# Full PBC" << std::endl;
    std::cout << "# Full PBC" << std::endl;
  }
};

struct no_pbc : public bounds {
  no_pbc(int nx, int ny, int nz) : bounds(nx, ny, nz){};
  int find_site_index(int ix, int iy, int iz) {
    if (ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz)
      return iz * ny * nx + iy * nx + ix;
    else
      return -1;
  }
  void report() {
    std::cerr << "# No PBC!!" << std::endl;
    std::cout << "# No PBC!!" << std::endl;
  }
};

struct slab : public bounds {
  slab(int nx, int ny, int nz) : bounds(nx, ny, nz){};
  int find_site_index(int ix, int iy, int iz) {
    if (iz >= 0 && iz < nz) // within the system (that is not replicated in z)
      return iz * ny * nx + pos_mod(iy, ny) * nx + pos_mod(ix, nx);
    else
      return -1;
  }
  void report() {
    std::cerr << "# Slab geometry: PBC only in XY" << std::endl;
    std::cout << "# Slab geometry: PBC only in XY" << std::endl;
  }
};

struct nodes {
  inline nodes();
  ~nodes();

  void init(int, int, int, double, std::unique_ptr<bounds>);
  void init(int, int, int, double);
  void init(int, int, int);
  void report();
  void determine_connectivity();
  int find_site_index(int, int, int);
  void attempt_to_add(std::vector<int> &, std::vector<bool> &, int, int, int,
                      bool);

  double calcE_for_one_site(int &);
  double calcE_for_two_sites(int &, int &);

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
