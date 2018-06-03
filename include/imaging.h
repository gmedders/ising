#ifndef IMAGING_H
#define IMAGING_H

#include <cstdlib>
#include <iostream>
#include <memory>

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

} // namespace ising

#endif // IMAGING_H
