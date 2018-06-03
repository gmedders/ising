#include "nodes.h"

#define PBC true
//#define NO_PBC true

////////////////////////////////////////////////////////////////////////////////

namespace ising {

int pos_mod(int &a, int &b) { return ((a % b + b) % b); }

//----------------------------------------------------------------------------//

nodes::~nodes() {
  if (frozen != 0)
    delete[] frozen;
  if (spin != 0)
    delete[] spin;
}

//----------------------------------------------------------------------------//

void nodes::init(int my_nx, int my_ny, int my_nz, double my_kInteraction) {
  nx = my_nx;
  ny = my_ny;
  nz = my_nz;
  kInteraction = my_kInteraction;

  nsites = nx * ny * nz;

  if (nx < 1 || ny < 1 || nz < 1) {
    std::cerr << "The number of nodes in each dimension must be at least 1"
              << " you did: nx = " << nx << ", "
              << " ny = " << ny << ", "
              << " nz = " << nz << std::endl;
    exit(1);
  }

  // Allocate arrays
  frozen = new bool[nsites];
  std::fill(frozen, frozen + nsites, false); // Init to unfrozen occupancy

  spin = new int[nsites];
  std::fill(spin, spin + nsites, 0); // Initialize to unoccupied

  // Determine neighbor connectivity
  determine_connectivity();
}

//----------------------------------------------------------------------------//

void nodes::init(int my_nx, int my_ny, int my_nz) {
  nodes::init(my_nx, my_ny, my_nz, 0.0);
}

//----------------------------------------------------------------------------//

void nodes::report() {
#ifdef PBC
  std::cerr << "# Full PBC" << std::endl;
  std::cout << "# Full PBC" << std::endl;
#elif defined(NO_PBC)
  std::cerr << "# No PBC!!" << std::endl;
  std::cout << "# No PBC!!" << std::endl;
#else
  std::cerr << "# Slab geometry: PBC only in XY" << std::endl;
  std::cout << "# Slab geometry: PBC only in XY" << std::endl;
#endif
  if (kInteraction == 0) {
    std::cerr << "# Using standard Ising model Hamiltonian" << std::endl;
    std::cout << "# Using standard Ising model Hamiltonian" << std::endl;
  } else {
    std::cerr << "# Using Ising model with Surface Tension parameter k = "
              << kInteraction << std::endl;
    std::cout << "# Using Ising model with Surface Tension parameter k = "
              << kInteraction << std::endl;
  }
}

//----------------------------------------------------------------------------//

int nodes::find_site_index(int ix, int iy, int iz) {
#ifdef PBC
  return pos_mod(iz, nz) * ny * nx + pos_mod(iy, ny) * nx + pos_mod(ix, nx);
#elif defined(NO_PBC)
  // within the system (that is not replicated in x,y,z)
  if (ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz)
    return iz * ny * nx + iy * nx + ix;
  else
    return -1;
#else
  if (iz >= 0 && iz < nz) // within the system (that is not replicated in z)
    return iz * ny * nx + pos_mod(iy, ny) * nx + pos_mod(ix, nx);
  else
    return -1;
#endif
}

//----------------------------------------------------------------------------//

void nodes::attempt_to_add(std::vector<int> &vint, std::vector<bool> &vbool,
                           int i, int j, int k, bool vert) {
  int site = ising::nodes::find_site_index(i, j, k);
  if (site >= 0) {
    vint.push_back(site);
    vbool.push_back(vert);
  }
}

//----------------------------------------------------------------------------//

void nodes::determine_connectivity() {

  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {

        // First find the nearest neighbors
        {
          std::vector<int> connections;
          std::vector<bool> vertical;

          // +/- x
          attempt_to_add(connections, vertical, ix + 1, iy, iz, false);
          attempt_to_add(connections, vertical, ix - 1, iy, iz, false);

          // +/- y
          attempt_to_add(connections, vertical, ix, iy + 1, iz, false);
          attempt_to_add(connections, vertical, ix, iy - 1, iz, false);

          // +/- z
          // For 2D lattice, don't add +/- z
          if (nz > 1) {
            attempt_to_add(connections, vertical, ix, iy, iz + 1, true);
            attempt_to_add(connections, vertical, ix, iy, iz - 1, true);
          }

          neighbors.push_back(connections);
          neighborVertical.push_back(vertical);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------//

double nodes::calcE_for_one_site(int &site) {
  double E(0);
  for (size_t ibr = 0; ibr < neighbors[site].size(); ++ibr) {
    int nbr = neighbors[site][ibr];
    double ss = spin[site] * spin[nbr];
    E -= (ss + kInteraction * ss * ss);
  }
  return E;
}

//----------------------------------------------------------------------------//

double nodes::calcE_for_two_sites(int &site_a, int &site_b) {
  double E0_a = calcE_for_one_site(site_a);
  double E0_b = calcE_for_one_site(site_b);

  return E0_a + E0_b;
}

//----------------------------------------------------------------------------//

} // namespace ising
