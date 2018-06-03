#include "nodes.h"

////////////////////////////////////////////////////////////////////////////////

namespace ising {

//----------------------------------------------------------------------------//

nodes::~nodes() {
  if (frozen != 0)
    delete[] frozen;
  if (spin != 0)
    delete[] spin;
}

//----------------------------------------------------------------------------//

void nodes::init(int my_nx, int my_ny, int my_nz) {
  nodes::init(my_nx, my_ny, my_nz, 0.0);
}

//----------------------------------------------------------------------------//

void nodes::init(int my_nx, int my_ny, int my_nz, double my_kInteraction) {
  std::unique_ptr<bounds> imaging = std::make_unique<pbc>(my_nx, my_ny, my_nz);
  init(my_nx, my_ny, my_nz, my_kInteraction, std::move(imaging));
}

//----------------------------------------------------------------------------//

void nodes::init(int my_nx, int my_ny, int my_nz, double my_kInteraction,
                 std::unique_ptr<bounds> imaging) {
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

  Handler = std::move(imaging);

  // Determine neighbor connectivity
  determine_connectivity();
}

//----------------------------------------------------------------------------//

void nodes::generate_random_spins(std::default_random_engine &generator) {
  std::uniform_int_distribution<int> distribution(0, 1);
  for (int i = 0; i < nsites; ++i) {
    if (distribution(generator) == 0)
      spin[i] = -1;
    else
      spin[i] = 1;
  }
}


int nodes::calculate_noccupied() {
  int noccupied = 0;
  for (int i = 0; i < nsites; ++i)
    if (spin[i] != 0)
      ++noccupied;
  return noccupied;
}

//----------------------------------------------------------------------------//

void nodes::report() {
  Handler->report();
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
  return Handler->find_site_index(ix, iy, iz);
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
