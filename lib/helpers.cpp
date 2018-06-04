#include "helpers.h"

#include <sstream>
#include <vector>

namespace ising {

//----------------------------------------------------------------------------//

int read_command_line_int(const char *argv) {
  int val;
  std::istringstream iss(argv);
  iss >> val;
  if (!iss || !iss.eof()) {
    std::ostringstream os;
    os << "could not convert '" << argv << "' to int";
    throw std::invalid_argument(os.str());
  }
  return val;
}

//----------------------------------------------------------------------------//

double read_command_line_double(const char *argv) {
  double val;
  std::istringstream iss(argv);
  iss >> val;
  if (!iss || !iss.eof()) {
    std::ostringstream os;
    os << "could not convert '" << argv << "' to double";
    throw std::invalid_argument(os.str());
  }
  return val;
}

//----------------------------------------------------------------------------//

void swap_spins(int *spin, int &orig, int &dest) {
  int tmp = spin[dest];
  spin[dest] = spin[orig];
  spin[orig] = tmp;
}

//----------------------------------------------------------------------------//

void print_cell(ising::nodes &lattice, std::string comment) {
  std::cerr << comment << std::endl;
  int *spin = lattice.spin;
  bool *frozen = lattice.frozen;
  for (int iz = 0; iz < lattice.nz; ++iz) {
    std::cerr << " <<< Layer " << iz << " >>>" << std::endl;
    for (int ix = 0; ix < lattice.nx; ++ix) {
      for (int iy = 0; iy < lattice.ny; ++iy) {
        int site = lattice.find_site_index(ix, iy, iz);
        char front = ' ';
        char back = ' ';
        if (frozen[site]) {
          front = '[';
          back = ']';
        }
        if (spin[site] < 0)
          std::cerr << front << '-' << back << ' ';
        else if (spin[site] == 0)
          std::cerr << front << 'O' << back << ' ';
        else
          std::cerr << front << '+' << back << ' ';
      }
      std::cerr << std::endl;
    }
  }
}

//----------------------------------------------------------------------------//

void collect_stats(ising::nodes &lattice, int &n_av, double &M_av,
                   double &numNeighbor_av, double &numVertNeighbor_av) {
  // Calculate magnetization
  int M(0);
  for (int i = 0; i < lattice.nsites; ++i)
    M += lattice.spin[i];

  M_av += ((double)std::abs(M));
  ++n_av;

  // Calculate number density
  int numNeighbor(0);
  int numVertNeighbor(0);
  for (int i = 0; i < lattice.nsites; ++i) {
    if (lattice.spin[i] != 0) {
      for (size_t ibr = 0; ibr < lattice.neighbors[i].size(); ++ibr) {
        int nbr = lattice.neighbors[i][ibr];
        if (lattice.spin[nbr] != 0) {
          ++numNeighbor;
          if (lattice.neighborVertical[i][ibr])
            ++numVertNeighbor;
        }
      }
    }
  }

  numNeighbor_av += ((double)numNeighbor);
  numVertNeighbor_av += ((double)numVertNeighbor);
}

//----------------------------------------------------------------------------//

int generated_desired_occupancy(ising::nodes &lattice, int &ndesiredOccupied) {

  int noccupied = lattice.calculate_noccupied();
  if (noccupied <= (ndesiredOccupied))
    return noccupied;

  std::vector<std::vector<int>::size_type> indices;
  for (int i = 0; i < lattice.nsites; ++i)
    indices.push_back(i);
  std::random_shuffle(indices.begin(), indices.end());

  while (noccupied > ndesiredOccupied && !indices.empty()) {

    int delete_this_site = indices.back();
    indices.pop_back();
    if (not lattice.frozen[delete_this_site]) {
      lattice.spin[delete_this_site] = 0;
      noccupied = lattice.calculate_noccupied();
    }
  }

  return noccupied;
}

//----------------------------------------------------------------------------//

} // namespace ising
