#include <cstdlib>

#include <cmath>
#include <random>

#include "helpers.h"
#include "nodes.h"

////////////////////////////////////////////////////////////////////////////////

namespace ising {

////////////////////////////////////////////////////////////////////////////////

struct mc_move {
  // std::uniform_real_distribution<double> rand_01(0.0, 1.0);
  mc_move() { rand_01 = std::uniform_real_distribution<double>(0.0, 1.0); };
  virtual void step(std::shared_ptr<ising::nodes>, std::default_random_engine &,
                    std::uniform_int_distribution<int> &) = 0;
  std::uniform_real_distribution<double> rand_01;
};

//----------------------------------------------------------------------------//

struct cluster_mc_spin_flip : public mc_move {

  double pCluster;
  cluster_mc_spin_flip(double T) : pCluster(1.0 - std::exp(-2.0 / T)){};

  void step(std::shared_ptr<ising::nodes>, std::default_random_engine &,
            std::uniform_int_distribution<int> &);
};

//----------------------------------------------------------------------------//

struct particle_swap : public mc_move {

  double beta;
  particle_swap(double T) : beta(1.0 / T){};

  void step(std::shared_ptr<ising::nodes>, std::default_random_engine &,
            std::uniform_int_distribution<int> &);
};

//----------------------------------------------------------------------------//

struct monte_carlo {

  double M_av;
  double numNeighbor_av;
  double numVertNeighbor_av;

  std::default_random_engine generator;
  std::shared_ptr<ising::nodes> lattice;

  std::vector<std::unique_ptr<mc_move>> mc_moves;

  void add_mc_move(std::string &, const double);
  int do_n_steps(int);
};

////////////////////////////////////////////////////////////////////////////////

} // namespace ising
