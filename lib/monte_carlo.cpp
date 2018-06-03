#include <cstdlib>

#include <algorithm>
#include <cmath>
#include <random>

#include "helpers.h"
#include "nodes.h"

////////////////////////////////////////////////////////////////////////////////

namespace ising {

////////////////////////////////////////////////////////////////////////////////

class action {
  // std::shared_ptr<ising::nodes> lattice;
  // action(std::shared_ptr<ising::nodes> the_lattice) : lattice(the_lattice){};
  std::uniform_int_distribution<int> rand_lattice_site(0, lattice.nsites - 1);
  std::uniform_real_distribution<double> rand_01(0.0, 1.0);
};

//----------------------------------------------------------------------------//

class cluster_step : public action {

  double pCluster;
  cluster_step(double T) : pCluster(1.0 - std::exp(-2.0 / T)){};

  void step(ising::nodes &lattice, std::default_random_engine &generator) {
    const int n_spins_to_flip(lattice.nsites);
    int n_flipped_spins(0);
    // Sample until we have flipped the required number of spins
    do {

      // Draw random sites until an occupied one is found
      int initial_site(-1);
      do {
        int tmp = rand_lattice_site(generator);
        if (lattice.spin[tmp] != 0)
          initial_site = tmp;
      } while (initial_site < 0);

      // The "pocket" will contain the list of candidate sites whose spins
      // we want to attempt to fliip
      std::vector<int> pocket;
      pocket.push_back(initial_site);

      // The "cluster" is the MC cluster that is grown with acceptance
      // probability P to ensure detailed balance
      std::vector<int> cluster;
      cluster.push_back(initial_site);

      // Grow cluster until it can no longer expand
      do {

        // Randomize the ordering of sites within the pocket
        std::random_shuffle(pocket.begin(), pocket.end());

        // Select the last element of the pocket and delete it
        int active_site = pocket.back();
        pocket.pop_back();

        // Now loop over the neighbors of the active_site
        for (size_t ibr = 0; ibr < lattice.neighbors[active_site].size();
             ++ibr) {
          int nbr = lattice.neighbors[active_site][ibr];

          // If the active site and neighbor have the same spin,
          // && current neighbor not already in cluster
          // && detailed balance is satisified ( rand < pCluster )
          if ((lattice.spin[active_site] == lattice.spin[nbr]) &&
              (std::find(cluster.begin(), cluster.end(), nbr) ==
               cluster.end()) &&
              (rand_01(generator) <
               pCluster)) { // Then add the neighbor to the cluster and pocket
            cluster.push_back(nbr);
            pocket.push_back(nbr);
          }
        }

      } while (pocket.size() != 0);

      // Now, flip the spins of the entire cluster
      for (size_t i = 0; i < cluster.size(); ++i)
        lattice.spin[cluster[i]] *= -1;

      n_flipped_spins += cluster.size();
    } while (n_flipped_spins < n_spins_to_flip);
  }
};

//----------------------------------------------------------------------------//

class particle_swap : public action {

  double beta;
  particle_swap(double T) : beta(1.0 / T){};

  void step(ising::nodes &lattice, std::default_random_engine &generator) {
    bool do_swap = true;
    // For swap, select two random sites (orig, dest) of different spin
    int orig_site = rand_lattice_site(generator);
    int dest_site = rand_lattice_site(generator);
    if (orig_site == dest_site ||
        (lattice.spin[orig_site] == lattice.spin[dest_site])) {
      do_swap = false;
    }

    // Do not allow swaps between sites when one is frozen
    if (lattice.frozen[orig_site] || lattice.frozen[dest_site])
      do_swap = false;

    // If swap is allowed
    if (do_swap) {
      // Calculate the energy before the swap
      const double E0 = lattice.calcE_for_two_sites(orig_site, dest_site);

      // Create trial swap by swapping the spins. Recalc Energy
      ising::swap_spins(lattice.spin, orig_site, dest_site);
      const double E1 = lattice.calcE_for_two_sites(orig_site, dest_site);

      const double dE = E1 - E0;

      // Monte-Carlo acceptance criteria
      double p = std::exp(-beta * dE);
      double r = rand_01(generator);
      // Test if the swap was REJECTED. If so, re-swap spins to undo swap
      if (r > p) {
        ising::swap_spins(lattice.spin, orig_site, dest_site);
      }
    }
  }
};

//----------------------------------------------------------------------------//

struct monte_carlo {

  double M_av;
  double numNeighbor_av;
  double numVertNeighbor_av;

  std::default_random_engine generator;
  std::shared_ptr<ising::nodes> lattice;

  std::vector<std::unique_ptr<action>> actions;

  int do_n_steps(int steps, const double T) {
    // (Re-)initialize the statistics
    int n_av(0);
    M_av = 0;
    numNeighbor_av = 0;
    numVertNeighbor_av = 0;

    int noccupied = lattice.calculate_noccupied();

    // Perform the requested number of steps
    for (size_t istep = 0; istep < steps; ++istep) {

      for (auto &mc_move : actions) {
        mc_move->step(lattice, generator);
      };

      ising::collect_stats(lattice, n_av, M_av, numNeighbor_av,
                           numVertNeighbor_av);
    }

    return n_av;
  }
};

////////////////////////////////////////////////////////////////////////////////

} // namespace ising
