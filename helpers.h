#include <cstdlib>
#include <algorithm>
#include <random>

#include "nodes.h"

namespace ising{

    int read_command_line_int(char*);

    void swap_spins(int*, int&, int&);

    void print_cell(ising::nodes&, std::string);

    int calcE_for_two_connected_sites(ising::nodes&, int&, int&);

    void collect_stats(ising::nodes&, int&, double&, double&, double&);

    int generated_desired_occupancy(ising::nodes&,
                                    std::default_random_engine&,
                                    std::uniform_int_distribution<int>&,
                                    int&);

} // namespace ising
