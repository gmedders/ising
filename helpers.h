#include <cstdlib>

#include "nodes.h"

namespace ising{

    int read_command_line_int(char*);

    void swap_spins(int*, int&, int&);

    void print_cell(ising::nodes&, std::string);

    int calcE_for_two_connected_sites(ising::nodes&, int&, int&);

    void collect_stats(ising::nodes&, int&, double&, double&, double&);

} // namespace ising
