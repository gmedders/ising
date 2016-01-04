#include <cstdlib>

#include "nodes.h"

namespace ising{

    int read_command_line_int(char*);

    void swap_spins(int*, int&, int&);

    void print_cell(ising::nodes&, std::string);

    int calcE_for_two_connected_sites(ising::nodes&, int&, int&);

} // namespace ising
