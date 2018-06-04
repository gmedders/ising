#ifndef HELPERS_H
#define HELPERS_H

#include <algorithm>
#include <cstdlib>

#include "nodes.h"

namespace ising {

int read_command_line_int(const char *);
double read_command_line_double(const char *);

void swap_spins(int *, int &, int &);

void print_cell(ising::nodes &, std::string);

void collect_stats(ising::nodes &, int &, double &, double &, double &);

} // namespace ising
#endif // HELPERS_H
