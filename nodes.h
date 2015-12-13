#ifndef ISING_H
#define ISING_H

#include <cstdlib>
#include <iostream>
#include <vector>

namespace ising {

struct nodes {
    inline nodes();

    void init(int, int, int);
    void determine_connectivity();
    int find_site_index(int, int, int);
    void attempt_to_add(std::vector<int>&, int, int, int);

    int nx;
    int ny;
    int nz;
    int nsites;

    std::vector< std::vector<int> > neighbors;
    std::vector< std::vector<int> > diagNeighbors;
};

inline nodes::nodes()
: nx(0), ny(0), nz(0), nsites(0), neighbors(0), diagNeighbors(0)
{
}

} // namespace ising

#endif // ISING_H
