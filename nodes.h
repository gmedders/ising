#ifndef ISING_H
#define ISING_H

#include <cstdlib>
#include <vector>

namespace ising {

struct nodes {
    inline nodes();

    void init(int, int, int);
    void determine_connectivity();
    int find_site_index(int&, int&, int&);

    int nx;
    int ny;
    int nz;
    int nsites;

    static std::vector< std::vector<int> > neighbors;
};

inline nodes::nodes()
: nx(0), ny(0), nz(0), nsites(0), neighbors(0)
{
}

void nodes::init(int my_nx, int my_ny, int my_nz)
{
    nx = my_nx;
    ny = my_ny;
    nz = my_nz;

    nsites = nx*ny*nz;

    if(nx < 1 || ny < 1 || nz < 1){
	std::cerr << "The number of nodes in each dimension must be at least 1"
	          << " you did: nx = " << nx << ", "
	          << " ny = " << ny << ", "
		  << " nz = " << nz << std::endl;
	exit(1);
    }

    determine_connectivity();
}

} // namespace ising

#endif // ISING_H
