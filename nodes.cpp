#include "nodes.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

int pos_mod(int& a, int& b){ return ( (a%b + b)%b ); }

//----------------------------------------------------------------------------//

}

////////////////////////////////////////////////////////////////////////////////

namespace ising {
   
//----------------------------------------------------------------------------//

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

//----------------------------------------------------------------------------//

int nodes::find_site_index(int ix, int iy, int iz)
{
#ifdef PBC_Z
    return pos_mod(iz, nz)*ny*nx + pos_mod(iy, ny)*nx + pos_mod(ix, nx);
#else
    if(iz >= 0 && iz < nz) // within the system (that is not replicated in z)
	return iz*ny*nx + pos_mod(iy, ny)*nx + pos_mod(ix, nx);
    else
	return -1;
#endif
}

//----------------------------------------------------------------------------//

void nodes::attempt_to_add(std::vector<int>& vint, int i, int j, int k)
{
    int site = ising::nodes::find_site_index(i, j, k);
    if(site >= 0)
	vint.push_back(site);
}

//----------------------------------------------------------------------------//

void nodes::determine_connectivity()
{

    for(int iz = 0; iz < nz; ++iz){
        for(int iy = 0; iy < ny; ++iy){
            for(int ix = 0; ix < nx; ++ix){

		std::vector<int> connections;

		// +/- x
		attempt_to_add(connections, ix + 1, iy, iz);
		attempt_to_add(connections, ix - 1, iy, iz);
		//connections.push_back( find_site_index( ix + 1, iy, iz) );
		//connections.push_back( find_site_index( ix - 1, iy, iz) );

		// +/- y
		attempt_to_add(connections, ix, iy + 1, iz);
		attempt_to_add(connections, ix, iy - 1, iz);
		//connections.push_back( find_site_index( ix, iy + 1, iz) );
		//connections.push_back( find_site_index( ix, iy - 1, iz) );

		// +/- z
		// For 2D lattice, don't add +/- z
		if (nz > 1){
		    attempt_to_add(connections, ix, iy, iz + 1);
		    attempt_to_add(connections, ix, iy, iz - 1);
		    //connections.push_back( find_site_index( ix, iy, iz + 1) );
		    //connections.push_back( find_site_index( ix, iy, iz - 1) );
		}

		neighbors.push_back(connections);

	    }
	}
    }

}

//----------------------------------------------------------------------------//

} // namespace ising
