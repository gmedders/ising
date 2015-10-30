#include "neighbors.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

int pos_mod(int& a, int& b){ return ( (a%b + b)%b ); }

//----------------------------------------------------------------------------//

void attempt_to_add(std::vector<int>& vint, int& i, int& j, int& k)
{
    int site = ising::find_site_index(i, j, k);
    if(site >= 0)
	vint.push_back(site);
}

//----------------------------------------------------------------------------//

}

////////////////////////////////////////////////////////////////////////////////

namespace ising {
   
//----------------------------------------------------------------------------//

int find_site_index(int& ix, int& iy, int& iz)
{
#ifdef PBC_Z
    return pos_mod(iz%nz)*ny*nx + pos_mod(iy%ny)*nx + pos_mod(ix%nx);
#else
    if(iz >= 0 && iz < nz) // within the system (that is not replicated in z)
	return iz*ny*nx + pos_mod(iy%ny)*nx + pos_mod(ix%nx);
    else
	return -1;
#endif
}

//----------------------------------------------------------------------------//

void nodes::determine_connectivity()
{

    for(size_t ix = 0; ix < nx; ++ix){
	for(size_t iy = 0; iy < ny; ++iy){
	    for(size_t iz = 0; iz < nz; ++iz){

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
