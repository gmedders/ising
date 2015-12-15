#include "nodes.h"

//#define PBC true
//#define NO_PBC true

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

#ifdef PBC
    std::cerr << "# Full PBC" << std::endl;
    std::cout << "# Full PBC" << std::endl;
#elif defined(NO_PBC)
    std::cerr << "# No PBC!!" << std::endl;
    std::cout << "# No PBC!!" << std::endl;
#else
    std::cerr << "# Slab geometry: PBC only in XY" << std::endl;
    std::cout << "# Slab geometry: PBC only in XY" << std::endl;
#endif

    determine_connectivity();
}

//----------------------------------------------------------------------------//

int nodes::find_site_index(int ix, int iy, int iz)
{
#ifdef PBC
    return pos_mod(iz, nz)*ny*nx + pos_mod(iy, ny)*nx + pos_mod(ix, nx);
#elif defined(NO_PBC)
    // within the system (that is not replicated in x,y,z)
    if(ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz)
	return iz*ny*nx + iy*nx + ix;
    else
	return -1;
#else
    if(iz >= 0 && iz < nz) // within the system (that is not replicated in z)
	return iz*ny*nx + pos_mod(iy, ny)*nx + pos_mod(ix, nx);
    else
	return -1;
#endif
}

//----------------------------------------------------------------------------//

void nodes::attempt_to_add(std::vector<int>& vint, std::vector<bool>& vbool,
                           int i, int j, int k, bool vert)
{
    int site = ising::nodes::find_site_index(i, j, k);
    if(site >= 0){
	vint.push_back(site);
        vbool.push_back(vert);
    }
}

//----------------------------------------------------------------------------//

void nodes::determine_connectivity()
{

    for(int iz = 0; iz < nz; ++iz){
        for(int iy = 0; iy < ny; ++iy){
            for(int ix = 0; ix < nx; ++ix){

              // First find the nearest neighbors
              {
	        std::vector<int> connections;
	        std::vector<bool> vertical;

	        // +/- x
	        attempt_to_add(connections, vertical, ix + 1, iy, iz, false);
	        attempt_to_add(connections, vertical, ix - 1, iy, iz, false);

	        // +/- y
	        attempt_to_add(connections, vertical, ix, iy + 1, iz, false);
	        attempt_to_add(connections, vertical, ix, iy - 1, iz, false);

	        // +/- z
	        // For 2D lattice, don't add +/- z
	        if (nz > 1){
	            attempt_to_add(connections, vertical, ix, iy, iz + 1, true);
	            attempt_to_add(connections, vertical, ix, iy, iz - 1, true);
	        }

	        neighbors.push_back(connections);
                neighborVertical.push_back(vertical);
              }

#if 0
                // Also construct list including NN and diagonal neighbors
                {
		    std::vector<int> connections;

		    // +/- x
                    for(int dx = -1; dx < 2; ++dx){
                    for(int dy = -1; dy < 2; ++dy){
                        // For 2D lattice, don't add +/- z
                        if(nz > 1){
                            for(int dz = -1; dz < 2; ++dz){
                                if(!(dx == 0 && dy == 0 && dz == 0)){
                                    attempt_to_add(connections, ix + dx,
                                                                iy + dy,
                                                                iz + dz);
                                }
                            }
                        }else{
                            int dz = 0;
                            if(!(dx == 0 && dy == 0 && dz == 0)){
                                attempt_to_add(connections, ix + dx,
                                                            iy + dy,
                                                            iz + dz);
                            }
                        }

		    }
                    }

		    diagNeighbors.push_back(connections);
                }
#endif

	    }
	}
    }

}

//----------------------------------------------------------------------------//

} // namespace ising
