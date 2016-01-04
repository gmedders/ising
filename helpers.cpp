#include "helpers.h"

#include <sstream>

namespace ising{

//----------------------------------------------------------------------------//

int read_command_line_int(char* argv)
{
    int val;
    std::istringstream iss(argv);
    iss >> val; 
    if (!iss || !iss.eof()) {
        std::cerr << "could not convert '" << argv
            << "' to int" << std::endl;
        return EXIT_FAILURE;
    }
    return val;
}

//----------------------------------------------------------------------------//

void swap_spins(int* spin, int& orig, int& dest)
{
    int tmp = spin[dest];
    spin[dest] = spin[orig];
    spin[orig] = tmp;
}

//----------------------------------------------------------------------------//

void print_cell(ising::nodes& lattice, std::string comment)
{    
    std::cerr << comment << std::endl;
    int* spin = lattice.spin;
    bool* frozen = lattice.frozen;
    for(int iz = 0; iz < lattice.nz; ++iz){
        std::cerr << " <<< Layer " << iz << " >>>" << std::endl;
        for(int ix = 0; ix < lattice.nx; ++ix){
            for(int iy = 0; iy < lattice.ny; ++iy){
                int site = lattice.find_site_index(ix, iy, iz);
                char front = ' ';
                char back  = ' ';
                if(frozen[site]){
                    front = '[';
                    back  = ']';
                }
                if(spin[site] < 0)
                    std::cerr << front << '-' << back << ' ';
                else if(spin[site] == 0)
                    std::cerr << front << 'O' << back << ' ';
                else
                    std::cerr << front << '+' << back << ' ';
            }
            std::cerr << std::endl;
        }
    }
}

//----------------------------------------------------------------------------//

int calcE_for_two_connected_sites(ising::nodes& lattice,
                                  int& orig_site, int& dest_site)
{
    int* spin = lattice.spin;
    int E0_orig(0);
    for(size_t ibr = 0; ibr < lattice.neighbors[orig_site].size();++ibr){
        int nbr = lattice.neighbors[orig_site][ibr];
        if(nbr != dest_site){
            E0_orig -= spin[orig_site] * spin[nbr];
        }
    }

    int E0_dest(0);
    for(size_t ibr = 0; ibr < lattice.neighbors[dest_site].size();++ibr){
        int nbr = lattice.neighbors[dest_site][ibr];
        if(nbr != orig_site){
            E0_dest -= spin[dest_site] * spin[nbr];
        }
    }

    return E0_orig + E0_dest;
}

} // namespace ising
