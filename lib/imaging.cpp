#include "imaging.h"

namespace ising {

int pos_mod(int &a, int &b) { return ((a % b + b) % b); }

} // namespace ising
