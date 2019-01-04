#include "fHE/context.hpp"
namespace fHE {
context::gauss_t context::fg_prng(context::sigma, 128, 1L << 14);
}

