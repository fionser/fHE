#pragma once
#include "yell/ops.hpp"
#include <cmath>
namespace yell {
namespace math {
// compute a^{-1} mod prime
// i.e., a^{prime - 2}
template <typename T>
T inv_mod_prime(T a, size_t cm)
{
  yell::ops::mulmod mulmod;
  T ret{1};
  T y = yell::params::P[cm] - 2;
  while (y > 0) {
    if (y & 1)
      ret = mulmod(ret, a, cm);
    y >>= 1u;
    a = mulmod(a, a, cm);
  }
  return ret;
}

template <typename T, typename F>
inline T round(F a)
{
  return a >= 0 ? (T)(a + 0.5f) : (T) (a - 0.5f);
}

/* Reverse 32-bit */
uint32_t reverse_bits(uint32_t operand);
uint32_t reverse_bits(uint32_t operand, int32_t bit_count);
}} // namespace yell::math

