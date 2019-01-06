#pragma once
#include "yell/poly.hpp"
#include <cstdint>
#include <cstddef>
namespace fHE {
struct context {
  static constexpr size_t degree = 8192; 
  static constexpr size_t nr_ctxt_moduli = 4;
  static constexpr size_t nr_sp_primes = nr_ctxt_moduli;
  static constexpr double sigma = 3.2;
  static constexpr double encoder_scale = (double) (1ULL << 40);
  static constexpr double plain_scale = (double) (1ULL << 20);
  static constexpr size_t HWT = 32;// hamming weight for the secret key

  static size_t index_sp_prime(size_t k) { 
    assert(k < nr_sp_primes); 
    return k + nr_ctxt_moduli; 
  }

  using poly_t = yell::poly<degree>;
  using gauss_struct = yell::gaussian<uint16_t, yell::params::value_type, 2>;
  using gauss_t = yell::FastGaussianNoise<uint16_t, yell::params::value_type, 2>;
  static gauss_t fg_prng;
};
} // namespace fHE
