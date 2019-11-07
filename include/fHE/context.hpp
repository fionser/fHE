#pragma once
#include "yell/poly.hpp"
#include <cstdint>
#include <cstddef>
namespace fHE {
struct context {
  static constexpr size_t degree = 1L << 14;
  static constexpr size_t max_moduli_bits = 620;
  static constexpr double sigma = 3.2;
  static constexpr double encoder_scale = (double) (1ULL << 62);
  static constexpr double plain_scale = (double) (1ULL << 20);
  static constexpr size_t HWT = 64;// hamming weight for the secret key

  using poly_t = yell::poly<degree>;
  using gauss_struct = yell::gaussian<uint16_t, yell::params::value_type, 2>;
  using gauss_t = yell::FastGaussianNoise<uint16_t, yell::params::value_type, 2>;
  static gauss_t fg_prng;
};
} // namespace fHE
