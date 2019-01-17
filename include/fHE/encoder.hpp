#pragma once
#include "fHE/context.hpp"
#include <array>
#include <vector>
#include <complex>
namespace fHE {
struct Encoder {
private:
  static constexpr size_t degree = context::degree;
  std::array<uint32_t, degree> matrix_reps_index_map_;
  std::array<std::complex<double>, degree> roots_;
  std::array<std::complex<double>, degree> inv_roots_;

  const size_t nr_slots = degree >> 1u;
  const size_t logn = yell::static_log2<context::degree>::value;
  const size_t generator = 3;

  void FFT_inplace(std::vector<std::complex<double>> &res) const;

  void iFFT_inplace(std::vector<std::complex<double>> &res) const;

  std::vector<std::complex<double>> apply_ifft(
   std::vector<std::complex<double>> const& values) const;

  std::vector<std::complex<double>> apply_ifft(
    std::vector<double> const& values) const;
public:
  explicit Encoder();

  ~Encoder();

  void encode(context::poly_t *msg,
              std::vector<std::complex<double>> const& values,
              const double scale = context::encoder_scale) const;

  void encode(context::poly_t *msg,
              std::vector<double> const& values,
              const double scale = context::encoder_scale) const;

  void decode(std::vector<std::complex<double>> *rop,
              context::poly_t const& msg,
              double scale_) const;

  void decode(std::vector<double> *rop,
              context::poly_t const& msg,
              double scale_) const;
};
} // namespace fHE
