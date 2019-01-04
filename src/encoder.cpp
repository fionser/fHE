#include "fHE/encoder.hpp"
#include <yell/utils/math.hpp>
#include <cassert>
#include <gmp.h>
namespace fHE {
Encoder::Encoder()
{
  const size_t m = degree << 1u;
  uint64_t pos = 1;
  for (size_t i = 0; i < nr_slots; i++) {
    // Position in normal bit order
    uint64_t index1 = (pos - 1) >> 1;
    uint64_t index2 = (m - pos - 1) >> 1;

    // Set the bit-reversed locations
    matrix_reps_index_map_.at(i) = yell::math::reverse_bits(index1, logn);
    matrix_reps_index_map_.at(nr_slots | i) = yell::math::reverse_bits(index2, logn);
    assert(matrix_reps_index_map_[i] < degree);
    assert(matrix_reps_index_map_[nr_slots | i] < degree);

    // Next primitive root
    pos *= generator;
    pos &= (m - 1); // mod m
  }

  const double PI_ = 3.14159265358979323846;
  std::complex<double> psi{ std::cos((2 * PI_) / (double)(m)), 
    std::sin((2 * PI_) / (double)(m)) };

  for (size_t i = 0; i < degree; i++) {
    roots_[i] = std::pow(psi, (double) (yell::math::reverse_bits(i, logn)));
    inv_roots_[i] = 1.0 / roots_[i];
  }
}

Encoder::~Encoder() {}

void Encoder::FFT_inplace(std::vector<std::complex<double>> &res) const {
  size_t tt = degree;
  for (size_t i = 0; i < logn; i++) {
    size_t mm = 1UL << i;
    tt >>= 1;

    for (size_t j = 0; j < mm; j++) {
      size_t j1 = 2 * j * tt;
      size_t j2 = j1 + tt - 1;
      auto s = roots_[mm + j];

      for (size_t k = j1; k < j2 + 1; k++) {
        assert(k < degree && k + tt < degree);
        auto u = res[k];
        auto v = res[k + tt] * s;
        res[k] = u + v;
        res[k + tt] = u - v;
      }
    }
  } // FFT
}

void Encoder::iFFT_inplace(std::vector<std::complex<double>> &res) const {
  size_t tt = 1;
  for (size_t i = 0; i < logn; i++) {
    size_t mm = 1ULL << (logn - i);
    size_t k_start = 0;
    size_t h = mm / 2;

    for (size_t j = 0; j < h; j++) {
      size_t k_end = k_start + tt;
      auto s = inv_roots_.at(h + j);

      for (size_t k = k_start; k < k_end; k++) {
        auto u = res[k];
        auto v = res[k + tt];
        res[k] = u + v;
        res[k + tt] = (u - v) * s;
      }

      k_start += 2 * tt;
    }
    tt *= 2;
  } // iFFT
}
std::vector<std::complex<double>> Encoder::apply_ifft(
  std::vector<std::complex<double>> const& values) const 
{
  const size_t input_sizes = std::min(values.size(), nr_slots);
  if (input_sizes != values.size()) {
    std::cerr << "WARN: Too many values to encode."
      << "Just handle the first " << nr_slots << " values.\n";
  }
  std::vector<std::complex<double>> conj_values(degree, 0.);
  for (size_t i = 0; i < input_sizes; ++i) {
    conj_values[matrix_reps_index_map_[i]] = values[i];
    conj_values[matrix_reps_index_map_[i + nr_slots]] = std::conj(values[i]);
  }
  iFFT_inplace(conj_values);
  return conj_values;
}

std::vector<std::complex<double>> Encoder::apply_ifft(
  std::vector<double> const& values) const 
{
  const size_t input_sizes = std::min(values.size(), nr_slots);
  if (input_sizes != values.size()) {
    std::cerr << "WARN: Too many values to encode."
              << "Just handle the first " << nr_slots << " values.\n";
  }
  std::vector<std::complex<double>> conj_values(degree, 0.);
  for (size_t i = 0; i < input_sizes; ++i) {
    conj_values[matrix_reps_index_map_[i]] = values[i];
    conj_values[matrix_reps_index_map_[i + nr_slots]] = values[i];
  }
  iFFT_inplace(conj_values);
  return conj_values;
}

void Encoder::encode(context::poly_t *msg,
                     std::vector<std::complex<double>> const& values,
                     const double scale) const 
{
  if (!msg) return;
  std::vector<std::complex<double>> conj_values{apply_ifft(values)};
  assert(conj_values.size() == context::degree);
  const double delta = scale / degree;
  double max_value = 1e-30;
  // TODO: avoid using mpz_t
  mpz_t coeff;
  mpz_init(coeff);
  const size_t nmoduli = msg->moduli_count();
  for (size_t cm = 0; cm < nmoduli; ++cm) {
    const auto P = (yell::params::signed_type) yell::params::P[cm];
    auto msg_ptr = msg->ptr_at(cm);
    for(size_t i = 0; i < degree; ++i) {
      double scaled = conj_values[i].real() * delta;
      mpz_set_d(coeff, scaled);
      mpz_mod_ui(coeff, coeff, P);
      max_value = std::max(max_value, std::abs(scaled));
      // i.e., msg(cm, i) = round<Int>(scaled) \bmod P[cm]
      auto t = (yell::params::signed_type) mpz_get_ui(coeff);
      assert(std::abs(t) <= P);
      *msg_ptr++ = (yell::params::value_type) (t < 0 ? P + t : t);
    } 
  }
  mpz_clear(coeff);
  if (max_value > 0)
    assert((int) std::log2(max_value) < nmoduli * yell::params::kModulusBitsize);
}

void Encoder::encode(context::poly_t *msg,
                     std::vector<double> const& values,
                     const double scale) const 
{
  std::vector<std::complex<double>> com_values(values.size(), 0.);
  std::transform(values.cbegin(), values.cend(), com_values.begin(),
                 [](double v) { return std::complex<double>(v); });
  encode(msg, com_values, scale);
}

void Encoder::decode(std::vector<std::complex<double>> *rop, 
                     context::poly_t const& msg,
                     double scale_) const 
{
  const size_t nmoduli = msg.moduli_count();
  assert(scale_ > 0 && (int) std::log2(scale_) < nmoduli * yell::params::kModulusBitsize);
  if (!rop) return;
  // TODO: avoid to use mpz
  mpz_t const& Q = yell::gmp::init_table(nmoduli)->moduli_product_;
  mpz_t Q2;
  mpz_init_set(Q2, Q);
  mpz_tdiv_q_ui(Q2, Q2, 2);
  auto coeffs_over_Q = msg.poly2mpz();
  const double inv_prec = 1. / scale_;
  std::vector<std::complex<double>> res(degree, 0.);
  for (size_t i = 0; i < degree; ++i) {
    if (mpz_cmp(coeffs_over_Q[i], Q2) >= 0)
      mpz_sub(coeffs_over_Q[i], coeffs_over_Q[i], Q);
    res[i].real(mpz_get_d(coeffs_over_Q[i]) * inv_prec);
    mpz_clear(coeffs_over_Q[i]);
  }

  mpz_clear(Q2);
  FFT_inplace(res);
  rop->resize(nr_slots);
  std::transform(matrix_reps_index_map_.cbegin(), 
                 matrix_reps_index_map_.cbegin() + nr_slots,
                 rop->begin(),
                 [&res](uint32_t i) { return res[i]; });
}

void Encoder::decode(std::vector<double> *rop, 
                     context::poly_t const &msg,
                     double scale_) const {
  if (!rop) return;
  std::vector<std::complex<double>> tmp;
  decode(&tmp, msg, scale_);
  rop->resize(tmp.size());
  std::transform(tmp.cbegin(), tmp.cend(), rop->begin(), 
                 [](std::complex<double> const& v) { return v.real(); });
}

} // namespace fHE
