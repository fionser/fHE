#include "fHE/yell.hpp"
#include <yell/utils/math.hpp>
namespace fHE {
void apply_galois(context::poly_t *op, size_t galois)
{
  if (!op) return;
  assert((galois & 1) && galois < 2 * context::degree);
  if (galois == 1) 
    return;
  const size_t logn = op->logn;
  const size_t mod_mask = ((2 << logn) - 1);
  std::array<yell::params::value_type, context::degree> tmp;
  const size_t nmoduli = op->moduli_count();
  for (size_t cm = 0; cm < nmoduli; ++cm) {
    auto src_ptr = op->cptr_at(cm);
    for (size_t d = 0; d < context::degree; ++d) {
      uint64_t reversed = yell::math::reverse_bits(d, logn);
      uint64_t index_raw = galois * (2 * reversed + 1);
      index_raw &= mod_mask;
      auto dd = yell::math::reverse_bits(d, logn);
      tmp[d] = src_ptr[yell::math::reverse_bits((index_raw - 1) >> 1, logn)];
    }
    std::memcpy(op->ptr_at(cm), tmp.begin(), sizeof(tmp));
  }
}

//! f(x) -> f(x^g)
//! make sure that the polynomial is in the power-basis
void apply_galois_non_ntt(fHE::context::poly_t *op, size_t galois)
{
  if (!op) return;
  constexpr size_t degree = context::degree;
  assert((galois & 1) && galois < 2 * degree);
  if (galois == 1) 
    return;
  const size_t logn = op->logn;
  const size_t mod_mask = ((1 << logn) - 1);
  std::array<yell::params::value_type, degree> tmp;
  const size_t nmoduli = op->moduli_count();

  for (size_t cm = 0; cm < nmoduli; ++cm) {
    auto input = op->cptr_at(cm);
    for (size_t d = 0; d < degree; ++d) {
      auto raw_index = d * galois;
      auto index = raw_index & mod_mask;
      auto t = *input++;
      if ((raw_index >> logn) & 1)
        tmp[index] = t > 0 ? yell::params::P[cm] - t : 0;
      else
        tmp[index] = t;
    }
    std::memcpy(op->ptr_at(cm), tmp.begin(), sizeof(tmp));
  }
}

//! rop += op0 * op1
void lazy_muladd(yell::params::gt_value_type *rop,
                 const yell::params::value_type* op0,
                 const yell::params::value_type* op1,
                 const size_t n)
{
  if (!rop || !op0 || !op1)
    return;
  for (size_t i = 0; i != n; ++i, ++rop)
    *rop += (yell::params::gt_value_type) (*op0++) * (*op1++);
}

} // namespace fHE

