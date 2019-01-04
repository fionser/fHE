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
      tmp[d] = src_ptr[yell::math::reverse_bits((index_raw - 1) >> 1, logn)];
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

