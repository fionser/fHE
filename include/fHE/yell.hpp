#pragma once
#include "fHE/context.hpp"
#include <yell/ops.hpp>
namespace fHE {
void apply_galois(context::poly_t *, size_t galois);

void apply_galois_non_ntt(context::poly_t *, size_t galois);

template <class Itr>
yell::params::value_type product_of(Itr begin, Itr end, size_t moduli_index)
{
  assert(moduli_index <= yell::params::kMaxNbModuli);
  yell::params::value_type product{1L};
  yell::ops::mulmod mulmod;
  while (begin != end) {
    auto i = *begin++;
    assert(i < yell::params::kMaxNbModuli);
    mulmod.compute(product, yell::params::P[i], moduli_index);
  }
  return product;
}

//! rop += op0 * op1
void lazy_muladd(yell::params::gt_value_type *rop,
                 const yell::params::value_type* op0,
                 const yell::params::value_type* op1,
                 const size_t n);
} // namespace fHE
