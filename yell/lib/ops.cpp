#include "yell/params.hpp"
#include <cassert>
namespace yell{ namespace ops {
void barret_reduction(params::gt_value_type *x, size_t cm)
{
  if (!x) return;
  // NOTE: x should less than 2^{126}.
  using T  = params::value_type;
  using gT = params::gt_value_type;
  using ST = params::signed_type;

  const T p = params::P[cm];
  const T pm = params::Pn[cm];
  constexpr size_t shift = params::kModulusRepresentationBitsize;
  constexpr size_t delta = shift - params::kModulusBitsize;
  gT q = ((gT) pm * ((*x) >> shift)) + ((*x) << delta) ;

  T r = (*x) - (q >> shift) * p;
  r -= (p & static_cast<T>(-static_cast<ST>(r >= p)));
  *x = (gT) r;
} 

yell::params::value_type shoupify(yell::params::value_type x, size_t cm)
{
  using gT = params::gt_value_type;
  return ((gT) x << params::kModulusRepresentationBitsize) / params::P[cm];
}
}} // namespace yell::ops
