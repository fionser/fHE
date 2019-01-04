#pragma once
#include "yell/params.hpp"
#include <type_traits>
namespace yell {
#define YELL_BINARY_OPERATOR(Op, Name) \
template<size_t degree> \
poly<degree> Op(poly<degree> const& op0, \
                poly<degree> const& op1) {\
  assert(op0.moduli_count() == op1.moduli_count()); \
  size_t nmoduli = std::min(op0.moduli_count(), \
                            op1.moduli_count()); \
  ops::Name op; \
  poly<degree> rop(nmoduli); \
  for (size_t cm = 0; cm < nmoduli; ++cm) { \
    auto dst = rop.ptr_at(cm); \
    auto op0_ptr = op0.cptr_at(cm); \
    auto op1_ptr = op1.cptr_at(cm); \
    for (size_t i = 0; i < degree; ++i) \
      *dst++ = op(*op0_ptr++, *op1_ptr++, cm);\
  } \
  return rop; \
}

#define YELL_SELF_BINARY_OPERATOR(Op, Name) \
template<size_t degree> \
poly<degree>& Op(poly<degree> &op0, \
                 poly<degree> const& op1) {\
  assert(op0.moduli_count() == op1.moduli_count()); \
  size_t nmoduli = std::min(op0.moduli_count(), \
                            op1.moduli_count()); \
  ops::Name op; \
  for (size_t cm = 0; cm < nmoduli; ++cm) { \
    auto dst = op0.ptr_at(cm); \
    auto op1_ptr = op1.cptr_at(cm); \
    for (size_t i = 0; i < degree; ++i) \
      op.compute(*dst++, *op1_ptr++, cm);\
  } \
  return op0; \
}

namespace ops
{

/* a -= (a >= p) ? p : 0 */
template<typename T>
static inline void mod_correct(T &a, const T p)
{
  using ST = typename std::make_signed<T>::type;
  a -= (p & static_cast<T>(-static_cast<ST>(a >= p)));
}

struct addmod {
  using T = typename params::value_type;

  T operator()(T x, T y, size_t cm) const {
    compute(x, y, cm);
    return x;
  }

  inline void compute(T &x, T y, size_t cm) const {
    auto const p = params::P[cm];
    assert(x < p && y < p);
    x += y;
    mod_correct(x, p);
  }
};

struct submod {
  using T = typename params::value_type;
  inline T operator()(T x, T y, size_t cm) const {
    compute(x, y, cm);
    return x;
  }

  inline void compute(T &x, T y, size_t cm) const {
    auto const p = params::P[cm];
    assert(x < p && y < p);
    x += (p - y);
    mod_correct(x, p);
  }
};

/*
   Barret Reduction.
   Reduce (*x) % params::P[cm].

   X * r = X * (4*2^64 + pm) = 4*X * 2^64 + X * pm is a 196-bit value.
   represent as two parts X := X1 * 2^64 + X0; then
   X * pm = X1 * pm * 2^64 + X0 * pm is a 196-bit value. We only care
   the highest 64-bit, i.e., (X1 * pm) >> 64;

   Final value is (4 * X * 2^64 + X * pm) >> 128. We do not have 192-bit word.
   So we compute the highest 64-bit as
      (4 * X + (X * pm) >> 64) >> 64
   -> (4*X + ((X >> 64) * pm)) >> 64;
*/
void barret_reduction(params::gt_value_type *x, size_t cm);

yell::params::value_type shoupify(yell::params::value_type x, size_t cm);

struct mulmod {
  using T = typename params::value_type;
  using gt_value_type = typename params::gt_value_type;
  inline T operator()(T x, T y, size_t cm) const 
  {
    assert(cm < params::kMaxNbModuli);
    gt_value_type res = (gt_value_type) x * y;
    barret_reduction(&res, cm);
    return (T) res;
  }

  inline void compute(T &x, T y, size_t cm) const {
    assert(cm < params::kMaxNbModuli);
    gt_value_type res = (gt_value_type) x * y;
    barret_reduction(&res, cm);
    x = (T) res;
  }
};

struct mulmod_shoup {
  using T = typename params::value_type;
  inline T operator()(T x, T y, T yprime, size_t cm) const
  {
    using gt_value_type = typename params::gt_value_type;
    auto const p = params::P[cm];
    T q = (((gt_value_type) x * yprime) >> params::kModulusRepresentationBitsize) * p;
    T res = (T) (x * y - q);
    mod_correct(res, p);
    return res;
  }

  inline void compute(T &x, T y, T yprime, size_t cm) const
  {
    using gt_value_type = typename params::gt_value_type;
    auto const p = params::P[cm];
    T q = (((gt_value_type) x * yprime) >> params::kModulusRepresentationBitsize) * p;
    x = x * y - q;
    mod_correct(x, p);
  }
};

struct muladd {
  using T = typename params::value_type;
  using gt_value_type = typename params::gt_value_type;
  inline T operator()(T rop, T x, T y, size_t cm) const {
    compute(rop, x, y, cm);
    return rop;
  }

  inline void compute(T &rop, T x, T y, size_t cm) const {
    gt_value_type res = (gt_value_type) x * y;
    res += rop;
    barret_reduction(&res, cm);
    rop = (T) res;
  }
};
} // namespace ops
} // namespace yell
