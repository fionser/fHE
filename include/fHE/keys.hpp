#pragma once
#include "fHE/context.hpp"
namespace fHE {
struct SK {
public:
  context::poly_t sx; // sx over normal moduli
  context::poly_t ext_sx; // sx over special moduli

  explicit SK(const size_t hwt = context::HWT) 
    : sx(context::nr_ctxt_moduli),
      ext_sx(context::nr_sp_primes)
  {
    constexpr size_t L = context::nr_ctxt_moduli;
    constexpr size_t K = context::nr_sp_primes; 
    constexpr size_t bytes = sizeof(yell::params::value_type) * context::degree;
    context::poly_t tmp(L + K, yell::hwt_dist(hwt));
    tmp.forward();
    for (size_t cm = 0; cm < L; ++cm)
      std::memcpy(sx.ptr_at(cm), tmp.cptr_at(cm), bytes);

    for (size_t cm = 0; cm < K; ++cm)
      std::memcpy(ext_sx.ptr_at(cm), tmp.cptr_at(L + cm), bytes);
    tmp.clear();
  }

  SK(SK const& sk) = delete;

  SK(SK && sk) = delete;

  SK& operator=(SK const& sk) = delete;

  ~SK() { }
};

struct PK {
  context::poly_t bx; 
  context::poly_t ax; 

  explicit PK(SK const& sk) 
    : bx(context::nr_ctxt_moduli, context::gauss_struct(&context::fg_prng)),
      ax(context::nr_ctxt_moduli, yell::uniform{})
  {
    //! bx = -(e + ax * sk.sx);
    bx.forward_lazy();
    bx.add_product_of(ax, sk.sx);
    bx.negate();
  }

  PK(PK const& pk) = delete;

  PK(PK && pk) = delete;

  PK& operator=(PK const& pk) = delete;

  ~PK() { }
};

} // namespace fHE

