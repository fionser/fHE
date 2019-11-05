#pragma once
#include "fHE/context.hpp"
namespace fHE {
struct SK {
public:
  context::poly_t sx; // sx over normal primes and special primes
  int n_nrml_primes_;
  int n_spcl_primes_;

  explicit SK(int n_nrml_primes, int n_spcl_primes)
    : sx(n_nrml_primes + n_spcl_primes), n_nrml_primes_(n_nrml_primes), n_spcl_primes_(n_spcl_primes)
  {
    if (n_nrml_primes < 1 || n_spcl_primes < 1)
        throw std::invalid_argument("Secret key constructor: number of primes should >= 1");
    constexpr size_t bytes = sizeof(yell::params::value_type) * context::degree;
    sx.set(yell::ZO_dist{}); // P(s = 1) = P(s = -1) = 0.25, P(s = 0) = 0.5
    sx.forward();
  }

  SK(SK const& sk) = delete;

  SK(SK && sk) = delete;

  SK& operator=(SK const& sk) = delete;

  int n_nrml_primes() const { return n_nrml_primes_; }

  int n_spcl_primes() const { return n_spcl_primes_; }

  ~SK() { }
};

struct PK {
  context::poly_t bx; 
  context::poly_t ax; 

  explicit PK(SK const& sk) 
    : bx(sk.n_nrml_primes(), context::gauss_struct(&context::fg_prng)),
      ax(sk.n_nrml_primes(), yell::uniform{})
  {
    //! bx = -(e + ax * sk.sx);
    bx.forward_lazy();
    yell::ops::muladd muladd;
    using T = yell::params::value_type;
    for (long cm = 0; cm < sk.n_nrml_primes(); ++cm) {
        auto bx_ptr = bx.ptr_at(cm);
        auto ax_ptr = ax.cptr_at(cm);
        auto sx_ptr = sk.sx.cptr_at(cm);
        for (long d = 0; d < sk.sx.degree; ++d) {
            muladd.compute(*bx_ptr++, *ax_ptr++, *sx_ptr++, cm);
        }
    }
    bx.negate();
  }

  PK(PK const& pk) = delete;

  PK(PK && pk) = delete;

  PK& operator=(PK const& pk) = delete;

  ~PK() { }
};

} // namespace fHE

