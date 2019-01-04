#pragma once
#include "fHE/context.hpp"
namespace fHE {
struct Cipher {
private:
  double scale_ = 1.;
public:
  context::poly_t bx;
  context::poly_t ax;

  explicit Cipher(size_t L = context::nr_ctxt_moduli) 
    : bx(L), ax(L) { assert(L < yell::params::kMaxNbModuli); }

  Cipher(const Cipher &oth) 
    : bx(oth.bx),
      ax(oth.ax),
      scale_(oth.scale_) 
  { }

  Cipher(Cipher &&oth) 
    : Cipher(0)
  { 
    std::swap(bx, oth.bx);
    std::swap(ax, oth.ax);
    scale_ = oth.scale_;
  }

  ~Cipher() { }

  Cipher& operator=(Cipher oth) = delete;

  double scale() const { return scale_; }

  size_t moduli_count() const { return bx.moduli_count(); }

  void scale(double s) { 
    assert(s > 0. && std::log2(s) < moduli_count() * yell::params::kModulusBitsize); 
    scale_ = s; 
  }
};

} // namespace fHE
