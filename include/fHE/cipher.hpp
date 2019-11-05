#pragma once
#include "fHE/context.hpp"
#include <memory>
namespace fHE {
struct Cipher {
private:
  double scale_ = 1.;
public:
  //! decryption logic: m = bx + ax * sx
  context::poly_t bx;
  context::poly_t ax;
  std::shared_ptr<context::poly_t> mul_aux = nullptr;

  explicit Cipher(size_t L) 
    : bx(L), ax(L) { 
        if (L == 0 || L >= yell::params::kMaxNbModuli)
            throw std::invalid_argument("Cipher:: invalid number of moduli");
    }

  Cipher(const Cipher &oth) 
    : scale_(oth.scale_),
      bx(oth.bx),
      ax(oth.ax),
      mul_aux(oth.mul_aux)
  { }

  Cipher(Cipher &&oth) 
    : scale_(oth.scale_),
      bx(std::move(oth.bx)),
      ax(std::move(oth.ax)),
      mul_aux(std::move(oth.mul_aux))
  { }

  ~Cipher() { }

  Cipher& operator=(Cipher oth) {
    std::swap(bx, oth.bx);
    std::swap(ax, oth.ax);
    std::swap(mul_aux, oth.mul_aux);
    scale_ = oth.scale_;
    return *this;
  }

  bool canonical() const { return mul_aux == nullptr; };

  double scale() const { return scale_; }

  size_t moduli_count() const { return bx.moduli_count(); }

  void scale(double s) { 
    assert(s > 0. && std::log2(s) < moduli_count() * yell::params::kModulusBitsize); 
    scale_ = s; 
  }
};

} // namespace fHE
