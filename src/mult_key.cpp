#include "fHE/mult_key.hpp"
namespace fHE {
MultKey::MultKey(SK const& sk)
  : alpha(std::vector<poly_t>(L, poly_t(L, yell::uniform{}))),
     beta(std::vector<poly_t>(L, poly_t(L)))
{
  constexpr size_t degree = context::degree;
  constexpr size_t L = context::nr_ctxt_moduli;
  auto sx_square(sk.sx);
  sx_square *= sk.sx;
  for (size_t j = 0; j < L; ++j) {
    poly_t *beta_ = get_beta_at(j);
    const poly_t *alpha_ = get_alpha_at(j); 
    beta_->set(context::gauss_struct(&context::fg_prng));
    beta_->forward();
    beta_->add_product_of(*alpha_, sk.sx);
    beta_->negate(); // -e - alpha * sx 

    //! Mathematically, we set all moduli except the j-th of sx_square to zero.
    //! Then add it the beta[j]. 
    //! We can just add the j-th moduli of sx_square to beta[j].
    //! beta = sx^2 - e - alpha * sx
    yell::ops::addmod addmod;
    auto sx2_ptr = sx_square.cptr_at(j);
    auto beta_ptr = beta_->ptr_at(j);
    for (size_t d = 0; d < degree; ++d)
      addmod.compute(*beta_ptr++, *sx2_ptr++, j);
  }
  sx_square.clear(); //! clear from RAM
}

MultKey::~MultKey() {}
} // namespace fHE
