#pragma once
#include "fHE/keys.hpp"
#include <vector>
namespace fHE {
struct MultKey {
private:
  using poly_t = context::poly_t;
  static constexpr size_t L = context::nr_ctxt_moduli;
  std::vector<poly_t> alpha; // L polynomials, each L moduli
  std::vector<poly_t> beta ; // L polynomials, each L moduli

public:
  explicit MultKey(SK const& sk);

  MultKey(MultKey const& oth) = delete;

  MultKey(MultKey && oth) = delete;

  MultKey& operator=(MultKey const& oth) = delete;

  ~MultKey();

  const poly_t* get_beta_at(size_t j) const {
    return &beta.at(j);
  }

  const poly_t* get_alpha_at(size_t j) const {
    return &alpha.at(j);
  }

  poly_t* get_beta_at(size_t j) {
    return &beta.at(j);
  }

  poly_t* get_alpha_at(size_t j) {
    return &alpha.at(j);
  }
};
} // namespace fHE

