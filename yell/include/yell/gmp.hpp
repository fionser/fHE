#pragma once
#include "yell/utils/safe_ptr.h"
#include "yell/params.hpp"
#include <gmp.h>
#include <array>
#include <unordered_map>

namespace yell {
template <size_t degree> class poly; // forward declaration
class gmp {
private:
  using T = params::value_type;

  gmp() {}
  gmp(gmp const& oth) = delete;
  gmp& operator=(gmp const& oth) = delete;

  ~gmp() {
    for (auto kv : tables) delete kv.second;
  }
  static gmp instance_; // singleton

  struct gmp_precomputed {
    explicit gmp_precomputed(size_t nmoduli);
    ~gmp_precomputed();
    mpz_t moduli_product_;
    mpz_t modulus_shoup;
    size_t bits_in_moduli_product;
    size_t bits_in_modulus_shoup;
    size_t shift_modulus_shoup;
    mpz_t *lifting_integers;
    const size_t nmoduli_;
  };

  sf::contention_free_shared_mutex<> guard;
  std::unordered_map<size_t, gmp_precomputed *> tables;

public:
  static const gmp_precomputed* init_table(size_t nmoduli);

  template <size_t degree_>
  static std::array<mpz_t, degree_> poly2mpz(poly<degree_> const& op) {
    const size_t nmoduli_ = op.moduli_count();
    auto tbl = init_table(nmoduli_);
    std::array<mpz_t, degree_> rop;
    for (size_t i = 0; i < degree_; ++i)
      mpz_init2(rop[i], tbl->shift_modulus_shoup - 1);

    mpz_t tmp;
    mpz_init2(tmp, tbl->shift_modulus_shoup - 1 + tbl->bits_in_modulus_shoup);
    // Loop on each coefficient
    for (size_t i = 0; i < degree_; i++) {
      mpz_set_ui(rop[i], 0);
      for (size_t j = 0; j < nmoduli_; j++) {
        if (op(j, i) != 0)
          mpz_addmul_ui(rop[i], tbl->lifting_integers[j], op(j, i));
      }

      // Modular reduction using Shoup
      mpz_mul(tmp, rop[i], tbl->modulus_shoup);
      mpz_tdiv_q_2exp(tmp, tmp, tbl->shift_modulus_shoup);
      mpz_submul(rop[i], tmp, tbl->moduli_product_);
      if (mpz_cmp(rop[i], tbl->moduli_product_) >= 0) {
        mpz_sub(rop[i], rop[i], tbl->moduli_product_);
      }
    }
    // Clean
    mpz_clear(tmp);
    return rop;
  }
}; 
} // namespace yell
