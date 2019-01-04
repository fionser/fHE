#pragma once
#include "yell/params.hpp"
#include "yell/utils/safe_ptr.h"
#include <unordered_map>
namespace yell {
template <size_t degree_>
class ntt {
private:
  static constexpr size_t degree = degree_;
  using T = params::value_type;
  using gt_value_type = params::gt_value_type;

  ntt() {}
  ntt(ntt const& oth) = delete;
  ntt(ntt && oth) = delete;
  ntt& operator=(ntt const& oth) = delete;

  ~ntt() {
    for (auto kv : tables) delete kv.second;
  }
  static ntt instance_; // singleton

  struct ntt_precomputed {
    ntt_precomputed() {}

    ~ntt_precomputed() {}

    void init(size_t cm);
    //! precomputed tables
    //! phi^(2 * degree) = 1 \bmod prime[cm]
    //! invphi = phi^(-1) \bmod prime[cm]
    //! shoupXX are Shoup's trick to accelerate the multiplication reduction.
    //! We merge the last layer of inv_ntt with the n^-1 step.
    //! Values in these tables are stored in the bit-reverse order.
    T phiTbl[degree],
      shoupPhiTbl[degree],
      invphiTbl[degree],
      shoupInvphiTbl[degree],
      invDegree,
      shoupInvDegree,
      invphiInvDegree, //! n^{-1} * invphi, used in the last layer of inv_ntt
      shoupInvphiInvDegree;
  };

  sf::contention_free_shared_mutex<> guard;
  std::unordered_map<size_t, ntt_precomputed *> tables;

public:
  static const ntt_precomputed* init_table(size_t cm);
  /* initialize the first K ntt tables
    */
  static void init_ntt_tables(size_t K) {
    assert(K < params::kMaxNbModuli);
    for (size_t cm = 0; cm < K; ++cm) init_table(cm);
  }
  /* apply forward ntt over the specified moduli
   */
  static void forward(T *op, size_t cm);
  //! skip the final modulus correction step.
  //! Result range in [0, 4 * p)
  static void forward_lazy(T *op, size_t cm);
  /* apply backward invntt over the specified moduli
   */
  static void backward(T* op ,size_t cm);
};
} // namespace yell

#include "yell/details/ntt.hpp"
