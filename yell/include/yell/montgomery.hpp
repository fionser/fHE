#pragma once
#include "yell/params.hpp"
#include "yell/utils/safe_ptr.h"
#include <unordered_map>
namespace yell {
class montgomery {
private:
  using T = params::value_type;
  using gt_value_type = params::gt_value_type;

  montgomery() {}
  montgomery(montgomery const& oth) = delete;
  montgomery(montgomery && oth) = delete;
  montgomery& operator=(montgomery const& oth) = delete;

  ~montgomery() {
    for (auto kv : tables) delete kv.second;
  }
  static montgomery instance_; // singleton

  struct mot_precomputed {
  public:
    explicit mot_precomputed(size_t cm);
    ~mot_precomputed() {}

    int w; //! R := 2^w
    size_t cm;
    T mu; //! -p^{-1} mod R
    T invR; //! R^{-1} mod p
    T shoupInvR;
    void to_montgomery(T *op, size_t degree);
    void reduce_from_montgomery(T *op, size_t degree);
  };

  sf::contention_free_shared_mutex<> guard;
  std::unordered_map<size_t, mot_precomputed *> tables;

public:
  static const mot_precomputed* init_table(size_t cm);
};
} // namespace yell

#include "yell/details/montgomery.hpp"
