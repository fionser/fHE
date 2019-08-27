#pragma once
#include "fHE/context.hpp"
#include <vector>
#include <memory>
#include <unordered_map>
namespace fHE {
struct SK; // forward declaration

struct RotKey {
private:
  static constexpr size_t degree = context::degree;
  static constexpr size_t L = context::nr_ctxt_moduli;
  static constexpr size_t K = context::nr_sp_primes;

  /*
     Multiply P (i.e., the product of all special primes).
     [x]_{q1, q2, ..., qL} --> [P * x]_{q1, q2, ..., qL}
   */
  void lift_normal_moduli(context::poly_t *op) const;

public:
  enum class Direction {
    LEFT,
    RIGHT
  };
  const size_t galois;
  context::poly_t beta; //! L + K moduli
  context::poly_t alpha;//! L + K moduli

  explicit RotKey(SK const &sk, uint16_t offset, Direction dir);

  RotKey(RotKey const&oth) = delete;

  RotKey(RotKey &&oth) = delete;

  RotKey& operator=(RotKey const&oth) = delete;

  ~RotKey();
};

struct MixedRotKey {
private:
  static constexpr size_t degree = context::degree;
  static constexpr size_t L = context::nr_ctxt_moduli;
  using poly_t = context::poly_t;
  std::vector<poly_t> alpha; // L polynomials, each L + 1 moduli
  std::vector<poly_t> beta ; // L polynomials, each L + 1 moduli

public:
  enum class Direction {
    LEFT,
    RIGHT
  };

  const size_t galois;

  explicit MixedRotKey(SK const &sk, uint16_t offset, Direction dir);

  MixedRotKey(MixedRotKey const&oth) = delete;

  MixedRotKey(MixedRotKey &&oth) = delete;

  MixedRotKey& operator=(MixedRotKey const&oth) = delete;

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

  ~MixedRotKey();
};

template <class RotKeyType>
struct RotKeySet {
private:
  static constexpr size_t logn = yell::static_log2<context::degree / 2>::value;
  using rotkey_ptr = std::shared_ptr<RotKeyType>;
  std::unordered_map<size_t, rotkey_ptr> rot_keys;

public:
  //! create log(N) (left) rotation key so that we can rotate with
  //! any offset in [0, N).
  explicit RotKeySet(SK const& sk) {
    for (size_t i = 0; i < logn; ++i) {
      int idx = 1 << i;
      auto obj = std::make_shared<RotKeyType>(sk, idx, RotKeyType::Direction::LEFT);
      rot_keys.insert({idx, std::move(obj)});
    }
  }

  RotKeySet(RotKeySet const&oth) = delete;

  RotKeySet(RotKeySet &&oth) = delete;

  RotKeySet& operator=(RotKeySet const&oth) = delete;

  const RotKeyType* get(uint16_t offset) const {
      if (offset >= (1LL << logn))
          return nullptr;
      auto kv = rot_keys.find(offset);
      if (kv == rot_keys.end())
          return nullptr;
      return kv->second.get();
  }

  ~RotKeySet() {}
};

} // namespace fHE

