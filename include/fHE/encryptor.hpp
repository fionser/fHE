#pragma once
#include "fHE/context.hpp"
#include <vector>
#include <complex>
#include <memory>
namespace fHE {
struct Cipher;
struct PK;
struct Encoder;

struct Encryptor {
  explicit Encryptor(std::shared_ptr<Encoder> encoder = nullptr);

  bool encrypt(Cipher *rop, context::poly_t const& msg, PK const& pk) const;

  bool encrypt(Cipher *rop, std::vector<double> const& values, PK const &pk) const;

  bool encrypt(Cipher *rop, std::vector<std::complex<double>> const& values, PK const &pk) const;

  std::shared_ptr<Encoder> encoder;
};

} // namespace fHE
