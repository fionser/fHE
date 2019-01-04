#include "fHE/keys.hpp"
#include "fHE/mult_key.hpp"
#include "fHE/rot_key.hpp"
#include "fHE/encryptor.hpp"
#include "fHE/decryptor.hpp"
#include "fHE/cipher.hpp"
#include "fHE/encoder.hpp"
#include "fHE/evaluator.hpp"
int main() {
  fHE::SK sk;
  fHE::PK pk(sk);
  fHE::MultKey mkey(sk);
  fHE::RotKey rkey(sk, 1, fHE::RotKey::Direction::LEFT);
  fHE::context::poly_t A(fHE::context::nr_ctxt_moduli, yell::uniform{});
  auto cpy(A);
  fHE::Cipher ctx;
  auto encoder = std::make_shared<fHE::Encoder>();
  fHE::Encryptor encryptor(encoder);
  fHE::Decryptor decryptor(encoder);
  fHE::Evaluator evaluator;
  std::vector<double> values{0.1, 0.2, 0.3, 0.4};

  encryptor.encrypt(&ctx, values, pk);
  evaluator.multiply_inplace(&ctx, ctx, mkey);
  evaluator.rotate_slots_inplace(&ctx, rkey);
  decryptor.decrypt(&values, ctx, sk);

  std::cout << values[0] << " " << values[1] << " " << values[4] << "\n";
  return 0;
}
