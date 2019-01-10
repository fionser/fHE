#include "fHE/keys.hpp"
#include "fHE/mult_key.hpp"
#include "fHE/rot_key.hpp"
#include "fHE/encryptor.hpp"
#include "fHE/decryptor.hpp"
#include "fHE/cipher.hpp"
#include "fHE/encoder.hpp"
#include "fHE/evaluator.hpp"
#include "yell/utils/timer.hpp"
int main() {
  fHE::SK sk;
  fHE::PK pk(sk);
  fHE::MultKey mkey(sk);
  fHE::RotKey rkey(sk, 1, fHE::RotKey::Direction::RIGHT);
  fHE::context::poly_t A(fHE::context::nr_ctxt_moduli, yell::uniform{});
  auto cpy(A);
  fHE::Cipher ctx;
  auto encoder = std::make_shared<fHE::Encoder>();
  fHE::Encryptor encryptor(encoder);
  fHE::Decryptor decryptor(encoder);
  fHE::Evaluator evaluator;
  std::vector<double> values(fHE::context::degree / 2, 0.123456);

  std::cout << fHE::context::nr_ctxt_moduli << " moduli with " << yell::params::kModulusBitsize << " bits\n";
  std::cout << fHE::context::degree << " poly degree\n";
  encryptor.encrypt(&ctx, values, pk);
  double time = 0.;
  for (int i = 0; i < (1<<7); ++i) {
    AutoTimer timer(&time);
    //evaluator.multiply_inplace(&ctx, ctx, mkey);
    evaluator.rotate_slots_inplace(&ctx, rkey);
  }
  decryptor.decrypt(&values, ctx, sk);

  std::cout << "time " << time << " " << values[4] << "\n";
  return 0;
}
