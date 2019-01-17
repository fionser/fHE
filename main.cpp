#include "fHE/keys.hpp"
#include "fHE/mult_key.hpp"
#include "fHE/rot_key.hpp"
#include "fHE/encryptor.hpp"
#include "fHE/decryptor.hpp"
#include "fHE/cipher.hpp"
#include "fHE/encoder.hpp"
#include "fHE/evaluator.hpp"
#include "yell/utils/timer.hpp"
#include "yell/montgomery.hpp"
#include "fHE/base_converter.hpp"
#include "fHE/yell.hpp"
#include <unordered_set>


int main() {
  fHE::SK sk(4);
  fHE::PK pk(sk);
  fHE::MultKey mkey(sk);
  fHE::RotKeySet rkey(sk);
  fHE::Cipher ctx, ctx0;
  auto encoder = std::make_shared<fHE::Encoder>();
  fHE::Encryptor encryptor(encoder);
  fHE::Decryptor decryptor(encoder);
  fHE::Evaluator evaluator(encoder);
  std::vector<double> values(fHE::context::degree / 2);
  for (size_t i = 0; i < values.size(); ++i)
    values[i] = 1.0 + i;

  std::cout << fHE::context::nr_ctxt_moduli << " moduli with " 
            << yell::params::kModulusBitsize << " bits\n";
  std::cout << fHE::context::degree << " poly degree\n";
  encryptor.encrypt(&ctx, values, pk);
  evaluator.replicate(&ctx, 4, rkey);
  decryptor.decrypt(&values, ctx, sk);
  std::cout << values[4] << "\n";
  return 0;
}
