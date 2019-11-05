#include <iostream>
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
#include "fHE/prime_bundles.hpp"
#include "fHE/special_prime_chain.hpp"
#include "fHE/context.hpp"
#include <unordered_set>
#include <gmp.h>

namespace fHE { struct context; }

void run(int n_spcl_primes) {
    size_t K = 50; // moduli size each
    int nrmal_bits = fHE::context::max_moduli_bits - K * n_spcl_primes;
    auto primeBundles = std::make_shared<fHE::PrimeBundles>((nrmal_bits + K - 1) / K, n_spcl_primes);

    const int n_nrml_primes = primeBundles->max_n_nrml_primes();
    if (n_spcl_primes > n_nrml_primes)
        return;
    std::cout << n_nrml_primes << " + " << n_spcl_primes << "\n";
    auto spChain = std::make_shared<fHE::SpecialPrimeChain>(n_nrml_primes, n_spcl_primes);

    fHE::SK sk(n_nrml_primes, n_spcl_primes);
    fHE::PK pk(sk);
    fHE::MultKey mkey(sk, primeBundles, spChain);

    fHE::RotKeySet<fHE::RotKey> rotkeys(sk, primeBundles, spChain);

    auto encoder = std::make_shared<fHE::Encoder>();
    fHE::Encryptor encryptor(encoder);
    fHE::Decryptor decryptor(encoder);
    fHE::Evaluator evaluator(encoder, primeBundles, spChain);

    std::vector<double> vals(fHE::context::degree/2);
    for (int i = 0; i < vals.size(); i++)
        vals[i] = 10. * (i + 1);

    fHE::Cipher rlwe(n_nrml_primes);
    encryptor.encrypt(&rlwe, vals, pk);

    double time = 0.;
    constexpr long NT = 100;
    yell::ntt<fHE::context::degree>::n_call_forward = 0;
    yell::ntt<fHE::context::degree>::n_call_backward = 0;
    for (long i = 0; i < NT; ++i) {

        AutoTimer timer(&time);
        evaluator.rotate_slots(&rlwe, 1, rotkeys);
    }
    std::cout << time / NT << "ms\n";
    std::cout << yell::ntt<fHE::context::degree>::n_call_forward << " NTT\n";
    std::cout << yell::ntt<fHE::context::degree>::n_call_backward << " invNTT\n";

    // if (n_spcl_primes == 1) {
        decryptor.decrypt(&vals, rlwe, sk);
        for (long i = 0; i < 8; ++i)
            std::cout << vals[i] << " ";
        std::cout << "\n";
    // }
}

int main(int argc, char *argv[]) {
    int n = argc < 1 ? 1 : std::atoi(argv[1]);
    // std::cout << fHE::context::max_moduli_bits << " bits\n";
    // for (long n = 13; n >= 1; --n)
    run(n);
    return 0;
}
