#pragma once
#include <immintrin.h>
//! multiply the packed 64 bit integers in A, B;
//! and store the high 64-bit result.
__m256i avx_mm256_mul64_hi(__m256i const& A, __m256i const& B);

__m256i avx_mm256_mul64_hi(__m256i const& A, 
                           __m256i const& Ahi, 
                           __m256i const& B);

//! multiply the packed 64 bit integers in A, B;
//! and store the low 64-bit result.
__m256i avx_mm256_mul64_lo(__m256i const& A, __m256i const& B);
__m256i avx_mm256_mul64_lo(__m256i const& A, 
                           __m256i const& Ahi, 
                           __m256i const& B);

