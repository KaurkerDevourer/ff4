#include "../benchmarking.h"
#include <vector>
#include <immintrin.h>

uint64_t sum_simple(const std::vector<uint32_t>& a) {
    uint64_t sum = 0;
    for (int i = 0; i < a.size(); i++) {
        sum += a[i];
    }
    return sum;
}

uint64_t sum_parallel(const std::vector<uint32_t>& a) {
    uint64_t sum = 0;
    for (int i = 0; i < a.size(); i += 16) {
        sum += a[i] + a[i + 1] + a[i + 2] + a[i + 3];
        sum += a[i + 4] + a[i + 5] + a[i + 6] + a[i + 7];
        sum += a[i + 8] + a[i + 9] + a[i + 10] + a[i + 11];
        sum += a[i + 12] + a[i + 13] + a[i + 14] + a[i + 15];
        __builtin_prefetch(&a[i + 16]);
    }
    return sum;
}

uint64_t sum_simd(const std::vector<uint32_t>& a) {
    uint64_t sum = 0;
    for (int i = 0; i < a.size(); i += 16) {
        sum += _mm512_reduce_add_epi32(_mm512_loadu_epi32(&a[i]));
    }
    return sum;
}

uint64_t sum_simd2(const std::vector<uint32_t>& a) {
    uint64_t sum = 0;
    for (int i = 0; i < a.size(); i += 16) {
        sum += _mm512_reduce_add_epi32(_mm512_loadu_epi32(&a[i]));
        __builtin_prefetch(&a[i + 16]);
    }
    return sum;
}

int main() {
    std::vector<uint32_t> a(65536); // 2 ^ 16
    for (int i = 0; i < 65536; i++) {
        a[i] = i;
    }
    FakeBenchmark(sum_simple, 10, "sum_simple", a);
    FakeBenchmark(sum_parallel, 10, "sum_parallel", a);
    FakeBenchmark(sum_simd, 10, "sum_simd", a);
    FakeBenchmark(sum_simd2, 10, "sum_simd2", a);
}
