#include "math.h"
#include <algorithm>

namespace NMath {
    // https://en.algorithmica.org/hpc/algorithms/gcd/
    int64_t gcd(int64_t a, int64_t b) {
        if (a < 0) {
            a = -a;
        }
        if (b < 0) {
            b = -b;
        }
        if (a == 0) return b;
        if (b == 0) return a;

        int64_t az = __builtin_ctz(a);
        int64_t bz = __builtin_ctz(b);
        int64_t shift = std::min(az, bz);
        b >>= bz;
        
        while (a != 0) {
            a >>= az;
            int64_t diff = b - a;
            az = __builtin_ctz(diff);
            b = std::min(a, b);
            a = std::abs(diff);
        }
        
        return b << shift;
    }
}