#pragma once
#include <vector>
#include <cstdint>
#include <iostream>

namespace FF4 {
    namespace NUtils {
        class Term : public std::vector<uint64_t> {
                void Normalize();
            public:
                using TBase = std::vector<uint64_t>;

                inline Term() = default;

                Term(std::initializer_list<uint64_t> il);

                bool IsDivisibleBy(const Term&) const noexcept;
                uint64_t TotalDegree() const noexcept;
                Term& operator*=(const Term&) noexcept;
                friend Term operator*(Term, const Term&) noexcept;
                Term& operator/=(const Term&) noexcept;
                friend Term operator/(Term, const Term&) noexcept;

                friend Term gcd(const Term&, const Term&) noexcept;
                friend Term lcm(const Term&, const Term&) noexcept;

                friend std::ostream& operator<<(std::ostream&, const Term&) noexcept;
            private:
                uint64_t sum_ = 0;
        };
    }
}
