#pragma once
#include <vector>
#include <cstdint>
#include <iostream>

namespace FF4 {
    namespace NUtils {
        class TTerm : public std::vector<uint64_t> {
                void Normalize();
            public:
                using TBase = std::vector<uint64_t>;
                using size_type = typename TBase::size_type;

                inline TTerm() = default;

                TTerm(std::initializer_list<uint64_t> il);

                bool IsDivisibleBy(const TTerm&) const noexcept;
                uint64_t TotalDegree() const noexcept;
                TTerm& operator*=(const TTerm&) noexcept;
                friend TTerm operator*(TTerm, const TTerm&) noexcept;
                TTerm& operator/=(const TTerm&) noexcept;
                friend TTerm operator/(TTerm, const TTerm&) noexcept;

                friend TTerm gcd(const TTerm&, const TTerm&) noexcept;
                friend TTerm lcm(const TTerm&, const TTerm&) noexcept;

                friend std::ostream& operator<<(std::ostream&, const TTerm&) noexcept;
            private:
                uint64_t sum_ = 0;
        };
    }
}
