#pragma once
#include <vector>
#include <cstdint>
#include <iostream>

namespace NUtils {
    class TTerm : public std::vector<uint64_t> {
        public:
            using TBase = std::vector<uint64_t>;
            using size_type = typename TBase::size_type;

            inline TTerm() : TBase() {}

            inline explicit TTerm(size_type sz) : TBase(sz) {}

            inline TTerm(std::initializer_list<uint64_t> il) : TBase(il) {
                Normalize();
            }

            bool IsDivisibleBy(const TTerm&) const noexcept;
            TTerm& operator*=(const TTerm&) noexcept;
            friend TTerm operator*(TTerm, const TTerm&) noexcept;
            TTerm& operator/=(const TTerm&) noexcept;
            friend TTerm operator/(TTerm, const TTerm&) noexcept;

            friend TTerm gcd(const TTerm&, const TTerm&) noexcept;
            friend TTerm lcm(const TTerm&, const TTerm&) noexcept;

            friend std::ostream& operator<<(std::ostream&, const TTerm&) noexcept;
        private:
            void Normalize();
    };
}
