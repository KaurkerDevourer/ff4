#pragma once
#include <vector>
#include <cstdint>
#include <iostream>

namespace NUtils {
    class TTerm : public std::vector<uint64_t> {
        private:
            void Normalize();
        public:
            using TBase = std::vector<uint64_t>;
            using size_type = typename TBase::size_type;

            inline TTerm() : TBase() {}

            inline TTerm(size_type sz) : TBase(sz) {}

            inline TTerm(std::initializer_list<uint64_t> il) : TBase(il) {
                Normalize();
                for (size_t i = 0; i < (*this).size(); i++) {
                    sum_ += (*this)[i];
                }
            }

            bool IsDivisibleBy(const TTerm&) const noexcept;
            uint64_t GetDegree() const noexcept;
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
