#pragma once
#include <cstdint>
#include <iostream>

namespace NUtils {
    uint64_t binpow(uint64_t value, int32_t pow, int32_t mod) {
        if (pow == 0) {
            return 1;
        }
        if (pow % 2 == 0) {
            uint64_t half = binpow(value, pow/2, mod);
            return (half * half) % mod;
        } else {
            return (value * binpow(value, pow - 1, mod)) % mod;
        }
    }

    constexpr bool IsPrime(uint32_t n) noexcept {
        if (n == 1) {
            return false;
        }
        for (uint32_t i = 2; i * i <= n; i++) {
            if (n % i == 0) {
                return false;
            }
        }
        return true;
    }

    template <uint32_t Mod>
    class PrimeField {
    public:
        PrimeField(uint64_t number)
        : number_(number % Mod)
        {
            static_assert(IsPrime(Mod));
        }

        friend bool operator==(const PrimeField& left, const PrimeField& right) noexcept {
            return left.number_ == right.number_;
        }

        friend bool operator!=(const PrimeField& left, const PrimeField& right) noexcept {
            return !(left == right);
        }

        PrimeField operator+() const noexcept {
            return *this;
        }

        PrimeField operator-() const noexcept {
            return PrimeField(Mod - number_);
        }

        PrimeField& operator+=(const PrimeField& other) noexcept {
            number_ += other.number_;
            if (number_ > Mod) {
                number_ -= Mod;
            }
            return *this;
        }

        friend PrimeField operator+(PrimeField left, const PrimeField& right) noexcept {
            left += right;
            return left;
        }

        PrimeField& operator-=(const PrimeField& other) noexcept {
            number_ -= other.number_;
            if (number_ < 0) {
                number_ += Mod;
            }
            return *this;
        }

        friend PrimeField operator-(PrimeField left, const PrimeField& right) noexcept {
            left -= right;
            return left;
        }

        PrimeField& operator*=(const PrimeField& other) noexcept {
            number_ *= other.number_;
            number_ %= Mod;
            return *this;
        }

        friend PrimeField operator*(PrimeField left, const PrimeField& right) noexcept {
            left *= right;
            return left;
        }

        PrimeField& operator/=(const PrimeField& other) {
            *this *= binpow(other.number_, Mod - 2, Mod);
            return *this;
        }

        friend PrimeField operator/(PrimeField left, const PrimeField& right) {
            left /= right;
            return left;
        }

        friend std::ostream& operator<<(std::ostream& out, const PrimeField& primeField) noexcept {
            return out << "(" << primeField.number_ << ")%" << primeField.mod;
        }

    private:
        uint64_t number_;
        static const uint32_t mod = Mod;
    };
}
