#pragma once
#include <cstdint>
#include <iostream>

namespace FF4 {
    namespace NUtils {

        constexpr bool IsPrime(int32_t n) noexcept {
            if (n <= 1 || n % 2 == 0) {
                return false;
            }
            for (int32_t i = 3; i * i <= n; i += 2) {
                if (n % i == 0) {
                    return false;
                }
            }
            return true;
        }

        template <int32_t Mod>
        class PrimeField {
            static int32_t binpow(int32_t value, int32_t pow) {
                int32_t res = 1;
                while(pow != 0) {
                    if (pow & 1) {
                        res = (res * 1ll * value) % Mod;
                    }
                    value = (value * 1ll * value) % Mod;
                    pow >>= 1;
                }
                return res;
            }

        public:
            PrimeField() {
                static_assert(IsPrime(Mod));
            }

            PrimeField(int32_t number)
            {
                number_ = number % Mod;
                if (number_ < 0) {
                    number_ += Mod;
                }
                static_assert(IsPrime(Mod));
            }

            friend bool operator==(const PrimeField& left, const PrimeField& right) noexcept {
                return left.number_ == right.number_;
            }

            friend bool operator!=(const PrimeField& left, const PrimeField& right) noexcept {
                return !(left == right);
            }

            bool IsPositive() const noexcept {
                return number_ != 0 && number_ != Mod - 1;
            }

            PrimeField operator+() const noexcept {
                return *this;
            }

            PrimeField operator-() const noexcept {
                return PrimeField(Mod - number_);
            }

            PrimeField& operator+=(const PrimeField& other) noexcept {
                number_ += other.number_;
                if (number_ >= Mod) {
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
                number_ = (number_ * 1ll * other.number_) % Mod;
                return *this;
            }

            friend PrimeField operator*(PrimeField left, const PrimeField& right) noexcept {
                left *= right;
                return left;
            }

            PrimeField& operator/=(const PrimeField& other) {
                assert(other.number_ != 0);
                *this *= binpow(other.number_, Mod - 2);
                return *this;
            }

            friend PrimeField operator/(PrimeField left, const PrimeField& right) {
                left /= right;
                return left;
            }

            friend std::ostream& operator<<(std::ostream& out, const PrimeField& primeField) noexcept {
                return out << primeField.number_;
            }

        private:
            int32_t number_ = 0;
        };
    }
}
