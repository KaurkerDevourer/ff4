#pragma once
#include <cstdint>
#include <iostream>

namespace FF4 {
    namespace NUtils {
        class Rational {
            using Integer = int64_t;
        public:
            Rational(Integer numerator = 0, Integer denominator = 1);

            Integer GetNumerator() const noexcept;
            Integer GetDenominator() const noexcept;
            bool MoreThanZero() const noexcept;

            Rational operator+() const noexcept;
            Rational operator-() const noexcept;

            friend bool operator==(const Rational&, const Rational&) noexcept;
            friend bool operator!=(const Rational&, const Rational&) noexcept;

            Rational& operator+=(const Rational&) noexcept;
            friend Rational operator+(Rational, const Rational&) noexcept;

            Rational& operator-=(const Rational&) noexcept;
            friend Rational operator-(Rational, const Rational&) noexcept;

            Rational& operator*=(const Rational&) noexcept;
            friend Rational operator*(Rational, const Rational&) noexcept;

            Rational& operator/=(const Rational&);
            friend Rational operator/(Rational, const Rational&);

            friend bool operator<(const Rational&, const Rational&) noexcept;
            friend bool operator>(const Rational&, const Rational&) noexcept;
            friend bool operator<=(const Rational&, const Rational&) noexcept;
            friend bool operator>=(const Rational&, const Rational&) noexcept;

            friend Rational pow(const Rational&, int64_t) noexcept;

            friend std::ostream& operator<<(std::ostream&, const Rational&) noexcept;

        private:
            void Normalize() noexcept;

            Integer numerator_;
            Integer denominator_;
        };
    }
}
