#pragma once
#include <cstdint>
#include <iostream>

namespace FF4 {
    namespace NUtils {
        using numerator = int64_t;
        using denominator = int64_t;
        class Rational {
        public:
            Rational(numerator numerator = 0, denominator denominator = 1);

            numerator GetNumerator() const noexcept;
            denominator GetDenominator() const noexcept;
            bool MoreThanZero() const noexcept ;

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

            numerator numerator_;
            denominator denominator_;
        };
    }
}
