#include "rational.h"
#include <algorithm>
#include <numeric>
#include <cassert>

namespace FF4 {
    namespace NUtils {
        Rational::Rational(Integer numerator)
        : numerator_(numerator)
        {}

        Rational::Rational(Integer numerator, Integer denominator)
        : numerator_(numerator)
        , denominator_(denominator)
        {
            Normalize();
        }

        int64_t Rational::GetNumerator() const noexcept {
            return numerator_;
        }

        int64_t Rational::GetDenominator() const noexcept {
            return denominator_;
        }

        Rational Rational::operator+() const noexcept {
            return *this;
        }

        Rational Rational::operator-() const noexcept {
            return Rational(-numerator_, denominator_);
        }

        Rational& Rational::operator+=(const Rational& other) noexcept {
            numerator_ = numerator_ * other.GetDenominator() + other.GetNumerator() * denominator_;
            denominator_ *= other.GetDenominator();
            Normalize();
            return *this;
        }

        Rational operator+(Rational left, const Rational& right) noexcept {
            left += right;
            return left;
        }

        Rational& Rational::operator-=(const Rational& other) noexcept {
            numerator_ = numerator_ * other.GetDenominator() - other.GetNumerator() * denominator_;
            denominator_ *= other.GetDenominator();
            Normalize();
            return *this;
        }

        Rational operator-(Rational left, const Rational& right) noexcept {
            left -= right;
            return left;
        }

        Rational& Rational::operator*=(const Rational& other) noexcept {
            numerator_ *= other.GetNumerator();
            denominator_ *= other.GetDenominator();
            Normalize();
            return *this;
        }

        Rational operator*(Rational left, const Rational& right) noexcept {
            left *= right;
            return left;
        }

        Rational& Rational::operator/=(const Rational& other) {
            assert(other != 0);
            numerator_ *= other.GetDenominator();
            denominator_ *= other.GetNumerator();
            Normalize();
            return *this;
        }

        Rational operator/(Rational left, const Rational& right) {
            left /= right;
            return left;
        }

        bool operator<(const Rational& left, const Rational& right) noexcept {
            return left.numerator_ * right.denominator_ < right.numerator_ * left.denominator_;
        }

        bool operator>(const Rational& left, const Rational& right) noexcept {
            return right < left;
        }

        bool operator<=(const Rational& left, const Rational& right) noexcept {
            return !(right > left);
        }

        bool operator>=(const Rational& left, const Rational& right) noexcept {
            return right <= left;
        }

        bool operator==(const Rational& left, const Rational& right) noexcept {
            return left.numerator_ == right.numerator_ && left.denominator_ == right.denominator_;
        }

        bool operator!=(const Rational& left, const Rational& right) noexcept {
            return !(left == right);
        }

        bool Rational::IsPositive() const noexcept {
            return numerator_ > 0;
        }

        std::ostream& operator<<(std::ostream& out, const Rational& rational) noexcept {
            if (rational.denominator_ == 1) {
                if (rational.numerator_ < 0) {
                    return  out << " - " << -rational.numerator_;
                } else {
                    return out << rational.numerator_;
                }
            }
            if (rational.numerator_ < 0) {
                out << " -(" << -rational.numerator_;
            } else {
                out << '(' << rational.numerator_;
            }
            return out << " / " << rational.denominator_ << ')';
        }

        void Rational::Normalize() noexcept {
            if (numerator_ == 0) {
                denominator_ = 1;
                return;
            }
            if (denominator_ < 0) {
                numerator_ = -numerator_;
                denominator_ = -denominator_;
            }
            Integer gcd = std::gcd(numerator_, denominator_);
            numerator_ /= gcd;
            denominator_ /= gcd;
        }
    }
}
