#include "rational.h"
#include <algorithm>
#include <numeric>

namespace NUtils {
    Rational::Rational(numerator numerator, denominator denominator)
    : numerator_(numerator)
    , denominator_(denominator)
    {
        Normalize();
    }

    numerator Rational::GetNumerator() const noexcept {
        return numerator_;
    }

    denominator Rational::GetDenominator() const noexcept {
        return denominator_;
    }

    bool Rational::MoreThanZero() const noexcept {
        return numerator_ > 0;
    }

    Rational Rational::operator+() const noexcept {
        return *this;
    }

    Rational Rational::operator-() const noexcept {
        return Rational(-numerator_, denominator_);
    }

    bool operator==(const Rational& left, const Rational& right) noexcept {
        return left.numerator_ == right.numerator_ && (left.numerator_ == 0 || left.denominator_ == right.denominator_);
    }

    bool operator!=(const Rational& left, const Rational& right) noexcept {
        return !(left == right);
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
        *this += -other;
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
        //assert(other != 0);
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
        return left.numerator_ * right.denominator_ <= right.numerator_ * left.denominator_;
    }

    bool operator>=(const Rational& left, const Rational& right) noexcept {
        return right <= left;
    }

    Rational pow(const Rational& number, int64_t power) noexcept {
        if (power < 0) {
            return pow(1 / number, -power);
        } else if (power == 0) {
            return 1;
        }
        if (power % 2 != 0) {
            return number * pow(number, power - 1);
        } else {
            Rational tmp = pow(number, power / 2);
            return tmp * tmp;
        }
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
        if (denominator_ < 0) {
            numerator_ = -numerator_;
            denominator_ = -denominator_;
        }
        int64_t gcd = std::gcd(abs(numerator_), abs(denominator_));
        numerator_ /= gcd;
        denominator_ /= gcd;
    }
}
