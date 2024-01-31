#include "rational.h"
#include "math.h"

Rational::Rational(int64_t numerator, int64_t denominator) {
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    // MakeIrreducible();
}

Rational::GetNumerator() const noexcept {
    return numerator_;
}

Rational::GetDenominator() const noexcept {
    return denominator_;
}

bool operator<(const Rational& left, const Rational& right) noexcept {
    return left.numerator_ * right.denominator_ < right.numerator_ * left.denominator_;
}

bool operator>(const Rational& left, const Rational& right) noexcept {
    return left.numerator_ * right.denominator_ > right.numerator_ * left.denominator_;
}

bool operator<=(const Rational& left, const Rational& right) noexcept {
    return left.numerator_ * right.denominator_ <= right.numerator_ * left.denominator_;
}

bool operator>=(const Rational& left, const Rational& right) noexcept {
    return left.numerator_ * right.denominator_ >= right.numerator_ * left.denominator_;
}

bool operator==(const Rational& left, const Rational& right) noexcept {
    return left.numerator_ == right.numerator_ && left.denominator_ == right.denominator_;
}

bool operator!=(const Rational& left, const Rational& right) noexcept {
    return left.numerator_ != right.numerator_ || left.denominator_ != right.denominator_;
}

/* Rational& Rational::operator+() const noexcept {
    return *this;
}

Rational Rational::operator-() const noexcept {
    return Rational(-GetNumerator(), GetDenominator());
} */

Rational& Rational::operator+=(const Rational& other) noexcept {
    numerator_ = GetNumerator() * other.GetDenominator() + other.GetNumerator() * GetDenominator();
    denominator_ *= other.GetDenominator();
    MakeIrreducible();
    return *this;
}

Rational operator+(Rational left, const Rational& right) noexcept {
    left += right;
    return left;
}

Rational& Rational::operator-=(const Rational& other) noexcept {
    numerator_ = GetNumerator() * other.GetDenominator() - other.GetNumerator() * GetDenominator();
    denominator_ *= other.GetDenominator();
    MakeIrreducible();
    return *this;
}

Rational operator-(Rational left, const Rational& right) noexcept {
    left -= right;
    return left;
}

Rational& Rational::operator*=(const Rational& other) noexcept {
    numerator_ *= other.GetNumerator();
    denominator_ *= other.GetDenominator();
    MakeIrreducible();
    return *this;
}

Rational operator*(Rational left, const Rational& right) noexcept {
    left *= right;
    return left;
}

Rational& Rational::operator/=(const Rational& other) {
    numerator_ *= other.GetDenominator();
    denominator_ *= other.GetNumerator();
    MakeIrreducible();
    return *this;
}

Rational operator/(Rational left, const Rational& right) {
    left /= right;
    return left;
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
        int64_t tmp = pow(number, power / 2);
        return tmp * tmp;
    }
}

std::ostream& operator<<(std::ostream& out, const Rational& rational) noexcept {
    if (rational.numerator_ < 0) {
        out << "-(" << -rational.numerator_;
    } else {
        out << '(' << rational.numerator_;
    }
    return out << " / " << rational.denominator_ << ')';
}


void Rational::MakeIrreducible() noexcept {
    int64_t gcd = NMath::gcd(numerator_, denominator_);
    numerator_ /= gcd;
    denominator_ /= gcd;
}