#pragma once
#include "rational.h"
#include <vector>

using TTerm = std::vector<int64_t>;
class Monomial {
public:
    explicit Monomial(TTerm&& term, Rational coef = 1);

    TTerm& GetTerm() const noexcept;
    Rational& GetCoef() const noexcept;

    bool IsDivisibleBy(const Monomial&) const noexcept __attribute__((pure));

    friend bool operator<(const Monomial&, const Monomial&) noexcept;
    friend bool operator>(const Monomial&, const Monomial&) noexcept;
    friend bool operator<=(const Monomial&, const Monomial&) noexcept;
    friend bool operator>=(const Monomial&, const Monomial&) noexcept;

    friend bool operator==(const Monomial&, const Monomial&) noexcept;
    friend bool operator!=(const Monomial&, const Monomial&) noexcept;

    Monomial operator+() const noexcept;
    Monomial operator-() const noexcept;

    Monomial& operator*=(const Monomial&) noexcept;
    friend Monomial operator*(Monomial, const Monomial&) noexcept;

    Monomial& operator*=(const Rational&) noexcept;
    friend Monomial operator*(Monomial, const Rational&) noexcept;

    Monomial& operator/=(const Monomial&);
    friend Monomial operator/(Monomial, const Monomial&);

    friend Monomial gcd(const Monomial&, const Monomial&) noexcept;

    friend Monomial lcm(const Monomial&, const Monomial&) noexcept;
private:
    static void Normalize(TTerm&);

private:
    TTerm term_;
    Rational coef_;
}