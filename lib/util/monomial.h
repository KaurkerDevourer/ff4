#pragma once
#include "rational.h"
#include <vector>

using TTerm = std::vector<int64_t>;
class Monomial {
    explicit Monomial(TTerm&& term, Rational coef = 1);

    TTerm& GetTerm() const noexcept;
    Rational& GetCoef() const noexcept;

    
private:
    TTerm term_;
    Rational coef_;
}