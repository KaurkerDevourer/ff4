#include "term.h"

namespace NUtils {
    TTerm gcd(const TTerm& left, const TTerm& right) noexcept {
        size_t sz = std::min(left.size(), right.size());
        TTerm term(sz);
        for (size_t i = 0; i < sz; i++) {
            term[i] = std::min(left[i], right[i]);
        }
        return std::move(term);
    }

    TTerm lcm(const TTerm& left, const TTerm& right) noexcept {
        size_t sz = std::min(left.size(), right.size());
        size_t lcm_sz = std::max(left.size(), right.size());
        TTerm term(lcm_sz);
        for (size_t i = 0; i < sz; i++) {
            term[i] = std::max(left[i], right[i]);
        }
        for (size_t i = left.size(); i < lcm_sz; i++) {
            term[i] = right[i];
        }
        for (size_t i = right.size(); i < lcm_sz; i++) {
            term[i] = left[i];
        }
        return std::move(term);
    }

    bool TTerm::IsDivisibleBy(const TTerm& other) const noexcept {
        if (other.size() > size()) {
            return false;
        }
        for (size_t i = 0; i < other.size(); i++) {
            if (other[i] > (*this)[i]) {
                return false;
            }
        }
        return true;
    };

    TTerm& TTerm::operator/=(const TTerm& other) noexcept {
        if (other.size() > size()) {
            resize(other.size());
        }
        for (size_t i = 0; i < other.size(); i++) {
            (*this)[i] -= other[i];
        }
        Normalize();
        return *this;
    }

    TTerm& TTerm::operator*=(const TTerm& other) noexcept {
        if (other.size() > size()) {
            resize(other.size());
        }
        for (size_t i = 0; i < other.size(); i++) {
            (*this)[i] += other[i];
        }
        return *this;
    }

    TTerm operator*(TTerm left, const TTerm& right) noexcept {
        left *= right;
        return left;
    }

    TTerm operator/(TTerm left, const TTerm& right) noexcept {
        left /= right;
        return left;
    }

    std::ostream& operator<<(std::ostream& out, const TTerm& term) noexcept {
        for (size_t i = 0; i < term.size(); i++) {
            if (term[i] == 1) {
                out << "x_" << i;
            } else if (term[i] != 0) {
                out << "x_" << i << "^{" << term[i] << "}";
            }
        }
        return out;
    }

    void TTerm::Normalize() {
        while(!empty() && back() == 0) {
            pop_back();
        }
    }
}
