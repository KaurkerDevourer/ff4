#include "term.h"

namespace FF4 {
    namespace NUtils {
        Term::Term(std::initializer_list<uint64_t> il) : TBase(il) {
            Normalize();
            for (size_t i = 0; i < (*this).size(); i++) {
                sum_ += (*this)[i];
            }
        }

        Term gcd(const Term& left, const Term& right) noexcept {
            size_t sz = std::min(left.size(), right.size());
            Term term;
            term.reserve(sz);
            for (size_t i = 0; i < sz; i++) {
                term.push_back(std::min(left[i], right[i]));
                term.sum_ += term[i];
            }
            term.Normalize();
            return term;
        }

        Term lcm(const Term& left, const Term& right) noexcept {
            size_t sz = std::min(left.size(), right.size());
            size_t lcm_sz = std::max(left.size(), right.size());
            Term term;
            term.reserve(lcm_sz);
            for (size_t i = 0; i < sz; i++) {
                term.push_back(std::max(left[i], right[i]));
                term.sum_ += term[i];
            }
            for (size_t i = left.size(); i < lcm_sz; i++) {
                term.push_back(right[i]);
                term.sum_ += term[i];
            }
            for (size_t i = right.size(); i < lcm_sz; i++) {
                term.push_back(left[i]);
                term.sum_ += term[i];
            }
            term.Normalize();
            return term;
        }

        bool Term::IsDivisibleBy(const Term& other) const noexcept {
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

        uint64_t Term::TotalDegree() const noexcept {
            return sum_;
        }

        Term& Term::operator/=(const Term& other) noexcept {
            assert(other.size() <= size());
            for (size_t i = 0; i < other.size(); i++) {
                assert((*this)[i] >= other[i]);
                (*this)[i] -= other[i];
            }
            sum_ -= other.sum_;
            Normalize();
            return *this;
        }

        Term& Term::operator*=(const Term& other) noexcept {
            if (other.size() > size()) {
                resize(other.size());
            }
            for (size_t i = 0; i < other.size(); i++) {
                (*this)[i] += other[i];
            }
            sum_ += other.sum_;
            return *this;
        }

        Term operator*(Term left, const Term& right) noexcept {
            left *= right;
            return left;
        }

        Term operator/(Term left, const Term& right) noexcept {
            left /= right;
            return left;
        }

        std::ostream& operator<<(std::ostream& out, const Term& term) noexcept {
            if (term.size() == 1 && term[0] == 0) {
                return out << "1";
            }
            for (size_t i = 0; i < term.size(); i++) {
                out << "x_" << i;
                if (term[i] != 0) {
                    out << "^{" << term[i] << "}";
                }
            }
            return out;
        }

        void Term::Normalize() {
            while(size() > 1 && back() == 0) {
                pop_back();
            }
        }
    }
}
