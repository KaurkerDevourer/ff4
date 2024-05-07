#include "term.h"

namespace FF4 {
    namespace NUtils {
        TTerm::TTerm(std::initializer_list<uint64_t> il) : TBase(il) {
            Normalize();
            for (size_t i = 0; i < (*this).size(); i++) {
                sum_ += (*this)[i];
            }
        }

        TTerm gcd(const TTerm& left, const TTerm& right) noexcept {
            size_t sz = std::min(left.size(), right.size());
            TTerm term;
            term.reserve(sz);
            for (size_t i = 0; i < sz; i++) {
                term.push_back(std::min(left[i], right[i]));
                term.sum_ += term[i];
            }
            term.Normalize();
            return term;
        }

        TTerm lcm(const TTerm& left, const TTerm& right) noexcept {
            size_t sz = std::min(left.size(), right.size());
            size_t lcm_sz = std::max(left.size(), right.size());
            TTerm term;
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

        uint64_t TTerm::TotalDegree() const noexcept {
            return sum_;
        }

        TTerm& TTerm::operator/=(const TTerm& other) noexcept {
            assert(other.size() <= size());
            for (size_t i = 0; i < other.size(); i++) {
                assert((*this)[i] >= other[i]);
                (*this)[i] -= other[i];
            }
            sum_ -= other.sum_;
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
            sum_ += other.sum_;
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
            if (term.size() == 1 && term[0] == 0) {
                return out << "1";
            }
            for (size_t i = 0; i < term.size(); i++) {
                out << "x_" << i;
                if (term[i] != 0) {
                    out << term[i] << "}";
                }
            }
            return out;
        }

        void TTerm::Normalize() {
            while(size() > 1 && back() == 0) {
                pop_back();
            }
        }
    }
}
