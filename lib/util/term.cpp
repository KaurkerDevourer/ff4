#include "term.h"
#include <cassert>

namespace FF4 {
    namespace NUtils {
        Term::Term(std::initializer_list<uint64_t> il) {
            data_.resize(il.size());
            size_t i = 0;
            for (uint64_t x : il) {
                data_[i] = x;
                sum_ += x;
                i++;
            }
            Normalize();
        }

        uint64_t& Term::operator[](size_t i) {
            return data_[i];
        }

        const uint64_t& Term::operator[](size_t i) const {
            return data_[i];
        }

        const std::vector<uint64_t>& Term::GetData() const {
            return data_;
        }

        void Term::resize(size_t sz) {
            data_.resize(sz);
        }

        size_t Term::size() const noexcept {
            return data_.size();
        }

        void Term::reserve(size_t sz) {
            data_.reserve(sz);
        }

        void Term::push_back(uint64_t value) {
            data_.push_back(value);
        }

        bool Term::IsOne() const noexcept {
            return sum_ == 0;
        }

        bool Term::IsDivisibleBy(const Term& other) const noexcept {
            if (other.size() > size()) {
                return false;
            }
            for (size_t i = 0; i < other.size(); i++) {
                if (other[i] > data_[i]) {
                    return false;
                }
            }
            return true;
        };

        uint64_t Term::TotalDegree() const noexcept {
            return sum_;
        }

        std::vector<Term> Term::GetAllDivisors() noexcept {
            std::vector<Term> divisors = {(*this)};
            FillDivisors(divisors, 0);
            for (auto& div : divisors) {
                div.Normalize();
            }
            return divisors;
        }

        Term& Term::operator/=(const Term& other) noexcept {
            assert(other.size() <= size());
            for (size_t i = 0; i < other.size(); i++) {
                assert(data_[i] >= other[i]);
                data_[i] -= other[i];
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
                data_[i] += other[i];
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
                if (term[i] == 0) {
                    continue;
                }
                out << "x_" << i;
                if (term[i] != 1) {
                    out << "^{" << term[i] << "}";
                }
            }
            return out;
        }

        bool operator<(const Term& a, const Term& b) noexcept {
            return a.data_ < b.data_;
        }

        bool operator>(const Term& a, const Term& b) noexcept {
            return a.data_ > b.data_;
        }

        bool operator<=(const Term& a, const Term& b) noexcept {
            return a.data_ <= b.data_;
        }

        bool operator>=(const Term& a, const Term& b) noexcept {
            return a.data_ >= b.data_;
        }

        bool operator==(const Term& a, const Term& b) noexcept {
            return a.sum_ == b.sum_ && a.data_ == b.data_;
        }

        bool operator!=(const Term& a, const Term& b) noexcept {
            return a.sum_ != b.sum_ || a.data_ != b.data_;
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

        void Term::Normalize() noexcept {
            while(data_.size() > 1 && data_.back() == 0) {
                data_.pop_back();
            }
        }

        void Term::FillDivisors(std::vector<Term>& divisors, size_t idx) noexcept {
            for (size_t i = idx; i < data_.size(); i++) {
                if (data_[i] != 0) {
                    data_[i]--;
                    sum_--;
                    divisors.push_back((*this));
                    FillDivisors(divisors, i);
                    data_[i]++;
                    sum_++;
                }
            }
        }
    }
}
