#pragma once
#include <vector>

namespace FF4 {
    namespace NUtils {
        template <typename TCoef>
        class Matrix {
        public:
            Matrix() = delete;

            Matrix(size_t n, size_t m, size_t pivots)
            : N_(n)
            , M_(m)
            , pivots_(pivots)
            {
                data_.resize(N_ * M_);
                offset_b_ = pivots_ * pivots_;
                offset_c_ = M_ * pivots_;
                offset_d_ = offset_c_ + (N_ - pivots_) * pivots_;
            }

            TCoef& operator()(size_t i, size_t j) {
                bool top = (i < pivots_);
                bool left = (j < pivots_);
                if (top && left) { // A
                    return data_[i * pivots_ + j];
                }
                if (top && !left) { // B
                    return data_[offset_b_ + i * (M_ - pivots_) + (j - pivots_)];
                }
                if (left) {
                    return data_[offset_c_ + (i - pivots_) * pivots_ + j];
                }
                return data_[offset_d_ + (i - pivots_) * (M_ - pivots_) + (j - pivots_)];
            }

            const TCoef& operator()(size_t i, size_t j) const {
                bool top = (i < pivots_);
                bool left = (j < pivots_);
                if (top && left) { // A
                    return data_[i * pivots_ + j];
                }
                if (top && !left) { // B
                    return data_[offset_b_ + i * (M_ - pivots_) + (j - pivots_)];
                }
                if (left) {
                    return data_[offset_c_ + (i - pivots_) * pivots_ + j];
                }
                return data_[offset_d_ + (i - pivots_) * (M_ - pivots_) + (j - pivots_)];
            }

            size_t N_;
            size_t M_;
            size_t pivots_;
        private:
            std::vector<TCoef> data_;
            size_t offset_b_;
            size_t offset_c_;
            size_t offset_d_;
        };
    }
}
