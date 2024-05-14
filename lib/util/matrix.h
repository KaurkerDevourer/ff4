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
                A_.resize(pivots_ * pivots_);
                B_.resize((M_ - pivots_) * pivots_);
                C_.resize((N_ - pivots_) * pivots_);
                D_.resize((M_ - pivots_) * (N_ - pivots_));
            }

            TCoef& operator()(size_t i, size_t j) {
                bool top = (i < pivots_);
                bool left = (j < pivots_);
                if (left && top) {
                    return A_[i * pivots_ + j];
                } else if (!left && top) {
                    return B_[i * (M_ - pivots_) + (j - pivots_)];
                } else if (left && !top) {
                    return C_[(i - pivots_) * pivots_ + j];
                } else {
                    return D_[(i - pivots_) * (N_ - pivots_) + (j - pivots_)];
                }
            }

            const TCoef& operator()(size_t i, size_t j) const {
                bool top = (i < pivots_);
                bool left = (j < pivots_);
                if (left && top) {
                    return A_[i * pivots_ + j];
                } else if (!left && top) {
                    return B_[i * (M_ - pivots_) + (j - pivots_)];
                } else if (left && !top) {
                    return C_[(i - pivots_) * pivots_ + j];
                } else {
                    return D_[(i - pivots_) * (N_ - pivots_) + (j - pivots_)];
                }
            }

            size_t N_;
            size_t M_;
            size_t pivots_;
        private:
            std::vector<TCoef> A_;
            std::vector<TCoef> B_;
            std::vector<TCoef> C_;
            std::vector<TCoef> D_;
        };
    }
}
