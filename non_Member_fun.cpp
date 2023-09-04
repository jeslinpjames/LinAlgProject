#include <vector>
#include <initializer_list>
#include <random> 
#include <fstream>   
#include <sstream> 
#include <algorithm>
#include<iostream>
#include<iomanip>
#include "matrix.cpp"
using namespace std;

static Matrix  diag_matrix(const vector<double>& diagValues) {
    size_t size = diagValues.size();
    Matrix diagonalMatrix(size, size);
    for (size_t i = 0; i < size; ++i) {
        diagonalMatrix.mat_data[i * size + i] = diagValues[i];
    }
    
    return diagonalMatrix;
}
bool is_triangular() const {
    if (r != c) return false;
    for (size_t i = 0; i < r; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (mat_data[i * c + j] != 0.0) return false;
        }
    }
    for (size_t i = 0; i < r; ++i) {
        for (size_t j = i + 1; j < c; ++j) {
            if (mat_data[i * c + j] != 0.0) return false;
        }
    }
    return true;
}
Matrix LU_Decomposition(const Matrix& matrix) {
    if (matrix.r != matrix.c) {
        throw runtime_error("LU decomposition requires a square matrix");
    }

    size_t n = matrix.r;
    Matrix L(n, n);
    Matrix U(n, n);

    for (size_t i = 0; i < n; ++i) {
        L.mat_data[i * n + i] = 1.0;
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L.mat_data[i * n + k] * U.mat_data[k * n + j];
            }
            U.mat_data[i * n + j] = matrix.mat_data[i * n + j] - sum;
        }

        for (size_t j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L.mat_data[j * n + k] * U.mat_data[k * n + i];
            }
            L.mat_data[j * n + i] = (matrix.mat_data[j * n + i] - sum) / U.mat_data[i * n + i];
        }
    }

    return L;
}
std::pair<Matrix, Matrix> QR_factorization() const {
        Matrix Q(r, c);
        Matrix R(r, c);
        std::vector<std::vector<double>> u(r, std::vector<double>(c));
        std::vector<std::vector<double>> v(r, std::vector<double>(c));
        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; j < c; ++j) {
                if (i == j) {
                    R.mat_data[i * c + j] = 1.0;
                } else {
                    R.mat_data[i * c + j] = 0.0;
                }
                Q.mat_data[i * c + j] = mat_data[i * c + j];
            }
        }
        for (size_t j = 0; j < c; ++j) {
            for (size_t k = 0; k < j; ++k) {
                double dot_product = 0.0;
                for (size_t i = 0; i < r; ++i) {
                    dot_product += Q.mat_data[i * c + j] * Q.mat_data[i * c + k];
                }

                for (size_t i = 0; i < r; ++i) {
                    u[i][j] = Q.mat_data[i * c + j] - dot_product * Q.mat_data[i * c + k];
                }
            }
            double norm = 0.0;
            for (size_t i = 0; i < r; ++i) {
                norm += u[i][j] * u[i][j];
            }
            norm = std::sqrt(norm);
            for (size_t i = 0; i < r; ++i) {
                Q.mat_data[i * c + j] = u[i][j] / norm;
            }
            for (size_t i = 0; i < r; ++i) {
                v[i][j] = Q.mat_data[i * c + j];
            }
            for (size_t k = j; k < c; ++k) {
                double dot_product = 0.0;
                for (size_t i = 0; i < r; ++i) {
                    dot_product += v[i][j] * mat_data[i * c + k];
                }
                for (size_t i = 0; i < r; ++i) {
                    R.mat_data[j * c + k] = dot_product;
                }
            }
        }
        return std::make_pair(Q, R);
    }
std::pair<Matrix, Matrix> eigen_decomposition() const {
    Matrix A = *this;
    Matrix Q(r, c);
    Matrix R(r, c);

    for (size_t i = 0; i < 50; ++i) {
        std::pair<Matrix, Matrix> qr = A.qr_factorization();
        Q = qr.first;
        R = qr.second;
        A = R * Q;
    }

    return std::make_pair(A, Q);
}
std::tuple<Matrix, Matrix, Matrix> singular_value_decomposition() const {
    Matrix A = *this;
    size_t m = A.shape().first;
    size_t n = A.shape().second;
    Matrix U(m, m, true);
    Matrix S(m, n);
    Matrix V(n, n, true);
    for (size_t i = 0; i < 50; ++i) {
        std::pair<Matrix, Matrix> qr = A.qr_factorization();
        A = qr.second * qr.first;
        Matrix Qt = qr.first.T();
        U *= Qt;
        if (i > 0) {
            S *= qr.second;
        } else {
            S = qr.second;
        }
    }
    return std::make_tuple(U, S, A.T());
}