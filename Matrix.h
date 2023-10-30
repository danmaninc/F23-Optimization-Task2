#ifndef OPTIMIZATION_MATRIX_H
#define OPTIMIZATION_MATRIX_H

#include <vector>
#include <string>
#include <istream>
#include <iomanip>
#include <cmath>

const double MIN_COUT = 0.05;

const double ZERO = 0;

const int NUMBERS_AFTER_DOT = 2;

const double EPS = 0.001;

struct Matrix {
    int n;
    int m;
    std::vector<std::vector<double>> table;

    Matrix(const int first, const int second) {
        n = first;
        m = second;

        table.resize(n);

        for (int i = 0; i < n; ++i)
            table[i].resize(m);
    }

    Matrix(const Matrix& other) = default;
    Matrix(Matrix&& other) = default;
    ~Matrix() = default;

    friend std::istream& operator>>(std::istream& in, Matrix& matrix) {
        for (int i = 0; i < matrix.n; ++i)
            for (int j = 0; j < matrix.m; ++j)
                in >> matrix.table[i][j];

        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
        for (auto& i : matrix.table) {
            for (int j = 0; j < i.size(); ++j) {
                if (fabs(i[j]) < MIN_COUT) {
                    out << std::setprecision(NUMBERS_AFTER_DOT) << std::fixed << ZERO;
                } else {
                    out << std::setprecision(NUMBERS_AFTER_DOT) << std::fixed << i[j];
                }
                if (j != i.size() - 1) {
                    out << " ";
                }
            }
            out << std::endl;
        }
        return out;
    }

    void resize(const int new_n, const int new_m) {
        if (n != new_n)
            table.resize(new_n);

        if (m != new_m)
            for (int i = 0; i < new_n; i++)
                table[i].resize(new_m);

        n = new_n;
        m = new_m;
    }

    bool check_sum_sub(const Matrix& A, const Matrix& B) {
        return A.n == B.n && A.m == B.m;
    }

    bool check_mul(const Matrix& A, const Matrix& B) {
        return A.m == B.n;
    }

    bool check_for_determintant(const Matrix& A) {
        return A.n == A.m;
    }

    double find_determinant(const Matrix& A) {
        if (A.n == 1)
            return A.table[0][0];

        if (A.n == 2)
            return A.table[0][0] * A.table[1][1] - A.table[0][1] * A.table[1][0];

        bool sign = true;
        double result = 0;

        for (int column = 0; column < A.n; ++column) {
            Matrix matrix(A.n - 1, A.n - 1);

            for (int i = 1; i < A.n; ++i) {
                int secondPosition = 0;

                for (int j = 0; j < A.n; ++j) {
                    if (j != column) {
                        matrix.table[i - 1][secondPosition] = A.table[i][j];
                        secondPosition++;
                    }
                }
            }

            if (sign) {
                result += A.table[0][column] * find_determinant(matrix);
                sign = false;
            } else {
                result -= A.table[0][column] * find_determinant(matrix);
                sign = true;
            }

        }

        return result;
    }

    Matrix& operator= (const Matrix& other) {
        resize(other.n, other.m);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                this->table[i][j] = other.table[i][j];

        return *this;
    }

    friend Matrix operator+ (const Matrix& first, const Matrix& second) {
        Matrix res(first.n, first.m);

        for (int i = 0; i < first.n; ++i)
            for (int j = 0; j < first.m; ++j)
                res.table[i][j] = first.table[i][j] + second.table[i][j];

        return res;
    }

    friend Matrix operator- (const Matrix& first, const Matrix& second) {
        Matrix res(first.n, first.m);

        for (int i = 0; i < first.n; ++i)
            for (int j = 0; j < first.m; ++j)
                res.table[i][j] = first.table[i][j] - second.table[i][j];

        return res;
    }

    friend Matrix operator* (const double number, const Matrix& matrix) {
        Matrix res(matrix.n, matrix.m);

        for (int i = 0; i < res.n; ++i)
            for (int j = 0; j < res.m; ++j)
                res.table[i][j] = number * matrix.table[i][j];

        return res;
    }

    friend Matrix operator* (const Matrix& first, const Matrix& second) {
        Matrix res(first.n, second.m);

        for (int i = 0; i < res.n; ++i) {
            for (int j = 0; j < res.m; ++j) {
                res.table[i][j] = 0;

                for (int k = 0; k < first.m; ++k)
                    res.table[i][j] += first.table[i][k] * second.table[k][j];
            }
        }

        return res;
    }

    friend bool operator== (const Matrix& first, const Matrix& second) {
        if (first.n != second.n || first.m != second.m)
            return false;

        for (int i = 0; i < first.n; ++i)
            for (int j = 0; j < first.m; ++j)
                if (std::fabs(first.table[i][j] - second.table[i][j]) > EPS)
                    return false;

        return true;
    }

    friend bool operator!= (const Matrix& first, const Matrix& second) {
        return !(first == second);
    }

    Matrix find_transpose_matrix() {
        Matrix transposed(m, n);

        for (int i = 0; i < transposed.n; ++i) {
            for (int j = 0; j < transposed.m; ++j) {
                transposed.table[i][j] = this->table[j][i];
            }
        }
        return transposed;
    }
};

struct SquareMatrix : public Matrix {
    explicit SquareMatrix(const int n) : Matrix(n, n) {};
    SquareMatrix(const SquareMatrix& other) = default;
    SquareMatrix(SquareMatrix&& other) = default;
    ~SquareMatrix() = default;

    friend SquareMatrix operator+ (const SquareMatrix& a, const SquareMatrix& b) {
        auto* ptrA = (const Matrix*) &a;
        auto* ptrB = (const Matrix*) &b;
        auto res = *ptrA + *ptrB;
        return * (SquareMatrix *) &res;
    }

    friend SquareMatrix operator- (const SquareMatrix& a, const SquareMatrix& b) {
        auto* ptrA = (const Matrix*) &a;
        auto* ptrB = (const Matrix*) &b;
        auto res = *ptrA - *ptrB;
        return *(SquareMatrix *) &res;
    }

    friend SquareMatrix operator* (SquareMatrix& a, SquareMatrix& b) {
        auto* ptrA = (const Matrix*) &a;
        auto* ptrB = (const Matrix*) &b;
        auto res = *ptrA * *ptrB;
        return *(SquareMatrix*) &res;
    }

    SquareMatrix find_transpose_matrix() {
        auto* ptr = (Matrix*) this;
        auto res = ptr->find_transpose_matrix();
        return *(SquareMatrix*) &res;
    }

    friend std::istream& operator>>(std::istream &in, const SquareMatrix &matrix) {
        auto* ptr = (Matrix*) &matrix;
        in >> *ptr;
        return in;
    }
    friend std::ostream& operator<<(std::ostream &out, const SquareMatrix &matrix) {
        auto* ptr = (Matrix*) &matrix;
        out << *ptr;
        return out;
    }

    void resize(const int new_size) {
        if (n == new_size)
            return;

        table.resize(new_size);

        for (int i = 0; i < new_size; i++)
            table[i].resize(new_size);

        n = new_size;
        m = new_size;
    }
};

struct IdentityMatrix : public SquareMatrix {
    explicit IdentityMatrix(const int n) : SquareMatrix(n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    this->table[i][j] = 1;
                } else {
                    this->table[i][j] = 0;
                }
            }
        }
    }

    IdentityMatrix(const IdentityMatrix& other) = default;
    IdentityMatrix(IdentityMatrix&& other) = default;
    ~IdentityMatrix() = default;

    friend std::ostream& operator<<(std::ostream& out, const IdentityMatrix& matrix) {
        auto* ptr = (Matrix *) &matrix;
        out << *ptr;
        return out;
    }

    friend IdentityMatrix operator- (const IdentityMatrix& a, const IdentityMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA - *ptrB;
        return *(IdentityMatrix*) &res;
    }

    friend IdentityMatrix operator+ (const IdentityMatrix& a, const IdentityMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA + *ptrB;
        return *(IdentityMatrix*) &res;
    }

    friend IdentityMatrix operator* (const IdentityMatrix& a, const IdentityMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA * *ptrB;
        return *(IdentityMatrix*) &res;
    }
};

struct ColumnVector : public Matrix {
    explicit ColumnVector(const int first) : Matrix(first, 1) {}
    ColumnVector(const ColumnVector& vector) = default;
    ColumnVector(ColumnVector&& vector) = default;
    ~ColumnVector() = default;

    friend ColumnVector operator+ (const ColumnVector& a, const ColumnVector& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA + *ptrB;
        return * (ColumnVector *) &res;
    }

    friend ColumnVector operator- (ColumnVector& a, ColumnVector& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA - *ptrB;
        return * (ColumnVector *) &res;
    }

    friend ColumnVector operator* (ColumnVector& a, ColumnVector& b) {
        auto* ptrA = (Matrix*) &a;
        auto* ptrB = (Matrix*) &b;
        auto res = *ptrA * *ptrB;
        return *(ColumnVector*) &res;
    }

    friend Matrix operator* (const IdentityMatrix& a, const ColumnVector& b) {
        auto* ptrA = (Matrix*) &a;
        auto* ptrB = (Matrix*) &b;
        auto res = *ptrA * *ptrB;
        return *(Matrix*) &res;
    }

    friend Matrix operator* (const Matrix& a, const ColumnVector& b) {
        auto* ptrA = &a;
        auto* ptrB = (Matrix*) &b;
        auto res = *ptrA * *ptrB;
        return *(Matrix*) &res;
    }

    friend std::istream& operator>>(std::istream &in, ColumnVector& matrix) {
        auto* ptr = (Matrix *) &matrix;
        in >> *ptr;
        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const ColumnVector& matrix) {
        auto* ptr = (Matrix *) &matrix;
        out << *ptr;
        return out;
    }

    void resize(const int newSize) {
        n = newSize;
        table.resize(newSize);
    }
};


#endif //OPTIMIZATION_MATRIX_H
