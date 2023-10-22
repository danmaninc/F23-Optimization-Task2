#ifndef OPTIMIZATION_MATRIX_H
#define OPTIMIZATION_MATRIX_H

#include <vector>
#include <string>
#include <istream>
#include <iomanip>

const double MIN_COUT = 0.05;

const double ZERO = 0;

const int NUMBERS_AFTER_DOT = 2;

struct Matrix {
    std::size_t n;
    std::size_t m;
    std::vector<std::vector<double>> table;

    std::vector<std::string> list_of_all_vars;
    std::vector<std::string> list_of_basic_vars;

    Matrix(const std::size_t n, const std::size_t m) {
        this->n = n;
        this->m = m;

        table.resize(n);

        for (auto& r : table)
            r.resize(m, 0);
    }

    friend std::istream& operator>>(std::istream& in, Matrix& matrix) {
        for (auto& r : matrix.table)
            for (auto& e : r)
                in >> e;

        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
        out << matrix.list_of_all_vars[0] << " " << "\t";

        for (int i = 1; i < matrix.list_of_all_vars.size(); i++)
            out << matrix.list_of_all_vars.at(i) << " " << "\t";

        out << std::endl;
        int counter = 0;

        for (auto& i : matrix.table) {
            out << matrix.list_of_basic_vars[counter++] << "\t";

            for (double j : i)
                out << std::setprecision(NUMBERS_AFTER_DOT)
                    << std::fixed
                    << (std::abs(j) < MIN_COUT ? ZERO : j)
                    << "\t";

            out << std::endl << std::endl;
        }

        return out;
    }

    Matrix& operator= (const Matrix& matrix) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                table[i][j] = matrix.table[i][j];

        return *this;
    }
};


class SquareMatrix : public Matrix {
public:
    explicit SquareMatrix(int n) : Matrix(n, n) {};
    friend SquareMatrix operator+ (SquareMatrix& a, SquareMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA + *ptrB;
        return * (SquareMatrix *) &res;
    }
    friend SquareMatrix operator- (SquareMatrix& a, SquareMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA - *ptrB;
        return * (SquareMatrix *) &res;
    }
    friend SquareMatrix operator* (SquareMatrix& a, SquareMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA * *ptrB;
        return * (SquareMatrix *) &res;
    }
    SquareMatrix findTransposeMatrix() {
        auto* ptr = (Matrix *) this;
        auto res = ptr->findTransposeMatrix();
        return * (SquareMatrix *) &res;
    }
    friend istream& operator>>(istream &in, SquareMatrix &matrix) {
        auto* ptr = (Matrix *) &matrix;
        in >> *ptr;
        return in;
    }
    friend ostream& operator<<(ostream &out, SquareMatrix &matrix) {
        auto* ptr = (Matrix *) &matrix;
        out << *ptr;
        return out;
    }
    void resize(int newSize) {
        vector<double> a;
        for (int i = 0; i < newSize; ++i) {
            a.push_back(0);
        }

        for (int i = 0; i < n; ++i) {
            for (int j = n; j < newSize; ++j) {
                table[i].push_back(0);
            }
        }
        for (int i = n; i < newSize; ++i) {
            table.push_back(a);
        }
        n = newSize;
        m = newSize;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    explicit IdentityMatrix(int n) : SquareMatrix(n) {
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
    friend ostream& operator<<(ostream &out, IdentityMatrix &matrix) {
        auto* ptr = (Matrix *) &matrix;
        out << *ptr;
        return out;
    }
    friend IdentityMatrix operator- (IdentityMatrix& a, IdentityMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA - *ptrB;
        return * (IdentityMatrix *) &res;
    }
    friend IdentityMatrix operator+ (IdentityMatrix& a, IdentityMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA + *ptrB;
        return * (IdentityMatrix *) &res;
    }
    friend IdentityMatrix operator* (IdentityMatrix& a, IdentityMatrix& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA * *ptrB;
        return * (IdentityMatrix *) &res;
    }
};

class ColumnVector : public Matrix {
public:
    explicit ColumnVector(int first) : Matrix(first, 1) {}
    friend ColumnVector operator+ (ColumnVector& a, ColumnVector& b) {
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
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA * *ptrB;
        return * (ColumnVector *) &res;
    }
    friend istream& operator>>(istream &in, ColumnVector &matrix) {
        auto* ptr = (Matrix *) &matrix;
        in >> *ptr;
        return in;
    }
    friend ostream& operator<<(ostream &out, ColumnVector &matrix) {
        auto* ptr = (Matrix *) &matrix;
        out << *ptr;
        return out;
    }
    void resize(int newSize) {
        vector<double> a;
        a.push_back(0);
        for (int i = n; i < newSize; ++i) {
            table.push_back(a);
        }
        n = newSize;
    }
};


#endif //OPTIMIZATION_MATRIX_H
