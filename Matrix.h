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

class Matrix {
public:
    int n;
    int m;
    std::vector<std::vector<double>> table;

    Matrix(int first, int second) {
        n = first;
        m = second;

        table.resize(n);
        for (int i = 0; i < n; ++i) {
            table[i].resize(m);
        }
    }


    friend std::istream& operator>>(std::istream &in, Matrix &matrix) {
        for (int i = 0; i < matrix.n; ++i) {
            for (int j = 0; j < matrix.m; ++j) {
                in >> matrix.table[i][j];
            }
        }
        return in;
    }

    friend std::ostream& operator<<(std::ostream &out, Matrix &matrix) {
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

    bool checkSumSub(Matrix& A, Matrix& B) {
        return A.n == B.n && A.m == B.m;
    }

    bool checkMul(Matrix& A, Matrix& B) {
        return A.m == B.n;
    }

    bool checkForDeterminant(Matrix& A) {
        return A.n == A.m;
    }
    double findDeterminant(Matrix& A) {
        if (A.n == 1) {
            return A.table[0][0];
        }
        if (A.n == 2) {
            return A.table[0][0] * A.table[1][1] - A.table[0][1] * A.table[1][0];
        }
        bool sign = true;
        double result = 0;

        for (int column = 0; column < A.n; ++column) {
            Matrix matrix(A.n-1, A.n-1);
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
                result += A.table[0][column] * findDeterminant(matrix);
                sign = false;
            } else {
                result -= A.table[0][column] * findDeterminant(matrix);
                sign = true;
            }

        }
        return result;
    }

    Matrix& operator= (Matrix const& matrix) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                this->table[i][j] = matrix.table[i][j];
            }
        }
        return *this;
    }
    friend Matrix operator+ (Matrix& first, Matrix& second) {
        Matrix res(first.n, first.m);
        for (int i = 0; i < first.n; ++i) {
            for (int j = 0; j < first.m; ++j) {
                res.table[i][j] = first.table[i][j] + second.table[i][j];
            }
        }
        return res;
    }
    friend Matrix operator- (Matrix& first, Matrix& second) {
        Matrix res(first.n, first.m);
        for (int i = 0; i < first.n; ++i) {
            for (int j = 0; j < first.m; ++j) {
                res.table[i][j] = first.table[i][j] - second.table[i][j];
            }
        }
        return res;
    }

    friend Matrix operator* (double number, Matrix& matrix) {
        Matrix res(matrix.n, matrix.m);

        for (int i = 0; i < res.n; ++i) {
            for (int j = 0; j < res.m; ++j) {
                res.table[i][j] = number * matrix.table[i][j];
            }
        }
        return res;
    }

    friend Matrix operator* (Matrix& first, Matrix& second) {
        Matrix res(first.n, second.m);

        for (int i = 0; i < res.n; ++i) {
            for (int j = 0; j < res.m; ++j) {
                res.table[i][j] = 0;
                for (int k = 0; k < first.m; ++k) {
                    res.table[i][j] += first.table[i][k] * second.table[k][j];
                }
            }
        }
        return res;
    }

    friend bool operator== (Matrix& first, Matrix& second) {
        if (first.n != second.n || first.m != second.m) {
            return false;
        }

        for (int i = 0; i < first.n; ++i) {
            for (int j = 0; j < first.m; ++j) {
                if (fabs(first.table[i][j] - second.table[i][j]) > EPS) return false;
            }
        }
        return true;
    }

    friend bool operator!= (Matrix& first, Matrix& second) {
        return !(first == second);
    }

    Matrix findTransposeMatrix() {
        Matrix transposed(m, n);

        for (int i = 0; i < transposed.n; ++i) {
            for (int j = 0; j < transposed.m; ++j) {
                transposed.table[i][j] = this->table[j][i];
            }
        }
        return transposed;
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
    friend std::istream& operator>>(std::istream &in, SquareMatrix &matrix) {
        auto* ptr = (Matrix *) &matrix;
        in >> *ptr;
        return in;
    }
    friend std::ostream& operator<<(std::ostream &out, SquareMatrix &matrix) {
        auto* ptr = (Matrix *) &matrix;
        out << *ptr;
        return out;
    }
    void resize(int newSize) {
        std::vector<double> a;
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
    friend std::ostream& operator<<(std::ostream &out, IdentityMatrix &matrix) {
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

    friend Matrix operator* (IdentityMatrix& a, ColumnVector& b) {
        auto* ptrA = (Matrix *) &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA * *ptrB;
        return * (Matrix *) &res;
    }

    friend Matrix operator* (Matrix& a, ColumnVector& b) {
        auto* ptrA = &a;
        auto* ptrB = (Matrix *) &b;
        auto res = *ptrA * *ptrB;
        return * (Matrix *) &res;
    }

    friend std::istream& operator>>(std::istream &in, ColumnVector &matrix) {
        auto* ptr = (Matrix *) &matrix;
        in >> *ptr;
        return in;
    }
    friend std::ostream& operator<<(std::ostream &out, ColumnVector &matrix) {
        auto* ptr = (Matrix *) &matrix;
        out << *ptr;
        return out;
    }
    void resize(int newSize) {
        std::vector<double> a;
        a.push_back(0);
        for (int i = n; i < newSize; ++i) {
            table.push_back(a);
        }
        n = newSize;
    }
};


#endif //OPTIMIZATION_MATRIX_H
