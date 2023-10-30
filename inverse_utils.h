#ifndef F23_OPTIMIZATION_TASK2_INVERSE_UTILS_H
#define F23_OPTIMIZATION_TASK2_INVERSE_UTILS_H

#include <iomanip>
#include <cmath>

#include "Matrix.h"

struct AugmentedMatrix : public Matrix {
    explicit AugmentedMatrix(const SquareMatrix& A) : Matrix(A.n, 2 * A.m) {
        IdentityMatrix identity(A.n);

        for (int i = 0; i < this->n; ++i) {
            for (int j = 0; j < A.m; ++j) {
                this->table[i][j] = A.table[i][j];
            }
            int index = A.m;
            for (int j = 0; j < A.m; ++j) {
                this->table[i][index] = identity.table[i][j];
                index++;
            }
        }
    }

    AugmentedMatrix(const SquareMatrix& A, const ColumnVector& b) : Matrix(A.n, A.m + 1) {
        for (int i = 0; i < this->n; ++i) {
            for (int j = 0; j < A.m; ++j)
                this->table[i][j] = A.table[i][j];
            this->table[i][A.m] = b.table[i][0];
        }
    }

    ~AugmentedMatrix() = default;

    SquareMatrix get_result() {
        SquareMatrix res(this->n);

        for (int i = 0; i < n; ++i) {
            int index = n;

            for (int j = 0; j < n; ++j) {
                res.table[i][j] = this->table[i][index];
                index++;
            }
        }

        return res;
    }

    ColumnVector get_vector() {
        ColumnVector res(this->n);

        for (int i = 0; i < n; ++i)
            res.table[i][0] = this->table[i][n];

        return res;
    }

    friend std::ostream& operator<<(std::ostream& out, const AugmentedMatrix& matrix) {
        for (auto& i : matrix.table) {
            for (int j = 0; j < i.size() - 1; ++j) {
                if (std::fabs(i[j]) < MIN_COUT) {
                    out << std::setprecision(NUMBERS_AFTER_DOT) << std::fixed << ZERO;
                } else {
                    out << std::setprecision(NUMBERS_AFTER_DOT) << std::fixed << i[j];
                }

                if (j != i.size() - 1)
                    out << " ";
            }

            out << std::endl;
        }

        for (int i = 0; i < matrix.n; ++i) {
            if (std::fabs(matrix.table[i][matrix.n]) < MIN_COUT) {
                out << std::setprecision(NUMBERS_AFTER_DOT) << std::fixed << ZERO << std::endl;
            } else {
                out << std::setprecision(NUMBERS_AFTER_DOT) << std::fixed << matrix.table[i][matrix.n] << std::endl;
            }
        }

        return out;
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    explicit EliminationMatrix(const int n) : IdentityMatrix(n) {}
    ~EliminationMatrix() = default;

    void eliminate(Matrix& A, int row, int column, int pivot) {
        double coefficient = (-1) * A.table[row][column] / A.table[pivot][column];
        this->table[row][column] = coefficient;
    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    explicit PermutationMatrix(const int n) : IdentityMatrix(n) {}
    ~PermutationMatrix() = default;

    void permutate(int first, int second) {
        for (int i = 0; i < n; ++i) {
            double temp = this->table[first][i];
            this->table[first][i] = this->table[second][i];
            this->table[second][i] = temp;
        }
    }
};


void reduce_matrix_to_diag(Matrix& A) {
    int last = A.n - 1;
    int pivot = last;
    for (int i = last; i >= 0; --i) {
        if (A.table[pivot][i] == 0) {
            pivot--;
            continue;
        }
        for (int j = pivot - 1; j >= 0; --j) {
            if (A.table[j][i] == 0) {
                continue;
            }
            EliminationMatrix E(A.n);
            E.eliminate(A, j, i, pivot);
            A = E * A;
        }
        pivot--;
    }
}

void reduce_matrix_to_upp(Matrix& A) {
    int pivot = 0;
    for (int i = 0; i < A.n; ++i) {
        double maxValue = abs(A.table[pivot][i]);
        int maxIndex = pivot;
        for (int j = pivot + 1; j < A.n; ++j) {
            if (abs(A.table[j][i]) > maxValue) {
                maxValue = abs(A.table[j][i]);
                maxIndex = j;
            }
        }
        if (maxValue == 0) {
            pivot++;
            continue;
        }
        if (maxIndex != pivot) {
            PermutationMatrix P(A.n);
            P.permutate(maxIndex, pivot);
            A = P * A;
        }
        for (int j = pivot + 1; j < A.n; ++j) {
            if (A.table[j][i] == 0) {
                continue;
            }
            EliminationMatrix E(A.n);
            E.eliminate(A, j, i, pivot);
            A = E * A;
        }
        pivot++;
    }
}


void normalize(Matrix& A) {
    for (int i = 0; i < A.n; ++i) {
        double value = A.table[i][i];

        for (int j = i; j < A.m; ++j)
            A.table[i][j] /= value;
    }
}


SquareMatrix find_inverse_by_gauss(const SquareMatrix& A) {
    if (A.n == 1) {
        SquareMatrix res(1);
        res.table[0][0] = 1 / A.table[0][0];
        return res;
    }

    AugmentedMatrix aug(A);
    reduce_matrix_to_upp(*(Matrix *) &aug);
    reduce_matrix_to_diag(*(Matrix *) &aug);
    normalize(*(Matrix *) &aug);
    return aug.get_result();
}

#endif //F23_OPTIMIZATION_TASK2_INVERSE_UTILS_H
