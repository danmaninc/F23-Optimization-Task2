#ifndef F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H
#define F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H

//TODO: write all functions

#include <iostream>
#include <optional>
#include <random>

#include "Matrix.h"
#include "Interior.h"
#include "inverse_utils.h"

int slack_vars = 0;
int number_of_vars = 0;
int number_of_equations = 0;

// Multiply all elements of vector by -1
void swap_to_max_problem(Matrix c) {
    for (auto& i : c.table)
        i[0] = -1 * i[0];
}

// Output error if the method is not applicable
void impossible_case(std::string& msg) {
    std::cout << "The method is not applicable!\n";
    std::cout << msg << std::endl;
}

[[nodiscard]] std::optional<std::tuple<ColumnVector, Matrix, ColumnVector, bool, int>> read_IP() {
    // Get the number from variables
    std::cout << "How many variables are in the task?\n";
    std::cin >> number_of_vars;

    // Get the number of constraints
    std::cout << "How many constraints are in the task (excluding x >= 0)?\n";
    std::cin >> number_of_equations;

    // Determine the type of optimization problem
    std::cout << "Is it maximization (enter 'max') or minimization (enter 'min') problem?\n";
    std::string problem;
    std::cin >> problem;

    // If the optimization problem is not max or min, then output error
    if (problem != "max" && problem != "min") {
        std::string msg = "Only 'max' and 'min'!";
        impossible_case(msg);
        return std::nullopt;
    }

    // Enter stage of coefficients for z-function
    std::cout << "Enter coefficients c(i-th) in z function (0 if absent)\n";
    ColumnVector c(number_of_vars + number_of_equations);
    for (int i = 0; i < number_of_vars; ++i) {
        std::cout << "c" << i + 1 << "=";
        std::cin >> c.table[i][0];
    }

    const bool is_max_problem = problem == "max";
    // Multiply all coefficients by -1 if it is max problem
    if (is_max_problem)
        swap_to_max_problem(c);

    // Matrix A: each row determines coefficients for each equation
    Matrix a(number_of_equations, number_of_vars+number_of_equations);

    ColumnVector b(number_of_equations);
    // Input stage of coefficients for constraints
    std::cout << "Enter coefficients (0 if absent) for all constraints (left hand side)\n";
    for (int i = 0; i < number_of_equations; ++i) {
        std::cout << i + 1 << " constraint\n";

        for (int j = 0; j < number_of_vars; ++j) {
            std::cout << "c" << j + 1 << "=";
            std::cin >> a.table[i][j];
        }

        std::cout << "Enter sign of inequality\n";
        std::string sign;
        std::cin >> sign;

        // Add slack variable to A
        if (sign == ">=") {
            a.table[i][number_of_vars + slack_vars] = -1;
            slack_vars++;
        } else if (sign == "<=") {
            a.table[i][number_of_vars + slack_vars] = 1;
            slack_vars++;
        }


        std::cout << "Enter right hand side of inequality\n";

        double RHS;
        std::cin >> RHS;

        if (RHS < 0) {
            std::string msg = std::to_string(RHS) + " is less than zero";
            impossible_case(msg);
            return std::nullopt;
        }

        b.table[i][0] = RHS;
    }

    std::cout << "Enter approximation accuracy (the number of values after the floating point, e.g., 2)\n";

    int eps;
    std::cin >> eps;

    if (eps < 0) {
        impossible_case((std::string &) "Approximation accuracy is less than zero");
        return std::nullopt;
    }

    //set_basic_vars(c, number_of_vars);
    //set_presentation(matrix, number_of_vars, number_of_equations);
    std::cout << std::endl;

    return std::make_optional(
            std::make_tuple(
                    c,
                    a,
                    b,
                    is_max_problem, eps
            )
    );
}

bool is_feasible(ColumnVector& x) {
    for (int i = 0; i < x.n; i++) {
        if (x.table[i][0] <= 0) {
            return false;
        }
    }
    return true;
}

ColumnVector set_initial_solution(
        ColumnVector c,
        Matrix &A,
        ColumnVector b) {
    ColumnVector x(number_of_vars + slack_vars);

    double maxB = 0;
    for (int i = 0; i < b.n; ++i) {
        if (maxB < b.table[i][0]) maxB = b.table[i][0];
    }
    std::cout << maxB << std::endl;
    while (!is_feasible(x)) {
        std::random_device rd;
        std::uniform_real_distribution<double> dist(0, maxB);
        for (int i = 0; i < c.n; i++)
        {
            if (c.table[i][0] != 0) x.table[i][0] = dist(rd);
            else x.table[i][0] = 0;
        }

        for (int i = 0; i < A.n; ++i) {
            double tempSum = 0;
            for (int j = 0; j < number_of_vars; j++) {
                tempSum += A.table[i][j] * x.table[j][0];
            }
            x.table[i + number_of_vars][0] = (b.table[i][0] - tempSum) * A.table[i][i+number_of_vars];
        }
    }
    return x;
}

auto calculateX_tilda(double alpha, Matrix& c_p)
{
    Matrix ones(number_of_vars + slack_vars, 1);
    for (int i = 0; i < ones.n; i++)
    {
        ones.table[i][0] = 1;
    }

    double v = 0;
    for (int i = 0; i < c_p.n; i++)
    {
        if (c_p.table[i][0] >= 0) continue;

        v = std::max(fabs(c_p.table[i][0]), v);
    }


    Matrix x_tilda = (alpha/v) * c_p;
    x_tilda = ones + x_tilda;

    return x_tilda;
}

[[nodiscard]] std::optional<Interior> perform_interior_method() {
    auto matrix_opt = read_IP();

    auto [c, A, b, is_max_problem, eps] = matrix_opt.value();

    std::cout << c << std::endl;
    std::cout << A << std::endl;
    std::cout << b << std::endl;
    std::cout << number_of_vars << std::endl;
    std::cout << number_of_equations << std::endl;

    IdentityMatrix I(number_of_vars + slack_vars);

    auto init = set_initial_solution(c, A, b);

    Matrix D = I;
    for (int i = 0; i < init.n; i++)
    {
        D.table[i][i] = init.table[i][0];
    }

    std::cout << D << std::endl;

    //Matrix* x;
    Matrix x(D.n, 1);
    Matrix x_tilda(D.n, 1);
    while (true)
    {
        Matrix A_tilda = A * D;
        Matrix c_tilda = D * c;
        Matrix A_tilda_t = A_tilda.findTransposeMatrix();

        Matrix arg = A_tilda * A_tilda_t;
        Matrix gauss = findInverseByGauss(*(SquareMatrix*)&arg);

        Matrix P = A_tilda_t * gauss;
        P = P * A_tilda;
        P = I - P;

        Matrix c_p = P * c_tilda;

        // alpha = 0.5
        x_tilda = calculateX_tilda(0.5, c_p);

        Matrix x_new = D * x_tilda;
        std::cout << "Decision bebra:\n" << x_tilda << std::endl;

        std::cout << "max/min val:\n" << x_new << std::endl;

        if (x == x_new) break;

        x = x_new;
        D = I;
        for (int i = 0; i < x.n; i++)
        {
            D.table[i][i] = x.table[i][0];
        }
    }

    std::cout << "Bebraversity 0.5" << std::endl;

    std::cout << "Decision variables:\n" << x_tilda << std::endl;

    // TODO: min/max
    std::cout << "max/min val:\n" << x << std::endl;

    Interior ans(2, eps);
    return ans;
}

#endif //F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H
