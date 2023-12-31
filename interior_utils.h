#ifndef F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H
#define F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H

#include <iostream>
#include <optional>
#include <random>

#include "Matrix.h"
#include "Interior.h"
#include "inverse_utils.h"

int slack_vars = 0;
int number_of_vars = 0;
int number_of_equations = 0;

/** Multiply all elements of vector by -1 */

void swap_to_min_problem(Matrix& c) {
    for (auto& i : c.table)
        i[0] = -1 * i[0];
}

/** Output error if the method is not applicable */

void impossible_case(const std::string& msg) {
    std::cout << "The method is not applicable!\n";
    std::cout << msg << std::endl;
}

/** Output error if the method is not applicable */

void impossible_case(std::string&& msg) {
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
        impossible_case(std::string("Only 'max' and 'min'!"));
        return std::nullopt;
    }

    // Enter stage of coefficients for z-function
    std::cout << "Enter coefficients c(i-th) in z function (0 if absent)\n";

    ColumnVector c(number_of_vars + number_of_equations);

    for (int i = 0; i < number_of_vars; ++i) {
        std::cout << "c" << i + 1 << "=";
        std::cin >> c.table[i][0];
    }

    const bool is_min_problem = problem == "min";

    // Multiply all coefficients by -1 if it is a min problem
    if (is_min_problem)
        swap_to_min_problem(c);

    // Matrix A: each row determines coefficients for each equation
    Matrix a(number_of_equations, number_of_vars+number_of_equations);

    // Matrix B: RHS
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
        } else if (sign != "=") {
            impossible_case("Invalid sign");
            return std::nullopt;
        }

        // Input stage of RHS
        std::cout << "Enter right hand side of inequality\n";
        double RHS;
        std::cin >> RHS;

        b.table[i][0] = RHS;
    }

    c.resize(number_of_vars + slack_vars);
    a.resize(number_of_equations, number_of_vars + slack_vars);

    std::cout << "Enter approximation accuracy (the number of values after the floating point, e.g., 2)\n";

    int eps;
    std::cin >> eps;

    if (eps < 0) {
        impossible_case(std::string("Approximation accuracy is less than zero"));
        return std::nullopt;
    }

    std::cout << std::endl;

    return std::make_optional(
            std::make_tuple(
                    c,
                    a,
                    b,
                    is_min_problem, eps
            )
    );
}

/** Checks if the x is a feasible solution */

bool is_feasible(const ColumnVector& x) {
    for (int i = 0; i < x.n; i++)
        if (x.table[i][0] < 0)
            return false;

    for (int i = 0; i < number_of_vars; i++)
        if (x.table[i][0] == 0)
            return false;

    return true;
}

/** Sets an initial trial solution */

std::optional<ColumnVector> set_initial_solution(
        const ColumnVector& c,
        const Matrix &A,
        const ColumnVector& b
) {
    ColumnVector x(number_of_vars + slack_vars);
    double maxB = 0;

    for (int i = 0; i < b.n; ++i)
        if (maxB < b.table[i][0])
            maxB = b.table[i][0];

    for (int counter = 0; !is_feasible(x); ++counter) {
        if (counter == 1000000)
            return std::nullopt;

        std::random_device rd;
        std::uniform_real_distribution<double> dist(0, maxB);

        for (int i = 0; i < c.n; i++) {
            if (c.table[i][0] != 0) x.table[i][0] = round(dist(rd) * 100)/100;
            else x.table[i][0] = 0;
        }

        const int number_of_equalities = number_of_equations - slack_vars;

        for (int i = 0; i < number_of_equations; ++i) {
            double tempSum = 0;

            for (int j = 0; j < number_of_vars - number_of_equalities; ++j)
                tempSum += A.table[i][j] * x.table[j][0];

            x.table[i + number_of_vars - number_of_equalities][0] =
                    (b.table[i][0] - tempSum) / A.table[i][i + number_of_vars - number_of_equalities];
        }
    }

    return std::make_optional(std::move(x));
}

/** Calculates x_tilda for the iteration */

Matrix calculateX_tilda(const double alpha, const Matrix& c_p) {
    Matrix ones(number_of_vars + slack_vars, 1);

    for (int i = 0; i < ones.n; i++)
        ones.table[i][0] = 1;

    double v = 0;

    for (int i = 0; i < c_p.n; i++) {
        if (c_p.table[i][0] >= 0) continue;
        v = std::max(fabs(c_p.table[i][0]), v);
    }

    Matrix x_tilda = (alpha / v) * c_p;
    x_tilda = ones + x_tilda;
    return x_tilda;
}

/** Compute an optimal solution via Interior-point Algorithm */

Matrix interior_main(double alpha, const Matrix& A, Matrix& D, const ColumnVector& c) {
    IdentityMatrix I(number_of_vars + slack_vars);

    Matrix x(D.n, 1);
    Matrix x_tilda(D.n, 1);

    while (true) {
        auto A_tilda = A * D;
        auto c_tilda = D * c;
        auto A_tilda_t = A_tilda.find_transpose_matrix();

        auto arg = A_tilda * A_tilda_t;
        auto gauss = find_inverse_by_gauss(*(SquareMatrix *) &arg);

        auto P = I - A_tilda_t * gauss * A_tilda;
        auto c_p = P * c_tilda;

        x_tilda = calculateX_tilda(alpha, c_p);

        auto x_new = D * x_tilda;
        if (x == x_new) break;

        x = x_new;
        D = I;

        for (int i = 0; i < x.n; i++)
            D.table[i][i] = x.table[i][0];
    }

    std::cout << "max/min solution:\n" << x << std::endl;
    return x;
}

/** Perform Interior-point algorithm */

[[nodiscard]] std::optional<Interior> perform_interior_method() {
    // Read user input
    auto matrix_opt = read_IP();

    if (!matrix_opt.has_value())
        return std::nullopt;

    auto [c, A, b, is_min_problem, eps] = matrix_opt.value();

    std::cout << "C:\n" << c << std::endl;
    std::cout << "A:\n" << A << std::endl;
    std::cout << "B:\n" << b << std::endl;
    std::cout << "Vars: " << number_of_vars << std::endl;
    std::cout << "Slack: " << slack_vars << std::endl;
    std::cout << "Number of equations: " << number_of_equations << std::endl;

    IdentityMatrix I(number_of_vars + slack_vars);

    // Get initial trial solution
    auto init_res = set_initial_solution(c, A, b);
    if (!init_res.has_value()) {
        impossible_case(std::string("No solution"));
        return std::nullopt;
    }

    auto init = init_res.value();
    auto D = I;

    // Fill matrix D
    for (int i = 0; i < init.n; i++)
        D.table[i][i] = init.table[i][0];

    std::cout << "D:\n" << D << std::endl;

    // Compute an optimal solution
    std::cout << "Alpha is 0.5\n";
    auto x1 = interior_main(0.5, A, D, c);

    std::cout << "Alpha is 0.9\n";
    auto x2 = interior_main(0.9, A, D, c);

    Interior ans(number_of_vars, eps);

    double result = 0;
    bool is_zero = true;

    for (int i = 0; i < number_of_vars; ++i) {
        result += x1.table[i][0] * c.table[i][0];
        ans.variables[i] = x1.table[i][0];
        if (ans.variables[i] != 0) is_zero = false;
    }

    ans.z = result;

    if (is_min_problem) ans.z *= -1;
    if (is_zero) return std::nullopt;
    return std::make_optional(std::move(ans));
}

#endif //F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H
