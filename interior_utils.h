#ifndef F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H
#define F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H

//TODO: write all functions

#include <iostream>
#include <optional>

#include "Matrix.h"
#include "Interior.h"

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

[[nodiscard]] std::optional<std::tuple<Matrix, Matrix, Matrix, int, int, bool, int>> read_IP() {
    // Get the number from variables
    std::cout << "How many variables are in the task?\n";
    int number_of_vars;
    std::cin >> number_of_vars;

    // Get the number of constraints
    std::cout << "How many constraints are in the task (excluding x >= 0)?\n";
    int number_of_equations;
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
    Matrix c(number_of_vars + number_of_equations, 1);
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

    Matrix b(number_of_equations, 1);
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
            a.table[i][number_of_vars] = 1;
            number_of_vars++;
        } else if (sign == "<=") {
            a.table[i][number_of_vars] = -1;
            number_of_vars++;
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
                    number_of_vars,
                    number_of_equations,
                    is_max_problem, eps
            )
    );
}

[[nodiscard]] std::optional<Interior> perform_interior_method() {
    auto matrix_opt = read_IP();

    auto [c, a, b, number_of_vars, number_of_equations, is_max_problem, eps] = matrix_opt.value();

    std::cout << c << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << number_of_vars << std::endl;
    std::cout << number_of_equations << std::endl;
}


bool is_feasible(std::vector<double> x) {
    for (int i = 0; i < x.size(); i++) {
        if (x[i] <= 0) {
            return false;
        }
    }
    return true;
}

std::vector<double> set_initial_solution(
        std::vector<double> c,
        Matrix &A,
        std::vector<double> b) {
    std::vector<double> x;
    for (int i = 0; i < c.size(); i++) {
        if (c[i] != 0) {
            x.push_back(rand());
        } else x.push_back(0);
    }

    int k = 0;

    while (!is_feasible(x)) {
        for (int i = 0; i < c.size(); i++) {
            if (x[i] == 0) {
                int sum = 0;
                for (int j = 0; j < A.n; j++) {
                    sum += A.table[i][j];
                }
                x[i] = c[k++] - sum;
            }
        }
    }
    return x;
}


#endif //F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H
