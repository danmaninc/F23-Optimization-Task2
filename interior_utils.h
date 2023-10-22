#ifndef F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H
#define F23_OPTIMIZATION_TASK2_INTERIOR_UTILS_H

//TODO: write all functions

#include <iostream>
#include <optional>

#include "Matrix.h"
#include "Interior.h"

[[nodiscard]] std::optional<Interior> perform_interior_method() {

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
