#ifndef F23_OPTIMIZATION_TASK2_INTERIOR_H
#define F23_OPTIMIZATION_TASK2_INTERIOR_H

//TODO: write a class for storing an answer for a task

#include <iostream>

#include "Matrix.h"

struct Interior {
    double z;
    int accuracy;
    std::vector<double> variables;

    Interior(const std::size_t size, const int eps) {
        z = 0;
        accuracy = eps;
        variables.resize(size, 0);
    }

    ~Interior() = default;

    friend std::ostream& operator<<(std::ostream& out, Interior& answer) {
        out << "Value of z: " << answer.z << std::endl;

        for (int i = 0; i < answer.variables.size(); ++i)
            out << std::setprecision(answer.accuracy)
                << std::fixed
                << "x" << i + 1 << " = " << answer.variables[i] << std::endl;

        out << std::endl;
        return out;
    }
};

#endif //F23_OPTIMIZATION_TASK2_INTERIOR_H
