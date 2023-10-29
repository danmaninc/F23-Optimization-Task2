#include <iostream>

#include "Interior.h"
#include "interior_utils.h"

int main(int argc, char** argv) {

    auto answer = perform_interior_method();

    if (answer.has_value())
        std::cout << answer.value();
    else
        std::cout << "No solution" << std::endl;
    return 0;
}