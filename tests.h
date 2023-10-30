#ifndef OPTIMIZATION_TESTS_H
#define OPTIMIZATION_TESTS_H

#include <fstream>
#include <string>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "json.hpp"
#include "interior_utils.h"

using json = nlohmann::json;

const double CMP_PRECISION = 0.5;

bool check_double_vector_equality(
        const std::vector<double>& first,
        const std::vector<double>& second
) {
    for (std::size_t i = 0; i < first.size(); ++i)
        if (std::abs(first[i] - second[i]) > CMP_PRECISION)
            return false;
    return true;
}

inline auto make_test_data(const char* const test_sample_path) {
    std::ifstream in(test_sample_path);
    const auto data = json::parse(in);

    slack_vars = 0;
    number_of_vars = data["number_of_vars"].get<int>();
    number_of_equations = data["number_of_equations"].get<int>();
    const auto problem = data["problem"].get<std::string>();

    auto z_row = data["z_row"].get<std::vector<double>>();
    auto table = data["matrix"].get<std::vector<std::vector<double>>>();
    auto signs = data["signs"].get<std::vector<std::string>>();
    auto rhs = data["rhs"].get<std::vector<double>>();

    ColumnVector c(number_of_vars + number_of_equations);

    for (int i = 0; i < number_of_vars; ++i)
        c.table[i][0] = z_row[i];

    const bool is_min_problem = problem == "min";

    if (is_min_problem)
        swap_to_min_problem(c);

    Matrix a(number_of_equations + 1, number_of_vars + number_of_equations + 2);
    ColumnVector b(number_of_equations);

    for (int i = 0; i < number_of_equations; ++i) {
        for (int j = 0; j < number_of_vars; ++j)
            a.table[i][j] = table[i][j];

        auto& sign = signs[i];

        if (sign == ">=") {
            a.table[i][number_of_vars + slack_vars] = -1;
            slack_vars++;
        } else if (sign == "<=") {
            a.table[i][number_of_vars + slack_vars] = 1;
            slack_vars++;
        }

        b.table[i][0] = rhs[i];
    }

    c.resize(number_of_vars + slack_vars);
    a.resize(number_of_equations, number_of_vars + slack_vars);

    return std::make_tuple(c, a, b, is_min_problem);
}

inline void perform_solution_test(
        const char* const test_sample_path,
        const double expected_z,
        std::vector<double>&& expected_vars
) {
    auto [c, A, b, is_min_problem] = make_test_data(test_sample_path);
    auto ans = calculate_answer(c, A, b, is_min_problem, 2, false).value();
    
    ASSERT_TRUE(std::abs(ans.z - expected_z) < CMP_PRECISION);
    ASSERT_TRUE(check_double_vector_equality(ans.variables, expected_vars));
}

inline void perform_no_solution_test(const char* const test_sample_path) {
    auto [c, A, b, is_min_problem] = make_test_data(test_sample_path);
    ASSERT_FALSE(calculate_answer(c, A, b, is_min_problem, 2, false).has_value());
}

TEST(InteriorTests, Test_1) {
    perform_solution_test(
            "../test_samples/1.json",
            -68,
            std::vector<double>{0, 0, 5.5, 35}
    );
}

TEST(InteriorTests, Test_2) {
    perform_solution_test(
            "../test_samples/2.json",
            -6,
            std::vector<double>{0, 0, 6}
    );
}

TEST(InteriorTests, Test_3) {
    perform_solution_test(
            "../test_samples/3.json",
            19.615,
            std::vector<double>{1.154, 7.308}
    );
}

TEST(InteriorTests, Test_4) {
    perform_solution_test(
            "../test_samples/4.json",
            -20.9,
            std::vector<double>{11.82, 0.9}
    );
}

/*TEST(InteriorTests, Test_5) {
    perform_solution_test(
            "../test_samples/5.json",
            34.96,
            std::vector<double>{0.38, 0, 0, 1.63, 2.33}
    );
}

TEST(InteriorTests, Test_6) {
    perform_solution_test(
            "../test_samples/6.json",
            -6.19,
            std::vector<double>{0, 0, 1.24, 0, 0}
    );
}

TEST(InteriorTests, Test_7) {
    perform_no_solution_test("../test_samples/7.json");
}*/

TEST(InteriorTests, Test_8) {
    perform_solution_test(
            "../test_samples/8.json",
            845.25,
            std::vector<double>{45, 52.5}
    );
}

int run_tests(int argc, char** const argv) {
    testing::InitGoogleTest(&argc, argv);
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif //OPTIMIZATION_TESTS_H
