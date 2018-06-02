#include <algorithm>
#include <stdexcept>
#include <string>

#include "helpers.h"
#include "gtest/gtest.h"

TEST(helpers, ReadCommandLineIntValidInput) {
  EXPECT_EQ(ising::read_command_line_int("12"), 12);
}

TEST(helpers, ReadCommandLineIntInvalidInput) {
  EXPECT_THROW(ising::read_command_line_int("12.0"), std::invalid_argument);
  EXPECT_THROW(ising::read_command_line_int("abc"), std::invalid_argument);
  EXPECT_THROW(ising::read_command_line_int(""), std::invalid_argument);
}

TEST(helpers, ReadCommandLineDoubleValidInput) {
  EXPECT_EQ(ising::read_command_line_double("12"), 12.0);
  EXPECT_EQ(ising::read_command_line_double("12.0"), 12.0);
}

TEST(helpers, ReadCommandLineDoubleInvalidInput) {
  EXPECT_THROW(ising::read_command_line_double("abc"), std::invalid_argument);
  EXPECT_THROW(ising::read_command_line_double(""), std::invalid_argument);
}
