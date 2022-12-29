#include <gtest/gtest.h>

// test init constructor with 0 parameters
TEST(GraphTesting, initConstructorZeroParam) {
  // Expect two strings not to be equal.
  
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}