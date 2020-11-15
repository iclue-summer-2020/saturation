// Copyright (c) 2020 ICLUE. All rights reserved.

#define CATCH_CONFIG_MAIN

#include <exception>
#include <iostream>
#include <numeric>
#include <vector>

#include <catch2/catch.hpp>
#include <saturation/inequalities.h>

using nlnum::Partition;
using saturation::Int;
using saturation::Set;
using saturation::check;
using saturation::complement;
using saturation::tau;


TEST_CASE("easy tau", "[tau]") {
  const Set I = {1, 2, 3, 4, 5};
  const Partition expected = {0, 0, 0, 0, 0};
  REQUIRE(tau(I) == expected);
}

TEST_CASE("empty tau", "[tau]") {
  const Set I = {};
  const Partition expected = {};
  REQUIRE(tau(I) == expected);
}

TEST_CASE("normal tau", "[tau]") {
  const Set I = {2, 4, 7};
  const Partition expected = {4, 2, 1};
  REQUIRE(tau(I) == expected);
}

TEST_CASE("empty complement", "[complement]") {
  const Int n = 4;
  const Set I = {};
  const Set expected = {};
  REQUIRE(complement(I, n) == expected);
}

TEST_CASE("easy complement", "[complement]") {
  const Int n = 4;
  const Set I = {1, 2, 3, 4, 5};
  const Set expected = {16, 15, 14, 13, 12};
  REQUIRE(complement(I, n) == expected);
}

TEST_CASE("normal complement", "[complement]") {
  const Int n = 4;
  const Set I = {1, 4, 12};
  const Set expected = {16, 13, 5};
  REQUIRE(complement(I, n) == expected);
}

TEST_CASE("invalid argument", "[complement]") {
  const Int n = 1;
  const Set I = {1, 2, 3, 4, 5};
  REQUIRE_THROWS_AS(complement(I, n), std::invalid_argument);
}

TEST_CASE("easy check", "[check]") {
  const Int a = 3;
  const Int b = 3;
  const Partition lam = {3, 2, 1};
  const Partition expected = {2, 1};
  REQUIRE(check(lam, a, b) == expected);
}

TEST_CASE("normal check", "[check]") {
  const Int a = 5;
  const Int b = 9;
  const Partition lam = {4, 4, 4, 3, 1, 1, 1};
  const Partition expected = {5, 5, 4, 4, 4, 2, 1, 1, 1};
  REQUIRE(check(lam, a, b) == expected);
}

TEST_CASE("tricky check", "[check]") {
  const Int a = 5;
  const Int b = 3;
  const Partition lam = {5, 4};
  const Partition expected = {5, 1};
  REQUIRE(check(lam, a, b) == expected);
}

TEST_CASE("not contained check", "[check]") {
  const Int a = 3;
  const Int b = 3;
  const Partition lam = {4, 3, 2, 1};
  REQUIRE_THROWS_AS(check(lam, a, b), std::invalid_argument);
}

TEST_CASE("invalid partition", "[check]") {
  const Int a = 8;
  const Int b = 7;
  // Needs to be weakly decreasing.
  const Partition lam = {1, 2, 3, 4};
  REQUIRE_THROWS_AS(check(lam, a, b), std::invalid_argument);
}
