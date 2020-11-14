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
