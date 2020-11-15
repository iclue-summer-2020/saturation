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
using saturation::Check;
using saturation::Chi;
using saturation::Complement;
using saturation::Disjoints;
using saturation::Set;
using saturation::Tau;

TEST_CASE("easy tau", "[tau]") {
  const Set I = {1, 2, 3, 4, 5};
  const Partition expected = {0, 0, 0, 0, 0};
  REQUIRE(Tau(I) == expected);
}

TEST_CASE("empty tau", "[tau]") {
  const Set I = {};
  const Partition expected = {};
  REQUIRE(Tau(I) == expected);
}

TEST_CASE("normal tau", "[tau]") {
  const Set I = {2, 4, 7};
  const Partition expected = {4, 2, 1};
  REQUIRE(Tau(I) == expected);
}

TEST_CASE("empty complement", "[complement]") {
  const Int n = 4;
  const Set I = {};
  const Set expected = {};
  REQUIRE(Complement(I, n) == expected);
}

TEST_CASE("easy complement", "[complement]") {
  const Int n = 4;
  const Set I = {1, 2, 3, 4, 5};
  const Set expected = {16, 15, 14, 13, 12};
  REQUIRE(Complement(I, n) == expected);
}

TEST_CASE("normal complement", "[complement]") {
  const Int n = 4;
  const Set I = {1, 4, 12};
  const Set expected = {16, 13, 5};
  REQUIRE(Complement(I, n) == expected);
}

TEST_CASE("invalid argument", "[complement]") {
  const Int n = 1;
  const Set I = {1, 2, 3, 4, 5};
  REQUIRE_THROWS(
      Complement(I, n),
      std::invalid_argument("I must be a subset of [4n]."));
}

TEST_CASE("easy check", "[check]") {
  const Int a = 3;
  const Int b = 3;
  const Partition lam = {3, 2, 1};
  const Partition expected = {2, 1};
  REQUIRE(Check(lam, a, b) == expected);
}

TEST_CASE("normal check", "[check]") {
  const Int a = 5;
  const Int b = 9;
  const Partition lam = {4, 4, 4, 3, 1, 1, 1};
  const Partition expected = {5, 5, 4, 4, 4, 2, 1, 1, 1};
  REQUIRE(Check(lam, a, b) == expected);
}

TEST_CASE("tricky check", "[check]") {
  const Int a = 5;
  const Int b = 3;
  const Partition lam = {5, 4};
  const Partition expected = {5, 1};
  REQUIRE(Check(lam, a, b) == expected);
}

TEST_CASE("not contained check", "[check]") {
  const Int a = 3;
  const Int b = 3;
  const Partition lam = {4, 3, 2, 1};
  REQUIRE_THROWS(
      Check(lam, a, b),
      std::invalid_argument("lam must be inside the rectangle (a^b)."));
}

TEST_CASE("invalid partition", "[check]") {
  const Int a = 8;
  const Int b = 7;
  // Needs to be weakly decreasing.
  const Partition lam = {1, 2, 3, 4};
  REQUIRE_THROWS(
      Check(lam, a, b),
      std::invalid_argument("Each partition must be weakly decreasing."));
}

TEST_CASE("easy chi", "[chi]") {
  const Int n = 2;
  const Set X = {1, 2, 3};
  const Set Y = {1, 2, 3, 4, 5};
  const auto res = Chi(X, Y, n, 2);

  REQUIRE(res.size() <= 4*n-Y.size()+X.size());
  REQUIRE(res.size() == 3);

  SECTION("containment of res in [4n-|Y|+|X|]") {
    for (const auto& e : res) {
      REQUIRE(e >= 1);
      REQUIRE(e <= 4 * n - Y.size() + X.size());
    }
  }

  REQUIRE(res == Set({1, 2, 3}));
}

TEST_CASE("tricky chi", "[chi]") {
  const Int n = 2;
  const Set X = {1, 3, 8};
  const Set Y = {1, 2, 3, 7, 8};
  const auto res = Chi(X, Y, n, 0);

  REQUIRE(res.size() <= Y.size());
  REQUIRE(res.size() == 3);

  SECTION("containment of res in [|Y|]") {
    for (const auto& e : res) {
      REQUIRE(e >= 1);
      REQUIRE(e <= Y.size());
    }
  }

  REQUIRE(res == Set({1, 3, 5}));
}

TEST_CASE("X not in Y", "[chi]") {
  const Int n = 2;
  const Set X = {1, 3, 6, 8};
  const Set Y = {1, 2, 3, 7, 8};
  REQUIRE_THROWS(
      Chi(X, Y, n, 0),
      std::invalid_argument("X must be a subset of Y."));
}

TEST_CASE("Y not in [4n]", "[chi]") {
  const Int n = 2;
  const Set X = {1, 3, 8};
  const Set Y = {1, 2, 3, 7, 8, 9};
  REQUIRE_THROWS(
      Chi(X, Y, n, 0),
      std::invalid_argument("Y must be a subset of [4n]."));
}

TEST_CASE("b is invalid", "[chi]") {
  const Int n = 2;
  const Set X = {1, 3, 8};
  const Set Y = {1, 2, 3, 7, 8};
  REQUIRE_THROWS(
      Chi(X, Y, n, 1),
      std::invalid_argument("b must be either 0 or 2."));
}

TEST_CASE("easy disjoints", "[disjoints]") {
  const Int n = 2;
  const Int r = 3;
  const auto& djs = Disjoints(n, r);

  // (4n choose r) * 2^{4n-r}.
  // So (8 choose 3) * 2^5 = 56 * 32 = 1792.
  REQUIRE(djs.size() == 1792);
}

TEST_CASE("edge case", "[disjoints]") {
  const Int n = 1;
  const Int r = 4;
  const auto& djs = Disjoints(n, r);

  // (4n choose r) * 2^{4n-r}.
  // So (4 choose 4) * 2^0 = 1.
  REQUIRE(djs.size() == 1);
}
