// Copyright 2020 ICLUE @ UIUC. All rights reserved.
#include <saturation/inequalities.h>
#include <nlnum/nlnum.h>
#include <combinations.hpp>
#include <permutations.hpp>
#include <powerset.hpp>
#include <prettyprint.hpp>
#include <product.hpp>
#include <range.hpp>

#include <algorithm>
#include <climits>
#include <deque>
#include <iterator>
#include <iostream>
#include <numeric>
#include <ostream>
#include <stack>
#include <utility>
#include <vector>

namespace saturation {

using nlnum::Partition;

Partition Tau(const Set& I) {
  const size_t r = I.size();

  Partition ss{I.begin(), I.end()};
  std::sort(ss.begin(), ss.end());

  const auto R = iter::range(1, static_cast<int32_t>(r)+1);
  const Partition rr{R.begin(), R.end()};

  Partition tt;
  for (size_t idx = 0; idx < r; ++idx) {
    tt.push_back(ss[idx]-rr[idx]);
  }
  std::reverse(tt.begin(), tt.end());

  return tt;
}

Set Complement(const Set& I, const Int n) {
  for (const auto& i : I) {
    if (i > 4*n) throw std::invalid_argument("I must be a subset of [4n].");
  }

  Set cI;
  std::transform(
      I.begin(), I.end(),
      std::inserter(cI, cI.begin()),
      [&](const Int& i) -> Int { return 4*n + 1 - i; });

  return cI;
}

// Returns true iff lam is completely inside of mu.
bool IsInside(const Partition& lam, const Partition& mu) {
  if (lam.size() > mu.size()) return false;
  for (size_t i = 0; i < lam.size(); ++i) {
    if (lam[i] > mu[i]) return false;
  }
  return true;
}

Partition Check(const Partition& lam, const Int a, const Int b) {
  nlnum::ValidatePartitions({ lam });

  Partition mu(b, a);
  if (!IsInside(lam, mu)) {
    throw std::invalid_argument("lam must be inside the rectangle (a^b).");
  }

  const Int k = lam.size();
  const auto diff = static_cast<int32_t>(b >= k ? b-k : 0);
  std::transform(
      lam.rbegin(), lam.rend(),
      mu.cbegin() + diff, mu.begin() + diff,
      [&](const Int& x, const Int& y) -> Int { return y - x; });

  while (mu[mu.size()-1] == 0) {
    mu.pop_back();
  }
  return mu;
}

// Throw an exception if it is *not* the case that X \subseteq Y \subseteq [4n].
void CheckContainment(const Set& X, const Set& Y, const Int n) {
  for (const auto& y : Y) {
    if (y < 1 || y > 4*n) {
      throw std::invalid_argument("Y must be a subset of [4n].");
    }
  }

  for (const auto& x : X) {
    if (Y.find(x) == Y.end()) {
      throw std::invalid_argument("X must be a subset of Y.");
    }
  }
}

Set Chi(const Set& X, const Set& Y, const Int n, const Int b) {
  CheckContainment(X, Y, n);
  if (b != 0 && b != 2) {
    throw std::invalid_argument("b must be either 0 or 2.");
  }

  Set chi;
  Int numElems = 0;

  // O(n*lg n).
  for (size_t i = 1; i <= 4*n; ++i) {
    if (X.find(i) != X.end()) {
      ++numElems;
      chi.insert(numElems);
    } else if (Y.find(i) != Y.end()) {
      if (b == 0) ++numElems;
    } else {
      if (b == 2) ++numElems;
    }
  }

  return chi;
}

std::vector<std::pair<Set, Set>> Disjoints(const Int n, const Int r) {
  const auto R = iter::range(1, static_cast<int32_t>(4*n+1));
  const auto RR = Set{R.begin(), R.end()};
  std::vector<std::pair<Set, Set>> ans;

  for (const auto& X : iter::combinations(R, r)) {
    Set Z;
    std::set_difference(
      RR.begin(), RR.end(),
      X.begin(), X.end(),
      std::inserter(Z, Z.begin()));

    for (const auto& Y : iter::powerset(Z)) {
      const auto XX = Set(X.begin(), X.end());
      const auto YY = Set(Y.begin(), Y.end());
      ans.emplace_back(XX, YY);
    }
  }

  return ans;
}


}  // namespace saturation
