// Copyright 2020 ICLUE @ UIUC. All rights reserved.
#include <saturation/inequalities.h>
#include <nlnum/nlnum.h>
#include <combinations.hpp>
#include <prettyprint.hpp>
#include <product.hpp>
#include <range.hpp>

#include <algorithm>
#include <cstdint>
#include <ostream>
#include <utility>
#include <vector>

namespace saturation {

using nlnum::Partition;

std::ostream& operator<<(std::ostream& os, const Triple& t) {
  os << "Triple{"
     << "la=" << t.la << ", "
     << "mu=" << t.mu << ", "
     << "nu=" << t.nu << "}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Sets& s) {
  os << "Sets{"
     << "I=" << s.I << ", "
     << "J=" << s.J << ", "
     << "K=" << s.K << ", "
     << "bI=" << s.bI << ", "
     << "bJ=" << s.bJ << ", "
     << "bK=" << s.bK << "}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const CounterExample& c) {
  return os << "CounterExample{[" << c.triple << "], [" << c.sets << "]}";
}

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

Set Bar(const Set& I, const Int n) {
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

Set Complement(const Set& X, const Int n) {
  const auto R = iter::range(1, static_cast<int32_t>(4*n+1));
  const auto RR = Set{R.begin(), R.end()};

  Set Z;
  std::set_difference(
      RR.begin(), RR.end(),
      X.begin(), X.end(),
      std::inserter(Z, Z.begin()));
  return Z;
}

std::vector<std::pair<Set, Set>> Disjoints(const Int n, const Int r) {
  const auto R = iter::range(1, static_cast<int32_t>(4*n+1));
  const auto RR = Set{R.begin(), R.end()};
  std::vector<std::pair<Set, Set>> ans;

  for (const auto& X : iter::combinations(R, r)) {
    const auto XX = Set(X.begin(), X.end());
    const auto bX = Bar(XX, n);
    Set Z;
    std::set_intersection(
        XX.begin(), XX.end(), bX.begin(), bX.end(),
                          std::inserter(Z, Z.begin()));

    if (Z.size() == 0) {
      ans.emplace_back(XX, bX);
    }
  }

  return ans;
}

bool RestrictedLengths(const Int n, const Int r, const Set& I, const Set& J,
                       const Set& K) {
  // |I\cap [2n]| + |J\cap[2n]| + |K\cap[2n]| = r
  // numAtMost = |I\cap [2n]| + |J\cap[2n]| + |K\cap[2n]|.
  Int numAtMost = 0;
  for (const auto& X : {I, J, K}) {
    for (const auto& x : X) {
      if (x <= 2 * n) ++numAtMost;
    }
  }
  return numAtMost == r;
}

bool Consecutive(const Int n, const Set& I, const Set& J, const Set& K) {
  for (const auto& X : {I, J, K}) {
    std::set<Int> sX;
    for (const auto& x : X) {
      if (n + 1 <= x && x <= 3 * n) sX.insert(x);
    }
    Int count = 0;
    // We know all the numbers in I,J,K are strictly greater than 0.
    for (Int k = 3 * n; k > 0; --k) {
      if (sX.find(k) == sX.end()) break;
      ++count;
    }
    if (count != sX.size()) return false;
  }
  return true;
}

bool GoodLrCoefs(const Int n, const Int r, const Set& I, const Set& J,
                 const Set& K, const Set& bI, const Set& bJ, const Set& bK) {
  // lrcoefs = 1.
  const auto bIc = Complement(bI, n);
  const auto bJc = Complement(bJ, n);
  const auto bKc = Complement(bK, n);

  const auto c1 = nlnum::lrcoef(
      Tau(Chi(K, bKc, n, 0)), Check(Tau(Chi(I, bIc, n, 0)), 4 * n - 2 * r, r),
      Check(Tau(Chi(J, bJc, n, 0)), 4 * n - 2 * r, r));
  if (c1 != 1) return false;

  const auto c2 =
      nlnum::lrcoef(Tau(Chi(K, bKc, n, 2)), Check(Tau(Chi(I, bIc, n, 2)), r, r),
                    Check(Tau(Chi(J, bJc, n, 2)), r, r));
  if (c2 != 1) return false;

  return true;
}

bool IsGood(const Int n, const Int r, const Set& I, const Set& J, const Set& K,
            const Set& bI, const Set& bJ, const Set& bK, Sets* s) {
  // Three conditions from the paper Shiliang sent.
  if (!RestrictedLengths(n, r, I, J, K)) return false;
  if (!Consecutive(n, I, J, K)) return false;
  if (!GoodLrCoefs(n, r, I, J, K, bI, bJ, bK)) return false;
  if (s != nullptr) {
    *s = {I, J, K, bI, bJ, bK};
  }
  return true;
}

std::vector<Sets> SatIneqs(const Int n, const Int r) {
  const auto djs = Disjoints(n, r);
  const auto RR = iter::range(1, static_cast<int32_t>(n)+1);
  const std::vector<int32_t> Rv{RR.begin(), RR.end()};
  const std::vector<Int> R{Rv.begin(), Rv.end()};
  const Set S{R.begin(), R.end()};

  std::vector<Sets> ans{};
#pragma omp parallel for schedule(dynamic)
  for (auto ita = djs.begin(); ita < djs.end(); ++ita) {
    const auto& Iq = *ita;
    const Set& I = Iq.first;
    const Set& bI = Iq.second;
    for (const auto& Jq : djs) {
      const Set& J = Jq.first;
      const Set& bJ = Jq.second;
      for (const auto& Kq : djs) {
        const Set& K = Kq.first;
        const Set& bK = Kq.second;

        Sets s;
        if (IsGood(n, r, I, J, K, bI, bJ, bK, &s)) {
#pragma omp critical
          ans.push_back(s);
        }
      }
    }
  }

  return ans;
}

// Computes the sum of pi on the indices X \cap [2n].
int64_t Sum(const nlnum::Partition& pi, const Set& X, const Int n) {
  Int ans = 0;
  for (Int ii : X) {
    if (ii > pi.size() || ii > 2 * n) continue;
    ans += pi[ii - 1];
  }
  return static_cast<int64_t>(ans);
}

bool Satisfies(const Triple& triple, const std::vector<Sets>& satIneqs,
               const Int n, Sets* sets) {
  const auto la = triple.la;
  const auto mu = triple.mu;
  const auto nu = triple.nu;

  for (const auto& ineq : satIneqs) {
    const auto lhs =
        (Sum(la, ineq.bI, n) - Sum(la, ineq.I, n) + Sum(mu, ineq.bJ, n) -
         Sum(mu, ineq.J, n) + Sum(nu, ineq.bK, n) - Sum(nu, ineq.K, n));

    if (lhs < 0) {
      if (sets != nullptr) {
        *sets = ineq;
      }
      return false;
    }
  }

  return true;
}

std::vector<CounterExample> Flagger(const Int n, const Int r) {
  const auto get = [](const nlnum::PartitionsIn& pi) -> std::vector<Partition> {
    std::vector<Partition> pars;
    for (const auto& par : pi) {
      pars.push_back(par);
    }
    return pars;
  };
  const auto satIneqs = SatIneqs(n, r);
  const auto iPars = nlnum::PartitionsIn(Partition(n, n), n);

  const auto ins = iter::product<3>(iter::range(1, static_cast<int32_t>(n)+1));
  std::vector<std::tuple<Int, Int, Int>> ns;
  for (const auto& t : ins) {
    const auto kl = std::get<0>(t);
    const auto km = std::get<1>(t);
    const auto kn = std::get<2>(t);
    if ((kl + km + kn) % 2 == 1) continue;
    ns.emplace_back(kl, km, kn);
  }

  std::vector<CounterExample> ans;
#pragma omp parallel for schedule(dynamic) shared(ins)
  for (auto it = ns.begin(); it < ns.end(); ++it) {
    const auto kl = std::get<0>(*it);
    const auto km = std::get<1>(*it);
    const auto kn = std::get<2>(*it);

    const auto nls = get(nlnum::PartitionsIn(Partition(kl, kl), kl));
    const auto nms = get(nlnum::PartitionsIn(Partition(km, km), km));
    const auto nns = get(nlnum::PartitionsIn(Partition(kn, kn), kn));

    for (auto la = nls.begin(); la < nls.end(); ++la) {
      for (auto mu = nms.begin(); mu < nms.end(); ++mu) {
        for (auto nu = nns.begin(); nu < nns.end(); ++nu) {
          Sets s;
          if (!Satisfies({*la, *mu, *nu}, satIneqs, n, &s)) {
#pragma omp critical
            ans.emplace_back(Triple{*la, *mu, *nu}, s);
          }
        }
      }
    }
  }

  return ans;
}

}  // namespace saturation
