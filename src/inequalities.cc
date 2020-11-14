// Copyright 2020 ICLUE @ UIUC. All rights reserved.
#include <saturation/inequalities.h>
#include <nlnum/nlnum.h>
#include <combinations.hpp>
#include <permutations.hpp>
#include <prettyprint.hpp>
#include <product.hpp>
#include <range.hpp>

#include <algorithm>
#include <climits>
#include <deque>
#include <iostream>
#include <numeric>
#include <ostream>
#include <stack>
#include <utility>
#include <vector>

namespace saturation {

std::ostream& operator<<(std::ostream& os, const Sets& s) {
  os << "Sets{"
     << "A=" << s.A << ", "
     << "B=" << s.B << ", "
     << "C=" << s.C << ", "
     << "A'=" << s.Ap << ", "
     << "B'=" << s.Bp << ", "
     << "C'=" << s.Cp << ", "
     << "A1=" << s.A1 << ", "
     << "B1=" << s.B1 << ", "
     << "C1=" << s.C1 << ", "
     << "A2=" << s.A2 << ", "
     << "B2=" << s.B2 << ", "
     << "C2=" << s.C2 << "}";
  return os;
}

std::vector<std::vector<Int>> product(const Set& s, Int repeat) {
  std::vector<std::vector<Int>> ans;
  std::stack<std::vector<Int>> call_stack;
  call_stack.push({});

  while (!call_stack.empty()) {
    const auto top = call_stack.top();
    call_stack.pop();

    if (top.size() >= repeat) {
      ans.push_back(top);
      continue;
    }

    for (const auto& e : s) {
      std::vector<Int> vt = top;
      vt.push_back(e);
      call_stack.push(vt);
    }
  }

  return ans;
}

std::vector<std::vector<Set>> disjoints(const Int n, const Int nsets) {
  std::vector<std::vector<Set>> ds{};
  const auto N = iter::range(nsets+1);
  for (const auto& bits : product(Set{N.begin(), N.end()}, n)) {
    std::deque<Set> S(nsets+1);
    for (size_t i = 0; i < bits.size(); ++i) {
      S[bits[i]].insert(i + 1);
    }
    S.pop_front();
    ds.emplace_back(S.begin(), S.end());
  }

  return ds;
}

std::vector<Int> tau(const Set& X) {
  const size_t d = X.size();

  std::vector<Int> ss{X.begin(), X.end()};
  std::sort(ss.begin(), ss.end());

  const auto R = iter::range(1, static_cast<int32_t>(d)+1);
  const std::vector<Int> rr{R.begin(), R.end()};

  std::vector<Int> tt;
  for (size_t idx = 0; idx < d; ++idx) {
    tt.push_back(ss[idx]-rr[idx]);
  }
  std::reverse(tt.begin(), tt.end());

  return tt;
}

template<typename T>
std::vector<Set> collect(const T& t) {
  std::vector<Set> ans{};
  for (const auto& e : t) {
    ans.push_back({e.begin(), e.end()});
  }
  return ans;
}

// Returns the empty set in the case that k == 0.
std::vector<Set> combinations(const std::vector<Int> R, const size_t k) {
  if (k == 0) return {{}};
  return collect(iter::combinations(R, k));
}

bool is_good(Int n, const std::vector<Int>& R, const Set& S,
             bool (*cond)(int64_t), const Set& A, const Set& B, const Set& C,
             const Set& Ap, const Set& Bp, const Set& Cp, Sets* s) {
  for (const auto& a : iter::product<2>(combinations(R, Ap.size()))) {
    const Set& A1 = std::get<0>(a);
    const Set& A2 = std::get<1>(a);
    if (!cond(nlnum::lrcoef(tau(Ap), tau(A1), tau(A2)))) continue;

    for (const auto& b : iter::product<2>(combinations(R, Bp.size()))) {
      const Set& B1 = std::get<0>(b);
      const Set& B2 = std::get<1>(b);
      if (!cond(nlnum::lrcoef(tau(Bp), tau(B1), tau(B2)))) continue;

      for (const auto& c : iter::product<2>(combinations(R, Cp.size()))) {
        const Set& C1 = std::get<0>(c);
        const Set& C2 = std::get<1>(c);
        if (!cond(nlnum::lrcoef(tau(Cp), tau(C1), tau(C2)))) continue;

        if (!cond(nlnum::lrcoef(tau(C), tau(A1), tau(B2))) ||
            !cond(nlnum::lrcoef(tau(B), tau(C1), tau(A2))) ||
            !cond(nlnum::lrcoef(tau(A), tau(B1), tau(C2))))
          continue;

        if (s != nullptr) {
          s->A = A;
          s->B = B;
          s->C = C;
          s->Ap = Ap;
          s->Bp = Bp;
          s->Cp = Cp;
          s->A1 = A1;
          s->B1 = B1;
          s->C1 = C1;
          s->A2 = A2;
          s->B2 = B2;
          s->C2 = C2;
        }
        return true;
      }
    }
  }
  return false;
}

std::vector<Sets> grand_ineqs(const Int n, bool (*cond)(int64_t)) {
  const auto djs = disjoints(n, 2);
  const auto RR = iter::range(1, static_cast<int32_t>(n)+1);
  const std::vector<int32_t> Rv{RR.begin(), RR.end()};
  const std::vector<Int> R{Rv.begin(), Rv.end()};
  const Set S{R.begin(), R.end()};

  std::vector<Sets> ans{};
#pragma omp parallel for schedule(dynamic)
  for (auto ita = djs.begin(); ita < djs.end(); ++ita) {
    const auto& Aq = *ita;
    const Set& A = Aq[0];
    const Set& Ap = Aq[1];
    for (const auto& Bq : djs) {
      const Set& B = Bq[0];
      const Set& Bp = Bq[1];
      for (const auto& Cq : djs) {
        const Set& C = Cq[0];
        const Set& Cp = Cq[1];

        if (A.size() != Bp.size() + Cp.size()
         || B.size() != Ap.size() + Cp.size()
         || C.size() != Ap.size() + Bp.size())
          continue;

        Sets s;
        if (is_good(n, R, S, cond, A, B, C, Ap, Bp, Cp, &s)) {
#pragma omp critical
          ans.push_back(s);
        }
      }
    }
  }

  return ans;
}

bool grand(const std::vector<Sets>& gi, const nlnum::Partition& lam,
           const nlnum::Partition& mu, const nlnum::Partition& nu) {
  const auto sum = [](const nlnum::Partition& pi, const Set& X) -> int64_t {
    Int ans = 0;
    for (Int ii : X) {
      if (ii > pi.size()) continue;
      ans += pi[ii - 1];
    }
    return static_cast<int64_t>(ans);
  };

  for (const auto& s : gi) {
    if (0 > (sum(mu, s.A) - sum(mu, s.Ap) + sum(nu, s.B) - sum(nu, s.Bp) +
             sum(lam, s.C) - sum(lam, s.Cp))) {
      return false;
    }
  }

  return true;
}

std::vector<std::vector<nlnum::Partition>> flagger(
    const Int n, const Int k, const std::vector<Sets>& gi) {
  const auto get_ks =
      [](const Int k) -> std::vector<std::tuple<Int, Int, Int>> {
    std::vector<std::tuple<Int, Int, Int>> ks;
    for (Int km = 0; km <= k; ++km) {
      for (Int kn = 0; kn <= km; ++kn) {
        const auto diff = km > kn ? km - kn : kn - km;
        const auto hi = std::min(kn, kn + km);
        for (Int kl = diff; kl <= hi; ++kl) {
          if ((km + kn + kl) % 2 == 1) continue;
          ks.push_back({km, kn, kl});
        }
      }
    }
    return ks;
  };

  std::vector<std::vector<nlnum::Partition>> ans;
  const auto ks = get_ks(k);

#pragma omp parallel for schedule(dynamic)
  for (Int idx = 0; idx < ks.size(); ++idx) {
    const auto km = std::get<0>(ks[idx]);
    const auto kn = std::get<1>(ks[idx]);
    const auto kl = std::get<2>(ks[idx]);

    const nlnum::PartitionsIn pim =
        nlnum::PartitionsIn(nlnum::Partition(n, km), km);
    const nlnum::PartitionsIn pin =
        nlnum::PartitionsIn(nlnum::Partition(n, kn), kn);
    const nlnum::PartitionsIn pil =
        nlnum::PartitionsIn(nlnum::Partition(n, kl), kl);

    const auto get =
        [](const nlnum::PartitionsIn pi) -> std::vector<nlnum::Partition> {
      std::vector<nlnum::Partition> acc;
      for (const auto& e : pi) {
        acc.push_back(e);
      }
      return acc;
    };

    const auto vim = get(pim);
    const auto vin = get(pin);
    const auto vil = get(pil);

    const auto pp = iter::product(vim, vin, vil);

    for (const auto& P : pp) {
      const auto& mu = std::get<0>(P);
      const auto& nu = std::get<1>(P);
      const auto& lam = std::get<2>(P);

      // If mu < nu or nu < lam, skip these.
      if (std::lexicographical_compare(mu.begin(), mu.end(), nu.begin(),
                                       nu.end()))
        continue;
      if (std::lexicographical_compare(nu.begin(), nu.end(), lam.begin(),
                                       lam.end()))
        continue;

      const auto pos = static_cast<bool>(nlnum::nlcoef(mu, nu, lam, true));
      bool satisfies = true;
      for (const auto& sx :
           iter::permutations(std::vector<nlnum::Partition>{mu, nu, lam})) {
        const nlnum::Partition& muu = sx[0];
        const nlnum::Partition& nuu = sx[1];
        const nlnum::Partition& lamm = sx[2];
        bool sat = false;
#pragma omp critical
        sat = grand(gi, lamm, muu, nuu);
        satisfies &= sat;
        if (pos && !sat) {
#pragma omp critical
          ans.push_back({muu, nuu, lamm});
          break;
        }
      }
      if (!pos && satisfies) {
#pragma omp critical
        ans.push_back({mu, nu, lam});
      }
    }
  }

  return ans;
}

}  // namespace saturation
