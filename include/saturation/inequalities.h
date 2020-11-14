// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#ifndef SATURATION_INEQUALITIES_H_
#define SATURATION_INEQUALITIES_H_

#include <cstdint>
#include <ostream>
#include <set>
#include <vector>

#include <nlnum/nlnum.h>

namespace saturation {

typedef uint64_t Int;
typedef std::set<Int> Set;

struct Sets {
  Set A, B, C, Ap, Bp, Cp, A1, B1, C1, A2, B2, C2;
  Sets() {}
  Sets(
      const Set& A,
      const Set& B,
      const Set& C,
      const Set& Ap,
      const Set& Bp,
      const Set& Cp,
      const Set& A1,
      const Set& B1,
      const Set& C1,
      const Set& A2,
      const Set& B2,
      const Set& C2
  ) : A{A}, B{B}, C{C}, Ap{Ap}, Bp{Bp}, Cp{Cp},
    A1{A1}, B1{B1}, C1{C1}, A2{A2}, B2{B2}, C2{C2} { }
};

std::ostream& operator<<(std::ostream&, const Sets&);

// Computes the generalized Cartesian product of a set with itself.
std::vector<std::vector<Int>> product(const Set& s, Int repeat);

std::vector<Sets> grand_ineqs(Int n, bool (*cond)(int64_t));

// Computes the tau function defined in the NL paper.
std::vector<Int> tau(const Set& X);

// Returns all `nsets`-tuple of sets of [n] that are disjoint.
std::vector<std::vector<Set>> disjoints(Int n, Int nsets);

bool grand(const std::vector<Sets>& gi, const nlnum::Partition& lam,
           const nlnum::Partition& mu, const nlnum::Partition& nu);

// Flags counterexamples.
std::vector<std::vector<nlnum::Partition>> flagger(
    const Int n, const Int k, const std::vector<Sets>& gi);

}  // namespace saturation

#endif  // SATURATION_INEQUALITIES_H_
