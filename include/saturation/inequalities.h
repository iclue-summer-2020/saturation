// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#ifndef SATURATION_INEQUALITIES_H_
#define SATURATION_INEQUALITIES_H_

#include <cstdint>
#include <ostream>
#include <set>
#include <vector>

#include <nlnum/nlnum.h>

namespace saturation {

// `Int` is a type alias for non-negative integers.
using Int = uint64_t;
using Set = std::set<Int>;

struct Sets {
    Set I, J, K, bI, bJ, bK;
    Sets() {}
    Sets(
      const Set& I, const Set& J, const Set& K,
      const Set& bI, const Set& bJ, const Set& bK
    ) : I{I}, J{J}, K{K}, bI{bI}, bJ{bJ}, bK{bK} { }
};

std::ostream& operator<<(std::ostream&, const Sets&);

// Computes the tau function defined in the pdf paper.
nlnum::Partition Tau(const Set& I);

// For I \subset [4n], bar(I) = {4n+1-i|i\in I}.
Set Bar(const Set& I, const Int n);

// For a partition lam = (lam_1, ..., lam_k) inside of the partition (a^b),
// define check(lam) to be the partition (a-lam_b, a-lam_{b-1}, ..., a-lam_1),
// where lam_i = 0 for i > k.
nlnum::Partition Check(const nlnum::Partition& lam, const Int a, const Int b);

// Refer to the pdf.
// We must ensure
//   * X \subseteq Y \subseteq [4n],
//   * b \in {0,2}.
Set Chi(const Set& X, const Set& Y, const Int n, const Int b);

// Returns all pairs (X, Bar(X)) of disjoint subsets of [4n] where |X|=r.
std::vector<std::pair<Set, Set>> Disjoints(const Int n, const Int r);

// Returns [4n] - X.
Set Complement(const Set& X, const Int n);

// Determines whether or not the sets satisfy the saturation properties.
// This function assumes the precomputation of the complement sets.
bool IsGood(
    const Int n, const Int r,
    const Set& I, const Set& J, const Set& K,
    const Set& bI, const Set& bJ, const Set& bK,
    Sets* s);

// Computes all of the sets in the saturation inequalities.
std::vector<Sets> SatIneqs(const Int n, const Int r);

}  // namespace saturation

#endif  // SATURATION_INEQUALITIES_H_
