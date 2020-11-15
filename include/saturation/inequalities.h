// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#ifndef SATURATION_INEQUALITIES_H_
#define SATURATION_INEQUALITIES_H_

#include <cstdint>
#include <ostream>
#include <set>
#include <vector>

#include <nlnum/nlnum.h>

namespace saturation {

using Int = uint64_t;
using Set = std::set<Int>;

// Computes the tau function defined in the pdf paper.
nlnum::Partition Tau(const Set& I);

// For I \subset [4n], complement(I) = {4n+1-i|i\in I}.
Set Complement(const Set& I, const Int n);

// For a partition lam = (lam_1, ..., lam_k) inside of the partition (a^b),
// define check(lam) to be the partition (a-lam_b, a-lam_{b-1}, ..., a-lam_1),
// where lam_i = 0 for i > k.
nlnum::Partition Check(const nlnum::Partition& lam, const Int a, const Int b);

}  // namespace saturation

#endif  // SATURATION_INEQUALITIES_H_
