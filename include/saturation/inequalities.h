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
nlnum::Partition tau(const Set& I);

// For I \subset [4n], complement(I) = {4n+1-i|i\in I}.
Set complement(const Set& I, const Int n);

}  // namespace saturation

#endif  // SATURATION_INEQUALITIES_H_
