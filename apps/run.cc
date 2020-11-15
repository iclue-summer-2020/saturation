// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#include <gflags/gflags.h>
#include <saturation/inequalities.h>

#include <iostream>
#include <map>
#include <prettyprint.hpp>
#include <tuple>
#include <vector>

DEFINE_uint32(n, 0, "n");
DEFINE_uint32(r, 0, "r");

using saturation::Flagger;
using saturation::SatIneqs;
using saturation::Set;
using saturation::Sets;

int main(int argc, char** argv) {
  gflags::SetUsageMessage("Finds counterexamples to the NL numbers claim");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  const auto n = FLAGS_n;
  const auto r = FLAGS_r;
  const auto ineqs = SatIneqs(n, r);

  std::cout << "Number of sets: " << ineqs.size() << std::endl;

  for (const auto& sets : ineqs) {
    std::cout << sets << std::endl;
  }

  const auto flagged = Flagger(n, r);
  for (const auto& counterexample : flagged) {
    std::cout << "Flagged: " << counterexample << std::endl;
  }

  std::cout << "==============================" << std::endl;
  return EXIT_SUCCESS;
}
