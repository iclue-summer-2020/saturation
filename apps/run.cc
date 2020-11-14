// Copyright (c) 2020 ICLUE @ UIUC. All rights reserved.

#include <gflags/gflags.h>
#include <saturation/inequalities.h>
#include <prettyprint.hpp>

#include <iostream>
#include <map>
#include <tuple>
#include <vector>

DEFINE_uint32(n, 0, "n");
DEFINE_uint32(k, 0, "k");

using saturation::flagger;
using saturation::grand_ineqs;
using saturation::Set;
using saturation::Sets;

int main(int argc, char** argv) {
  gflags::SetUsageMessage("Finds counterexamples to the NL numbers claim");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  const auto n = FLAGS_n;
  const auto k = FLAGS_k;

  const auto sets = grand_ineqs(n, [](int64_t c) -> bool { return c > 0; });
  std::cout << "Number of sets: " << sets.size() << std::endl;

  const auto flagged = flagger(n, k, sets);

  for (const auto& flag : flagged) {
    std::cout << flag << std::endl;
  }

  std::cout << "===========================" << std::endl;

  return EXIT_SUCCESS;
}
