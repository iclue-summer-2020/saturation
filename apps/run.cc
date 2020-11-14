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

using saturation::Set;

int main(int argc, char** argv) {
  gflags::SetUsageMessage("Finds counterexamples to the NL numbers claim");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  return EXIT_SUCCESS;
}
