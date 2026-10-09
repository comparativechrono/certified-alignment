#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace cba {

struct Prediction {
  std::vector<double> centers;
  std::uint64_t probe_cells = 0;
};

struct Options {
  std::string mode = "linear";
  std::uint64_t width = 8;
  int block = 64;
  int seed_k = 11;
  // A rational keeps the optional approximation certificate exact.
  std::uint64_t epsilon_num = 0;
  std::uint64_t epsilon_den = 1;
  std::uint64_t max_cells = 50000000;
};

struct RoundTrace {
  std::uint64_t width = 0;
  std::uint64_t cells = 0;
  int upper = -1;  // -1 denotes infinity/unavailable in external results.
  int lower = -1;
};

struct Result {
  int distance = -1;
  int lower_bound = -1;
  bool certified_exact = false;
  bool certified_approx = false;
  std::string cigar;
  std::string error;
  std::uint64_t rounds = 0;
  std::uint64_t dp_cells = 0;
  std::uint64_t certificate_edges = 0;
  std::uint64_t peak_cells = 0;
  double dp_ms = 0;
  double certificate_ms = 0;
  std::vector<RoundTrace> trace;
};

Result align(const std::string& x, const std::string& y,
             const Prediction& prediction, const Options& options);

}  // namespace cba
