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
  std::string bound = "length";
  // Retained suffix-bound index payload; the DP budget is independent.
  std::uint64_t bound_memory = 64ULL * 1024 * 1024;
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
  double bound_preprocess_ms = 0;
  std::uint64_t bound_bytes = 0;
  std::uint64_t bound_calls = 0;
  std::uint64_t composition_calls = 0;
  std::uint64_t q2_calls = 0;
  std::uint64_t q4_calls = 0;
  std::uint64_t bound_cursor_updates = 0;
  std::string bound_active = "length";
  std::string bound_fallback;
  std::vector<RoundTrace> trace;
};

Result align(const std::string& x, const std::string& y,
             const Prediction& prediction, const Options& options);

}  // namespace cba
