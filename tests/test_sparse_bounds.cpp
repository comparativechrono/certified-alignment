// Independent sparse threshold audit: the oracle minimizes over EVERY target
// chunk endpoint. It does not use the engine's three-candidate recurrence,
// packed matching, interval cache, or partial-block query shortcut.
#include "sparse_bounds.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
std::uint64_t pairs = 0, indices = 0, threshold_checks = 0, suffix_checks = 0;
std::uint64_t limited_checks = 0, interval_checks = 0, cap_checks = 0;
std::uint64_t partial_positions = 0, boundary_positions = 0, byte_pairs = 0;
std::uint64_t complete_proofs = 0, limited_proofs = 0, guard_rejections = 0;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

std::vector<std::vector<int>> edit_suffixes(const std::string& x,
                                           const std::string& y) {
  std::vector<std::vector<int>> d(x.size() + 1,
                                std::vector<int>(y.size() + 1));
  for (std::size_t i = x.size() + 1; i-- > 0;)
    for (std::size_t j = y.size() + 1; j-- > 0;) {
      if (i == x.size()) d[i][j] = static_cast<int>(y.size() - j);
      else if (j == y.size()) d[i][j] = static_cast<int>(x.size() - i);
      else d[i][j] = std::min({d[i + 1][j] + 1, d[i][j + 1] + 1,
                             d[i + 1][j + 1] + (x[i] != y[j])});
    }
  return d;
}

struct SlowBlocks {
  std::vector<std::size_t> boundaries;
  std::vector<std::vector<int>> table;
  std::size_t m;

  SlowBlocks(const std::string& x, const std::string& y, int k) : m(y.size()) {
    for (std::size_t s = 0; s < x.size(); s += k) boundaries.push_back(s);
    boundaries.push_back(x.size());
    table.assign(boundaries.size(), std::vector<int>(m + 1));
    for (std::size_t j = 0; j <= m; ++j) table.back()[j] = static_cast<int>(m - j);
    for (std::size_t b = boundaries.size() - 1; b-- > 0;) {
      const auto start = boundaries[b], count = boundaries[b + 1] - start;
      const std::string source = x.substr(start, count);
      for (std::size_t j = 0; j <= m; ++j) {
        int best = std::numeric_limits<int>::max();
        for (std::size_t end = j; end <= m; ++end) {
          const bool match = source == y.substr(j, end - j);
          const int difference = std::abs(static_cast<int>(count) -
                                          static_cast<int>(end - j));
          const int cost = match ? 0 : std::max(1, difference);
          best = std::min(best, cost + table[b + 1][end]);
        }
        table[b][j] = best;
      }
    }
  }

  int query(std::size_t i, std::size_t j) const {
    const auto it = std::lower_bound(boundaries.begin(), boundaries.end(), i);
    const auto b = static_cast<std::size_t>(it - boundaries.begin());
    if (*it == i) return table[b][j];
    int best = std::numeric_limits<int>::max();
    for (std::size_t end = j; end <= m; ++end)
      best = std::min(best, std::abs(static_cast<int>(*it - i) -
                                   static_cast<int>(end - j)) + table[b][end]);
    return best;
  }
};

void audit_cache(const cba::SparseOrderedBounds& index, const SlowBlocks& slow) {
  index.visit_cached([&](std::size_t b, std::size_t j, int lo, int hi) {
    require(b < slow.table.size() && j <= slow.m, "cache state outside table");
    const auto expected = slow.table[b][j];
    require(0 <= lo && lo <= expected && expected <= hi,
            "cached interval excludes slow relaxed optimum");
    ++interval_checks;
  });
}

std::vector<std::pair<std::size_t, std::size_t>> positions_for(
    const std::string& x, const std::string& y, std::mt19937& rng, bool sampled) {
  std::vector<std::pair<std::size_t, std::size_t>> positions;
  if (!sampled) {
    for (std::size_t i = 0; i <= x.size(); ++i)
      for (std::size_t j = 0; j <= y.size(); ++j) positions.emplace_back(i, j);
  } else {
    for (auto i : std::vector<std::size_t>{0, 1, 7, 8, 9, 31, 63, 64, 65, x.size()})
      if (i <= x.size())
        for (std::size_t j = 0; j <= y.size(); ++j) positions.emplace_back(i, j);
    for (int r = 0; r < 180; ++r)
      positions.emplace_back(rng() % (x.size() + 1), rng() % (y.size() + 1));
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
  }
  std::shuffle(positions.begin(), positions.end(), rng);
  return positions;
}

void check_pair(const std::string& x, const std::string& y, std::mt19937& rng,
                bool sampled = false) {
  ++pairs;
  const auto truth = edit_suffixes(x, y);
  const auto positions = positions_for(x, y, rng, sampled);
  for (int k = 1; k <= 8; ++k) {
    const SlowBlocks slow(x, y, k);
    cba::SparseOrderedBounds index(x, y, k, std::numeric_limits<std::uint64_t>::max(),
                                  64ULL * 1024 * 1024);
    ++indices;
    // Query different thresholds in a shuffled order and interleave unrelated
    // suffixes. This exercises both cached lower and cached upper intervals.
    struct Query { std::size_t i, j; int threshold; };
    std::vector<Query> queries;
    for (const auto& p : positions) {
      const int f = slow.query(p.first, p.second);
      require(0 <= f && f <= truth[p.first][p.second], "slow relaxation exceeds real ED");
      ++suffix_checks;
      if (p.first % k && p.first != x.size()) ++partial_positions;
      if (p.second == 63 || p.second == 64 || p.second == 65) ++boundary_positions;
      for (int t : {0, 1, std::max(0, f - 1), f, f + 1, f + 7})
        queries.push_back({p.first, p.second, t});
    }
    std::shuffle(queries.begin(), queries.end(), rng);
    for (const auto& q : queries) {
      const int f = slow.query(q.i, q.j);
      const int observed = index.query(q.i, q.j, q.threshold);
      const auto context = "n=" + std::to_string(x.size()) + ",m=" +
          std::to_string(y.size()) + ",k=" + std::to_string(k) + ",i=" +
          std::to_string(q.i) + ",j=" + std::to_string(q.j) + ",T=" +
          std::to_string(q.threshold);
      require(0 <= observed && observed <= f, "invalid sparse lower bound: " + context);
      require((observed >= q.threshold) == (f >= q.threshold),
              "uncapped proof predicate differs from slow oracle: " + context);
      ++threshold_checks;
      complete_proofs += observed >= q.threshold;
    }
    audit_cache(index, slow);
    require(index.stats().bytes <= 64ULL * 1024 * 1024, "uncapped index exceeded memory cap");
    guard_rejections += index.stats().guard_rejections;
  }
}

void check_limits(std::mt19937& rng) {
  std::string x, y;
  for (int i = 0; i < 129; ++i) x += "ACGT"[rng() % 4];
  y = x.substr(0, 19) + "TTTTTTGCGCGC" + x.substr(19, 44) + x.substr(75);
  const auto positions = positions_for(x, y, rng, true);
  for (int k : {1, 3, 8}) {
    const SlowBlocks slow(x, y, k);
    const auto layout = cba::SparseOrderedBounds::checked_size(x.size(), y.size(), k);
    for (std::uint64_t work : {0ULL, 1ULL, 2ULL, 7ULL, 64ULL, 1000ULL}) {
      for (std::uint64_t memory : {0ULL, 1ULL, 256ULL, 4096ULL, 1048576ULL}) {
        cba::SparseOrderedBounds index(x, y, k, work, memory);
        for (int repeat = 0; repeat < 2; ++repeat)
          for (const auto& p : positions) {
            const int expected = slow.query(p.first, p.second);
            const int t = repeat ? expected : std::max(1, expected - 1);
            const int v = index.query(p.first, p.second, t);
            require(0 <= v && v <= expected, "limited query returned invalid lower bound");
            limited_proofs += v >= t;
            ++limited_checks;
          }
        // Especially important after UNKNOWN propagation: incomplete proof is
        // not evidence for either a raised lower or a reduced upper endpoint.
        audit_cache(index, slow);
        require(index.stats().expansions <= work, "work cap exceeded");
        require(index.stats().bytes <= memory, "memory cap exceeded");
        cap_checks += 2;
      }
    }
    require(layout.initial_bytes > 0, "empty initial sparse payload");
    for (auto memory : {layout.initial_bytes - 1, layout.initial_bytes}) {
      cba::SparseOrderedBounds index(x, y, k, 1000000, memory);
      for (const auto& p : positions) {
        const auto expected = slow.query(p.first, p.second);
        require(index.query(p.first, p.second, expected) <= expected,
                "directory/tile allocation boundary changed soundness");
        ++limited_checks;
      }
      audit_cache(index, slow);
      require(index.stats().bytes <= memory, "initial allocation boundary exceeded cap");
      ++cap_checks;
    }
  }
  for (int k : {-1, 0, 9}) {
    bool caught = false;
    try { (void)cba::SparseOrderedBounds::checked_size(4, 4, k); }
    catch (const std::exception&) { caught = true; }
    require(caught, "invalid block size accepted");
    ++cap_checks;
  }
  bool overflow_caught = false;
  try {
    (void)cba::SparseOrderedBounds::checked_size(std::numeric_limits<std::size_t>::max(),
                                                std::numeric_limits<std::size_t>::max(), 1);
  } catch (const std::exception&) { overflow_caught = true; }
  require(overflow_caught, "overflow layout accepted");
  ++cap_checks;
  cba::SparseOrderedBounds index(x, y, 8, 100, 1048576);
  for (auto p : std::vector<std::pair<std::size_t, std::size_t>>{{x.size() + 1, 0}, {0, y.size() + 1}}) {
    bool caught = false;
    try { (void)index.query(p.first, p.second, 1); }
    catch (const std::out_of_range&) { caught = true; }
    require(caught, "out-of-range suffix accepted");
    ++cap_checks;
  }
}
}  // namespace

int main() {
  try {
    std::mt19937 rng(872491);
    std::vector<std::string> words;
    for (int n = 0; n <= 4; ++n)
      for (int bits = 0; bits < (1 << n); ++bits) {
        std::string word;
        for (int i = 0; i < n; ++i) word += (bits >> i) & 1 ? 'C' : 'A';
        words.push_back(word);
      }
    for (const auto& x : words) for (const auto& y : words) check_pair(x, y, rng);
    for (int sample = 0; sample < 36; ++sample) {
      std::string alphabet = sample < 18 ? "ACGT" : "aN?7-.*9#[]ACGTacgtNRYSwkmBDHV";
      if (sample >= 27) {
        alphabet.clear();
        for (int byte = 0; byte < 256; ++byte) alphabet += static_cast<char>(byte);
        ++byte_pairs;
      }
      std::string x, y;
      for (int n = rng() % 42; n > 0; --n) x += alphabet[rng() % alphabet.size()];
      for (int n = rng() % 42; n > 0; --n) y += alphabet[rng() % alphabet.size()];
      check_pair(x, y, rng);
    }
    for (const auto& p : std::vector<std::pair<std::string, std::string>>{
             {"A", "CC"}, {std::string(65, 'A'), std::string(64, 'A')},
             {std::string(64, 'A'), std::string(65, 'C')},
             {"AAAACCCCGGGGTTTT", "TTTTGGGGCCCCAAAA"},
             {std::string("\0\1\2\3\377\0\200", 7), std::string("\200\0\3\377\1", 5)}})
      check_pair(p.first, p.second, rng, true);
    std::string bytes;
    for (int b = 0; b < 256; ++b) bytes += static_cast<char>(b);
    check_pair(bytes, bytes.substr(63) + bytes.substr(0, 63), rng, true);
    ++byte_pairs;
    check_limits(rng);
    require(guard_rejections > 0 && partial_positions > 0 && boundary_positions > 0,
            "required guard, partial or tile-boundary paths were not exercised");
    std::cout << "{\"passed\":true,\"pairs_checked\":" << pairs
              << ",\"indices_checked\":" << indices
              << ",\"suffix_oracle_checks\":" << suffix_checks
              << ",\"uncapped_threshold_checks\":" << threshold_checks
              << ",\"capped_threshold_checks\":" << limited_checks
              << ",\"cached_interval_checks\":" << interval_checks
              << ",\"cap_and_validation_checks\":" << cap_checks
              << ",\"partial_suffix_positions\":" << partial_positions
              << ",\"tile_boundary_positions\":" << boundary_positions
              << ",\"arbitrary_byte_pairs\":" << byte_pairs
              << ",\"uncapped_successful_proofs\":" << complete_proofs
              << ",\"capped_successful_proofs\":" << limited_proofs
              << ",\"guard_rejections\":" << guard_rejections
              << ",\"allocation_failure_injection_tested\":false}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "sparse-bound audit failed: " << error.what() << '\n';
    return 1;
  }
}
