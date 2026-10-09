// Independent ordered-bound audit. The reference enumerates EVERY target
// chunk endpoint; it does not use the implementation's three-candidate
// recurrence, packed word matching, or shifted partial-block query formula.
// All checks remain enabled under -DNDEBUG.
#include "ordered_bounds.hpp"
#include "lower_bounds.hpp"

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
std::uint64_t pairs_checked = 0, indices_checked = 0, suffix_checks = 0;
std::uint64_t threshold_checks = 0, combined_checks = 0, lipschitz_checks = 0;
std::uint64_t fallback_checks = 0, strict_ordered_improvements = 0;

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

  SlowBlocks(const std::string& x, const std::string& y, int block) : m(y.size()) {
    for (std::size_t s = 0; s < x.size(); s += block) boundaries.push_back(s);
    boundaries.push_back(x.size());
    table.assign(boundaries.size(), std::vector<int>(m + 1));
    for (std::size_t j = 0; j <= m; ++j) table.back()[j] = static_cast<int>(m - j);
    for (std::size_t row = boundaries.size() - 1; row-- > 0;) {
      const auto start = boundaries[row], count = boundaries[row + 1] - start;
      const std::string source_block = x.substr(start, count);
      for (std::size_t j = 0; j <= m; ++j) {
        int best = std::numeric_limits<int>::max();
        for (std::size_t end = j; end <= m; ++end) {
          const bool match = source_block == y.substr(j, end - j);
          const int difference = std::abs(static_cast<int>(count) -
                                          static_cast<int>(end - j));
          const int cost = match ? 0 : std::max(1, difference);
          best = std::min(best, cost + table[row + 1][end]);
        }
        table[row][j] = best;
      }
    }
  }

  int query(std::size_t i, std::size_t j) const {
    const auto at = std::lower_bound(boundaries.begin(), boundaries.end(), i);
    const auto row = static_cast<std::size_t>(at - boundaries.begin());
    if (*at == i) return table[row][j];
    // The incomplete source block is relaxed to a length-only cost. Enumerate
    // every possible consumed target prefix rather than using Lipschitzness.
    int best = std::numeric_limits<int>::max();
    const auto remaining = *at - i;
    for (std::size_t end = j; end <= m; ++end)
      best = std::min(best, std::abs(static_cast<int>(remaining) -
                                   static_cast<int>(end - j)) + table[row][end]);
    return best;
  }
};

void check_pair(const std::string& x, const std::string& y, std::mt19937& rng) {
  ++pairs_checked;
  const auto true_ed = edit_suffixes(x, y);
  std::vector<std::pair<std::size_t, std::size_t>> positions;
  for (std::size_t i = 0; i <= x.size(); ++i)
    for (std::size_t j = 0; j <= y.size(); ++j) positions.emplace_back(i, j);
  std::shuffle(positions.begin(), positions.end(), rng);
  for (int k = 1; k <= 8; ++k) {
    cba::OrderedBounds index(x, y, k);
    const SlowBlocks slow(x, y, k);
    ++indices_checked;
    cba::Options options;
    options.bound = "ordered";
    options.order_block = k;
    cba::Result metrics;
    cba::SuffixBounds combined(x, y, options, metrics);
    cba::Options q_options;
    q_options.bound = "qgram";
    cba::Result q_metrics;
    cba::SuffixBounds qgram(x, y, q_options, q_metrics);
    for (const auto& position : positions) {
      const auto i = position.first, j = position.second;
      const int observed = index.query(i, j), expected = slow.query(i, j);
      const std::string context = "n=" + std::to_string(x.size()) + ", m=" +
          std::to_string(y.size()) + ", k=" + std::to_string(k) + ", i=" +
          std::to_string(i) + ", j=" + std::to_string(j);
      require(observed == expected, "slow block oracle mismatch: " + context);
      require(0 <= observed && observed <= true_ed[i][j],
              "ordered bound inadmissible: " + context);
      ++suffix_checks;
      if (j < y.size()) {
        require(std::abs(observed - index.query(i, j + 1)) <= 1,
                "ordered row is not 1-Lipschitz: " + context);
        ++lipschitz_checks;
      }
      const auto stream = static_cast<cba::ExitStream>((i + 3 * j) % 5);
      const int old = qgram.query(i, j, 100000, stream);
      const int value = combined.query(i, j, 100000, stream);
      require(value == std::max(old, observed), "combined bound mismatch: " + context);
      require(value <= true_ed[i][j], "combined bound inadmissible: " + context);
      strict_ordered_improvements += value > old;
      ++combined_checks;
      for (int threshold : {0, 1, std::max(0, value - 1), value, value + 1}) {
        const int capped = combined.query(i, j, threshold, stream);
        // A threshold query can return the uncapped length bound. It must
        // preserve the exact value when below threshold and the threshold
        // predicate otherwise, without ever exceeding the true distance.
        require(capped <= true_ed[i][j], "threshold returned inadmissible value");
        require(std::min(capped, threshold) == std::min(value, threshold),
                "threshold changed bound predicate");
        ++threshold_checks;
      }
    }
    combined.reset_round();
    qgram.reset_round();
    for (auto position : std::vector<std::pair<std::size_t, std::size_t>>{
             {x.size(), y.size()}, {0, 0}, {x.size(), 0}, {0, y.size()}}) {
      const auto stream = cba::ExitStream::Horizontal;
      require(combined.query(position.first, position.second, 100000, stream) ==
                  std::max(index.query(position.first, position.second),
                    qgram.query(position.first, position.second, 100000, stream)),
              "reset changed combined bound");
      ++combined_checks;
    }
  }
}

void check_fallbacks() {
  const std::string x = "AAAACCCCTTTTGGGGACGTACGTA";
  const std::string y = "TTTTAAAAGGGGCCCCACGTAC";
  const auto d = edit_suffixes(x, y);
  cba::Options q_options;
  q_options.bound = "qgram";
  cba::Result q_metrics;
  cba::SuffixBounds qgram(x, y, q_options, q_metrics);
  qgram.query(0, 0, 100000, cba::ExitStream::Horizontal);
  for (int k : {1, 3, 8}) {
    const auto layout = cba::OrderedBounds::checked_size(x.size(), y.size(), k);
    const auto required_memory = q_metrics.bound_bytes + layout.bytes;
    for (int scenario = 0; scenario < 5; ++scenario) {
      cba::Options options;
      options.bound = "ordered";
      options.order_block = k;
      if (scenario == 0) options.order_max_cells = 0;
      if (scenario == 1) options.order_max_cells = layout.cells - 1;
      if (scenario == 2) options.order_max_cells = layout.cells;
      if (scenario == 3) options.bound_memory = required_memory - 1;
      if (scenario == 4) options.bound_memory = required_memory;
      cba::Result metrics;
      cba::SuffixBounds bound(x, y, options, metrics);
      const bool should_build = scenario == 2 || scenario == 4;
      for (std::size_t i = 0; i <= x.size(); ++i)
        for (std::size_t j = 0; j <= y.size(); ++j) {
          const int value = bound.query(i, j, 100000, cba::ExitStream::Horizontal);
          const int old = qgram.query(i, j, 100000, cba::ExitStream::Horizontal);
          require(value >= old && value <= d[i][j], "fallback lost admissibility or qgram bound");
          if (!should_build) require(value == old, "omitted ordered bound changed qgram");
          ++fallback_checks;
        }
      require(metrics.bound_bytes <= options.bound_memory, "shared bound memory cap exceeded");
      require((metrics.ordered_cells != 0) == should_build, "incorrect ordered cell-cap boundary");
      require((metrics.bound_active.find("+ordered") != std::string::npos) == should_build,
              "incorrect active bound label");
      require(metrics.ordered_fallback.empty() == should_build, "missing or unexpected fallback reason");
    }
  }
  cba::Options none_options;
  none_options.bound = "ordered";
  none_options.bound_memory = 0;
  cba::Result none_metrics;
  cba::SuffixBounds none(x, y, none_options, none_metrics);
  for (std::size_t i = 0; i <= x.size(); ++i)
    for (std::size_t j = 0; j <= y.size(); ++j) {
      require(none.query(i, j, 100000, cba::ExitStream::Horizontal) ==
          std::abs(static_cast<int>(x.size() - i) - static_cast<int>(y.size() - j)),
          "zero shared memory cap must retain length-only bound");
      ++fallback_checks;
    }
  require(none_metrics.bound_bytes == 0, "zero cap retained index payload");
  for (int k : {-1, 0, 9}) {
    bool caught = false;
    try { (void)cba::OrderedBounds::checked_size(4, 4, k); }
    catch (const std::exception&) { caught = true; }
    require(caught, "invalid ordered block accepted");
  }
  bool overflow_caught = false;
  try {
    (void)cba::OrderedBounds::checked_size(std::numeric_limits<std::size_t>::max(),
                                         std::numeric_limits<std::size_t>::max(), 1);
  } catch (const std::exception&) { overflow_caught = true; }
  require(overflow_caught, "overflow layout was accepted");
}
}  // namespace

int main() {
  try {
    std::mt19937 rng(973214);
    std::vector<std::string> words;
    for (int n = 0; n <= 4; ++n)
      for (int bits = 0; bits < (1 << n); ++bits) {
        std::string word;
        for (int i = 0; i < n; ++i) word += (bits >> i) & 1 ? 'C' : 'A';
        words.push_back(word);
      }
    for (const auto& x : words) for (const auto& y : words) check_pair(x, y, rng);
    for (int sample = 0; sample < 60; ++sample) {
      std::string alphabet = sample < 30 ? "ACGT" : "aN?7-.*9#[]ACGTacgtNRYSwkmBDHV";
      if (sample >= 45) {
        alphabet.clear();
        for (int byte = 0; byte < 256; ++byte) alphabet += static_cast<char>(byte);
      }
      std::string x, y;
      for (int length = rng() % 49; length > 0; --length) x += alphabet[rng() % alphabet.size()];
      for (int length = rng() % 49; length > 0; --length) y += alphabet[rng() % alphabet.size()];
      check_pair(x, y, rng);
    }
    for (const auto& pair : std::vector<std::pair<std::string, std::string>>{
             {"A", "CC"}, {"ACGTACGTACGTACGT", "GTACGTACGTACGTAC"},
             {"AAAACCCCGGGGTTTT", "TTTTGGGGCCCCAAAA"},
             {"AAAAAAAAAAAAAA", "AAA"}, {"ACGTTGCACCATGG", "ACGTCATGG"},
             {"ACGTTGCACCATGG", "TTACGTTGCACC"},
             {std::string("\0\1\2\3\377\0\200", 7), std::string("\200\0\3\377\1", 5)}})
      check_pair(pair.first, pair.second, rng);
    // Even k=1 is a relaxation: an unequal length-2 target chunk costs one
    // here, while transforming A into CC costs two real edit operations.
    const std::string one = "A", two = "CC";
    cba::OrderedBounds k_one(one, two, 1);
    require(k_one.query(0, 0) == 1 && edit_suffixes(one, two)[0][0] == 2,
            "k=1 relaxation semantics changed");
    check_fallbacks();
    require(strict_ordered_improvements > 0, "ordered bound never strengthened qgram");
    std::cout << "{\"passed\":true,\"pairs_checked\":" << pairs_checked
              << ",\"indices_checked\":" << indices_checked
              << ",\"suffix_checks\":" << suffix_checks
              << ",\"combined_checks\":" << combined_checks
              << ",\"threshold_checks\":" << threshold_checks
              << ",\"lipschitz_checks\":" << lipschitz_checks
              << ",\"fallback_checks\":" << fallback_checks
              << ",\"strict_ordered_improvements\":" << strict_ordered_improvements
              << ",\"allocation_failure_injection_tested\":false}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ordered-bound audit failed: " << error.what() << '\n';
    return 1;
  }
}
