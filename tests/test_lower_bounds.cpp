// Independent audit of suffix lower bounds. This reference uses string-keyed
// q-gram maps and unrestricted suffix DP; it does not share cursor or hashing
// code with the implementation under test.
#include "lower_bounds.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
std::uint64_t pairs_checked = 0, queries_checked = 0, threshold_checks = 0;
std::uint64_t strict_composition_improvements = 0, strict_qgram_improvements = 0;
std::uint64_t backwards_steps = 0, resets = 0;

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

int gram_bound(const std::string& x, const std::string& y,
               std::size_t i, std::size_t j, int q) {
  std::map<std::string, int> difference;
  for (std::size_t at = i; at + q <= x.size(); ++at) ++difference[x.substr(at, q)];
  for (std::size_t at = j; at + q <= y.size(); ++at) --difference[y.substr(at, q)];
  int surplus = 0, deficit = 0;
  for (const auto& item : difference) {
    surplus += std::max(0, item.second);
    deficit += std::max(0, -item.second);
  }
  return (std::max(surplus, deficit) + q - 1) / q;
}

std::vector<std::vector<int>> suffix_distances(const std::string& x,
                                             const std::string& y) {
  std::vector<std::vector<int>> distance(x.size() + 1,
                                       std::vector<int>(y.size() + 1));
  for (std::size_t i = x.size() + 1; i-- > 0;) {
    for (std::size_t j = y.size() + 1; j-- > 0;) {
      if (i == x.size()) distance[i][j] = static_cast<int>(y.size() - j);
      else if (j == y.size()) distance[i][j] = static_cast<int>(x.size() - i);
      else distance[i][j] = std::min({distance[i + 1][j] + 1,
                                    distance[i][j + 1] + 1,
                                    distance[i + 1][j + 1] + (x[i] != y[j])});
    }
  }
  return distance;
}

bool small_alphabet(const std::string& x, const std::string& y) {
  std::array<bool, 256> used{};
  for (unsigned char symbol : x) used[symbol] = true;
  for (unsigned char symbol : y) used[symbol] = true;
  return std::count(used.begin(), used.end(), true) <= 4;
}

void check_pair(const std::string& x, const std::string& y, std::mt19937& rng,
                bool all_suffixes) {
  ++pairs_checked;
  const auto distance = suffix_distances(x, y);
  std::vector<std::pair<std::size_t, std::size_t>> positions;
  if (all_suffixes) {
    for (std::size_t i = 0; i <= x.size(); ++i)
      for (std::size_t j = 0; j <= y.size(); ++j) positions.emplace_back(i, j);
  } else {
    positions = {{0, 0}, {x.size(), y.size()}, {0, y.size()}, {x.size(), 0}};
    for (int sample = 0; sample < 192; ++sample)
      positions.emplace_back(rng() % (x.size() + 1), rng() % (y.size() + 1));
  }
  // Deliberately nonmonotone: callers may revisit shorter and longer suffixes.
  std::shuffle(positions.begin(), positions.end(), rng);
  positions.emplace_back(x.size(), y.size());
  positions.emplace_back(0, 0);
  const bool exact_grams = small_alphabet(x, y);
  for (const std::string kind : {"length", "composition", "qgram"}) {
    cba::Options options;
    options.bound = kind;
    cba::Result metrics;
    cba::SuffixBounds bounds(x, y, options, metrics);
    std::array<std::pair<std::size_t, std::size_t>, 5> previous{};
    for (int pass = 0; pass < 2; ++pass) {
      if (pass) {
        bounds.reset_round();
        ++resets;
        std::reverse(positions.begin(), positions.end());
      }
      for (std::size_t index = 0; index < positions.size(); ++index) {
        const auto [i, j] = positions[index];
        const auto stream_index = (index + pass * 2) % 5;
        const auto stream = static_cast<cba::ExitStream>(stream_index);
        if (i < previous[stream_index].first || j < previous[stream_index].second)
          ++backwards_steps;
        previous[stream_index] = {i, j};
        const int length = std::abs(static_cast<int>(x.size() - i) -
                                    static_cast<int>(y.size() - j));
        const int composition = std::max(length, gram_bound(x, y, i, j, 1));
        const int exact_qgrams = std::max({composition, gram_bound(x, y, i, j, 2),
                                           gram_bound(x, y, i, j, 4)});
        const int observed = bounds.query(i, j, static_cast<int>(x.size() + y.size() + 1), stream);
        ++queries_checked;
        const std::string context = "pair=" + std::to_string(pairs_checked) +
            " kind=" + kind + " suffix=(" + std::to_string(i) + "," + std::to_string(j) + ")";
        require(observed >= length, "weaker than length: " + context);
        require(observed <= distance[i][j], "inadmissible: " + context);
        if (kind == "length") require(observed == length, "length differs: " + context);
        if (kind == "composition")
          require(observed == composition, "composition differs: " + context);
        if (kind == "qgram") {
          require(observed >= composition, "qgram weaker than composition: " + context);
          require(observed <= exact_qgrams, "bucketed qgram exceeds exact histogram: " + context);
          if (exact_grams) require(observed == exact_qgrams, "exact qgram differs: " + context);
        }
        strict_composition_improvements += kind == "composition" && observed > length;
        strict_qgram_improvements += kind == "qgram" && observed > composition;
        // Include a threshold at zero and just above/below the true bound.
        // The implementation may short-circuit once it has reached threshold;
        // it must retain at least min(full bound, threshold) and remain safe.
        for (int threshold : {0, 1, std::max(0, observed - 1), observed, observed + 1}) {
          const int limited = bounds.query(i, j, threshold, stream);
          ++threshold_checks;
          require(limited >= std::min(observed, threshold), "threshold weakens too far: " + context);
          require(limited >= length, "threshold weakens length: " + context);
          require(limited <= observed, "threshold exceeds unlimited bound: " + context);
          require(limited <= distance[i][j], "threshold inadmissible: " + context);
        }
      }
    }
  }
}

void check_memory_fallback() {
  const std::string x = "ACGTTGCAACGT", y = "GTCAAGCTTGCA";
  for (std::uint64_t cap : {0ULL, 1ULL, 63ULL, 1024ULL}) {
    cba::Options options;
    options.bound = "qgram";
    options.bound_memory = cap;
    cba::Result metrics;
    cba::SuffixBounds bounds(x, y, options, metrics);
    for (std::size_t i = 0; i <= x.size(); ++i) {
      for (std::size_t j = 0; j <= y.size(); ++j) {
        const int observed = bounds.query(i, j, 1000, cba::ExitStream::Horizontal);
        const int length = std::abs(static_cast<int>(x.size() - i) - static_cast<int>(y.size() - j));
        require(observed == length, "small memory cap must fall back to length");
      }
    }
    require(metrics.bound_bytes <= cap, "bound memory accounting exceeds cap");
  }
  // Obtain the implementation's payload size without assuming an ABI. The
  // boundary cases permit one stage, omit later stages, or sit one byte below
  // an allocation; every permitted fallback must still remain admissible.
  cba::Options composition_options;
  composition_options.bound = "composition";
  cba::Result composition_metrics;
  cba::SuffixBounds composition(x, y, composition_options, composition_metrics);
  composition.query(0, 0, 1000, cba::ExitStream::Horizontal);
  cba::Options all_options;
  all_options.bound = "qgram";
  cba::Result all_metrics;
  cba::SuffixBounds all(x, y, all_options, all_metrics);
  all.query(0, 0, 1000, cba::ExitStream::Horizontal);
  const auto distance = suffix_distances(x, y);
  for (std::uint64_t cap : {composition_metrics.bound_bytes - 1,
                            composition_metrics.bound_bytes,
                            all_metrics.bound_bytes - 1, all_metrics.bound_bytes}) {
    cba::Options options;
    options.bound = "qgram";
    options.bound_memory = cap;
    cba::Result metrics;
    cba::SuffixBounds bounds(x, y, options, metrics);
    for (std::size_t i = 0; i <= x.size(); ++i)
      for (std::size_t j = 0; j <= y.size(); ++j) {
        const int observed = bounds.query(i, j, 1000, cba::ExitStream::Horizontal);
        const int length = std::abs(static_cast<int>(x.size() - i) - static_cast<int>(y.size() - j));
        require(length <= observed && observed <= distance[i][j], "partial fallback inadmissible");
        if (cap == composition_metrics.bound_bytes)
          require(observed == std::max(length, gram_bound(x, y, i, j, 1)),
                  "one-stage cap must retain the composition bound");
      }
    require(metrics.bound_bytes <= cap, "partial index memory exceeds cap");
  }
}
}  // namespace

int main() {
  try {
    std::mt19937 rng(430972);
    std::vector<std::string> words;
    for (int n = 0; n <= 4; ++n)
      for (int bits = 0; bits < (1 << n); ++bits) {
        std::string word;
        for (int i = 0; i < n; ++i) word.push_back((bits >> i) & 1 ? 'C' : 'A');
        words.push_back(word);
      }
    for (const auto& x : words) for (const auto& y : words) check_pair(x, y, rng, true);
    for (int sample = 0; sample < 80; ++sample) {
      std::string alphabet = sample < 40 ? "ACGT" : "aN?7-.*9#[]ACGTacgtNRYSwkmBDHV";
      if (sample >= 60) {
        alphabet.clear();
        for (int byte = 0; byte < 256; ++byte) alphabet.push_back(static_cast<char>(byte));
      }
      std::string x, y;
      for (int i = static_cast<int>(rng() % 75); i > 0; --i) x += alphabet[rng() % alphabet.size()];
      for (int j = static_cast<int>(rng() % 75); j > 0; --j) y += alphabet[rng() % alphabet.size()];
      check_pair(x, y, rng, false);
    }
    check_pair("AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", rng, true);
    check_pair("AAAACCCCGGGGTTTT", "TTTTGGGGCCCCAAAA", rng, true);
    // More than 256 distinct bigrams guarantees collisions in a 256-bucket
    // implementation. Exact string-key reference still bounds it from above.
    std::string collision_x, collision_y;
    for (int i = 32; i < 50; ++i) for (int j = 32; j < 50; ++j) {
      collision_x += static_cast<char>(i);
      collision_x += static_cast<char>(j);
      collision_y += static_cast<char>(j);
      collision_y += static_cast<char>(i);
    }
    check_pair(collision_x, collision_y, rng, false);
    check_memory_fallback();
    require(backwards_steps > 0 && resets > 0, "cursor exercise was empty");
    require(strict_composition_improvements > 0, "composition never improved length");
    require(strict_qgram_improvements > 0, "qgram never improved composition");
    std::cout << "{\"passed\":true,\"pairs_checked\":" << pairs_checked
              << ",\"queries_checked\":" << queries_checked
              << ",\"threshold_checks\":" << threshold_checks
              << ",\"backwards_steps\":" << backwards_steps
              << ",\"round_resets\":" << resets
              << ",\"strict_composition_improvements\":" << strict_composition_improvements
              << ",\"strict_qgram_improvements\":" << strict_qgram_improvements << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "lower-bound audit failed: " << error.what() << '\n';
    return 1;
  }
}
