#pragma once

#include "engine.hpp"
#include "ordered_bounds.hpp"
#include "sparse_bounds.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace cba {

// Separating the five boundary streams avoids repeatedly moving a suffix
// cursor between the left and right sides of a band. Backward moves remain
// valid, so no monotonicity assumption about the predictor is required.
enum class ExitStream {
  Horizontal = 0, VerticalLeft = 1, VerticalRight = 2,
  DiagonalLeft = 3, DiagonalRight = 4
};

class SuffixBounds {
 public:
  SuffixBounds(const std::string& x, const std::string& y,
               const Options& options, Result& metrics)
      : x_(x), y_(y), options_(options), metrics_(metrics) {}

  // Returns the combined lower bound, possibly capped at threshold. A value
  // >= threshold suffices when updating a running minimum. Length is checked
  // first and can be returned uncapped. Stages are constructed and queried
  // only when every cheaper bound is below threshold.
  int query(std::size_t i, std::size_t j, int threshold, ExitStream stream) {
    const auto rx = x_.size() - i, ry = y_.size() - j;
    int value = static_cast<int>(rx > ry ? rx - ry : ry - rx);
    if (value >= threshold || options_.bound == "length") return value;
    ++metrics_.bound_calls;
    for (int index = 0; index <
         (options_.bound == "qgram" || options_.bound == "ordered" ||
          options_.bound == "sparse" ? 3 : 1); ++index) {
      if (!ensure_stage(index)) continue;
      auto& stage = *stages_[index];
      auto& cursor = stage.cursors[static_cast<std::size_t>(stream)];
      if (index == 0) ++metrics_.composition_calls;
      else if (index == 1) ++metrics_.q2_calls;
      else ++metrics_.q4_calls;
      move_cursor(stage, cursor, i, j);
      const int mass = std::max(cursor.positive, cursor.negative);
      value = std::max(value, (mass + stage.q - 1) / stage.q);
      if (value >= threshold) return threshold;
    }
    if (options_.bound == "ordered" && ensure_ordered()) {
      ++metrics_.ordered_calls;
      value = std::max(value, ordered_->query(i, j));
      if (value >= threshold) return threshold;
    }
    if (options_.bound == "sparse") {
      value = std::max(value, sparse_query(i, j, threshold));
      if (value >= threshold) return threshold;
    }
    return value;
  }

  void reset_round() {
    for (auto& stage : stages_) {
      if (!stage) continue;
      for (auto& cursor : stage->cursors) initialize_cursor(*stage, cursor);
    }
  }

 private:
  struct Cursor {
    std::size_t i = 0, j = 0;
    std::array<int, 256> delta{};
    int positive = 0, negative = 0;
  };

  struct Stage {
    int q = 1;
    int bins = 0;
    std::vector<std::uint8_t> xgrams, ygrams;
    std::array<int, 256> initial_delta{};
    int initial_positive = 0, initial_negative = 0;
    std::array<Cursor, 5> cursors;
  };

  const std::string& x_;
  const std::string& y_;
  const Options& options_;
  Result& metrics_;
  std::array<std::unique_ptr<Stage>, 3> stages_;
  std::array<bool, 3> attempted_{};
  std::array<unsigned int, 256> symbol_{};
  unsigned int alphabet_ = 0;
  bool alphabet_ready_ = false;
  std::unique_ptr<OrderedBounds> ordered_;
  bool ordered_attempted_ = false;
  std::unique_ptr<SparseOrderedBounds> sparse_;
  bool sparse_attempted_ = false;

  int sparse_query(std::size_t i, std::size_t j, int threshold) {
    if (!sparse_attempted_) {
      sparse_attempted_ = true;
      const auto start = std::chrono::steady_clock::now();
      try {
        sparse_ = std::make_unique<SparseOrderedBounds>(x_, y_, options_.order_block,
            options_.sparse_max_work,
            metrics_.bound_bytes >= options_.bound_memory ? 0 :
                options_.bound_memory - metrics_.bound_bytes);
      } catch (const std::bad_alloc&) {
        metrics_.sparse_fallback = "allocation failed; omitted sparse evaluator";
        if (!metrics_.bound_fallback.empty()) metrics_.bound_fallback += "; ";
        metrics_.bound_fallback += metrics_.sparse_fallback;
      }
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - start).count();
      metrics_.sparse_preprocess_ms += elapsed;
      metrics_.bound_preprocess_ms += elapsed;
    }
    if (!sparse_) return 0;
    const auto prior_ms = sparse_->stats().preprocess_ms;
    const int value = sparse_->query(i, j, threshold);
    const auto& stats = sparse_->stats();
    if (!metrics_.sparse_bytes && stats.bytes) metrics_.bound_active += "+sparse";
    metrics_.bound_bytes += stats.bytes - metrics_.sparse_bytes;
    metrics_.sparse_bytes = stats.bytes;
    metrics_.sparse_queries = stats.queries;
    metrics_.sparse_expansions = stats.expansions;
    metrics_.sparse_cache_hits = stats.cache_hits;
    metrics_.sparse_guard_rejections = stats.guard_rejections;
    metrics_.sparse_cached_states = stats.cached_states;
    metrics_.sparse_tiles = stats.tiles;
    metrics_.sparse_peak_stack = stats.peak_stack;
    metrics_.sparse_preprocess_ms += stats.preprocess_ms - prior_ms;
    metrics_.bound_preprocess_ms += stats.preprocess_ms - prior_ms;
    if (metrics_.sparse_fallback.empty() && !stats.fallback.empty()) {
      metrics_.sparse_fallback = stats.fallback;
      if (!metrics_.bound_fallback.empty()) metrics_.bound_fallback += "; ";
      metrics_.bound_fallback += stats.fallback;
    }
    return value;
  }

  bool ensure_ordered() {
    if (ordered_attempted_) return static_cast<bool>(ordered_);
    ordered_attempted_ = true;
    const auto start = std::chrono::steady_clock::now();
    const auto size = OrderedBounds::checked_size(x_.size(), y_.size(), options_.order_block);
    if (size.cells > options_.order_max_cells) {
      metrics_.ordered_fallback = "order-max-cells limit omitted ordered index";
    } else if (metrics_.bound_bytes > options_.bound_memory ||
               size.bytes > options_.bound_memory - metrics_.bound_bytes) {
      metrics_.ordered_fallback = "shared bound-memory limit omitted ordered index";
    } else {
      try {
        ordered_ = std::make_unique<OrderedBounds>(x_, y_, options_.order_block);
        metrics_.ordered_cells = size.cells;
        metrics_.ordered_bytes = size.bytes;
        metrics_.bound_bytes += size.bytes;
        metrics_.bound_active += "+ordered";
      } catch (const std::bad_alloc&) {
        metrics_.ordered_fallback = "allocation failed; omitted ordered index";
      }
    }
    if (!metrics_.ordered_fallback.empty()) {
      if (!metrics_.bound_fallback.empty()) metrics_.bound_fallback += "; ";
      metrics_.bound_fallback += metrics_.ordered_fallback;
    }
    const auto elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
    metrics_.ordered_preprocess_ms += elapsed;
    metrics_.bound_preprocess_ms += elapsed;
    return static_cast<bool>(ordered_);
  }

  void initialize_alphabet() {
    if (alphabet_ready_) return;
    std::array<bool, 256> observed{};
    for (unsigned char c : x_) observed[c] = true;
    for (unsigned char c : y_) observed[c] = true;
    for (unsigned int c = 0; c < 256; ++c)
      if (observed[c]) symbol_[c] = alphabet_++;
    alphabet_ready_ = true;
  }

  std::uint8_t gram_id(const std::string& text, std::size_t offset, int q) const {
    if (q == 1) return static_cast<std::uint8_t>(symbol_[static_cast<unsigned char>(text[offset])]);
    if (alphabet_ <= 4) {
      unsigned int code = 0;
      for (int k = 0; k < q; ++k)
        code = code * alphabet_ + symbol_[static_cast<unsigned char>(text[offset + k])];
      return static_cast<std::uint8_t>(code);
    }
    // A fixed additive histogram sketch: collisions can cancel differences
    // but cannot increase positive or negative histogram mass. Thus hashing
    // can only weaken the lower bound, never invalidate it.
    std::uint32_t hash = 2166136261U;
    for (int k = 0; k < q; ++k) {
      hash ^= static_cast<unsigned char>(text[offset + k]);
      hash *= 16777619U;
    }
    hash ^= hash >> 16;
    hash *= 0x7feb352dU;
    hash ^= hash >> 15;
    return static_cast<std::uint8_t>(hash);
  }

  static void initialize_cursor(const Stage& stage, Cursor& cursor) {
    cursor.i = cursor.j = 0;
    cursor.delta = stage.initial_delta;
    cursor.positive = stage.initial_positive;
    cursor.negative = stage.initial_negative;
  }

  bool ensure_stage(int index) {
    if (attempted_[index]) return static_cast<bool>(stages_[index]);
    attempted_[index] = true;
    const auto start = std::chrono::steady_clock::now();
    const int q = index == 0 ? 1 : index == 1 ? 2 : 4;
    const auto nx = x_.size() >= static_cast<std::size_t>(q) ? x_.size() - q + 1 : 0;
    const auto ny = y_.size() >= static_cast<std::size_t>(q) ? y_.size() - q + 1 : 0;
    const std::uint64_t bytes = sizeof(Stage) + static_cast<std::uint64_t>(nx) + ny;
    if (metrics_.bound_bytes > options_.bound_memory ||
        bytes > options_.bound_memory - metrics_.bound_bytes) {
      if (!metrics_.bound_fallback.empty()) metrics_.bound_fallback += "; ";
      metrics_.bound_fallback += "memory limit omitted q=" + std::to_string(q);
      metrics_.bound_preprocess_ms += std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - start).count();
      return false;
    }
    initialize_alphabet();
    auto stage = std::make_unique<Stage>();
    stage->q = q;
    stage->bins = q == 1 ? static_cast<int>(alphabet_) : 256;
    if (q > 1 && alphabet_ <= 4) {
      stage->bins = 1;
      for (int k = 0; k < q; ++k) stage->bins *= alphabet_;
    }
    stage->xgrams.resize(nx);
    stage->ygrams.resize(ny);
    for (std::size_t i = 0; i < nx; ++i) {
      const auto id = gram_id(x_, i, q);
      stage->xgrams[i] = id;
      ++stage->initial_delta[id];
    }
    for (std::size_t j = 0; j < ny; ++j) {
      const auto id = gram_id(y_, j, q);
      stage->ygrams[j] = id;
      --stage->initial_delta[id];
    }
    for (int b = 0; b < stage->bins; ++b) {
      const int value = stage->initial_delta[b];
      stage->initial_positive += std::max(value, 0);
      stage->initial_negative += std::max(-value, 0);
    }
    for (auto& cursor : stage->cursors) initialize_cursor(*stage, cursor);
    metrics_.bound_bytes += sizeof(Stage) + stage->xgrams.capacity() + stage->ygrams.capacity();
    metrics_.bound_active += index == 0 ? "+composition" : index == 1 ? "+q2" : "+q4";
    stages_[index] = std::move(stage);
    metrics_.bound_preprocess_ms += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
    return true;
  }

  static void update(Cursor& cursor, std::uint8_t id, int amount) {
    const int old = cursor.delta[id];
    const int value = old + amount;
    cursor.delta[id] = value;
    cursor.positive += std::max(value, 0) - std::max(old, 0);
    cursor.negative += std::max(-value, 0) - std::max(-old, 0);
  }

  void move_cursor(const Stage& stage, Cursor& cursor, std::size_t i, std::size_t j) {
    metrics_.bound_cursor_updates += i > cursor.i ? i - cursor.i : cursor.i - i;
    metrics_.bound_cursor_updates += j > cursor.j ? j - cursor.j : cursor.j - j;
    while (cursor.i < i) {
      if (cursor.i < stage.xgrams.size()) update(cursor, stage.xgrams[cursor.i], -1);
      ++cursor.i;
    }
    while (cursor.i > i) {
      --cursor.i;
      if (cursor.i < stage.xgrams.size()) update(cursor, stage.xgrams[cursor.i], 1);
    }
    while (cursor.j < j) {
      if (cursor.j < stage.ygrams.size()) update(cursor, stage.ygrams[cursor.j], 1);
      ++cursor.j;
    }
    while (cursor.j > j) {
      --cursor.j;
      if (cursor.j < stage.ygrams.size()) update(cursor, stage.ygrams[cursor.j], -1);
    }
  }
};

}  // namespace cba
