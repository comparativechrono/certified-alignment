#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace cba {

// Demand-driven threshold proofs for the same relaxation as OrderedBounds.
// An interrupted proof is UNKNOWN, never a lower bound from a partial minimum.
class SparseOrderedBounds {
 public:
  struct Stats {
    std::uint64_t queries = 0, expansions = 0, cache_hits = 0;
    std::uint64_t guard_rejections = 0, cached_states = 0, tiles = 0;
    std::uint64_t bytes = 0, peak_stack = 0;
    double preprocess_ms = 0;
    std::string fallback;
  };

  struct Layout {
    std::size_t blocks = 0, columns = 0, tile_columns = 0;
    std::size_t full_grams = 0, tail_grams = 0;
    std::uint64_t directory_slots = 0, initial_bytes = 0;
  };

  // Payload includes the direct directory, exact packed words, and explicit
  // proof stack. Tiles are charged individually. Object/allocator overhead is
  // excluded. Every product is checked before an allocation can occur.
  static Layout checked_size(std::size_t n, std::size_t m, int block) {
    if (block < 1 || block > 8)
      throw std::invalid_argument("order-block must be between 1 and 8");
    constexpr auto score_limit = static_cast<std::size_t>(
        std::numeric_limits<int>::max() / 4);
    if (n >= score_limit || m >= score_limit || n + m >= score_limit)
      throw std::length_error("sparse bound exceeds the integer score limit");
    Layout size;
    const auto k = static_cast<std::size_t>(block);
    size.blocks = n / k + (n % k != 0);
    size.columns = m + 1;
    size.tile_columns = size.columns / kTileWidth + (size.columns % kTileWidth != 0);
    size.directory_slots = multiply(size.blocks, size.tile_columns);
    if (n >= k && m >= k) size.full_grams = m - k + 1;
    const auto tail = n % k;
    if (tail && m >= tail) size.tail_grams = m - tail + 1;
    size.initial_bytes = add(multiply(size.directory_slots, sizeof(Tile*)),
        add(multiply(add(add(size.full_grams, size.tail_grams), size.blocks),
                     sizeof(std::uint64_t)),
            multiply(size.blocks, sizeof(Frame))));
    if (size.initial_bytes > std::numeric_limits<std::size_t>::max())
      throw std::length_error("sparse bound allocation exceeds address space");
    return size;
  }

  SparseOrderedBounds(const std::string& x, const std::string& y, int block,
                      std::uint64_t max_work, std::uint64_t memory_limit)
      : x_(x), y_(y), n_(x.size()), m_(y.size()),
        block_(static_cast<std::size_t>(block)),
        layout_(checked_size(n_, m_, block)), max_work_(max_work),
        memory_limit_(memory_limit) {}

  ~SparseOrderedBounds() {
    if (directory_)
      for (std::uint64_t i = 0; i < layout_.directory_slots; ++i) delete directory_[i];
  }

  SparseOrderedBounds(const SparseOrderedBounds&) = delete;
  SparseOrderedBounds& operator=(const SparseOrderedBounds&) = delete;

  const Stats& stats() const { return stats_; }

  // A valid lower bound, capped at the requested threshold after a successful
  // proof. A failed/abandoned proof returns only a previously valid lower.
  int query(std::size_t i, std::size_t j, int threshold) {
    if (i > n_ || j > m_) throw std::out_of_range("sparse suffix position");
    ++stats_.queries;
    const auto b = i / block_ + (i % block_ != 0);
    const auto boundary = std::min(n_, b * block_);
    const auto shifted = j + boundary - i;
    const auto target = std::min(m_, shifted);
    const int offset = static_cast<int>(shifted > m_ ? shifted - m_ : 0);
    const int wanted = threshold <= offset ? 0 : threshold - offset;
    if (prove(b, target, wanted) == Verdict::True)
      return std::max(threshold, lower(b, target) + offset);
    return lower(b, target) + offset;
  }

  // Read-only audit hook. Untouched entries in an allocated tile are omitted.
  template <class Visitor> void visit_cached(Visitor visitor) const {
    if (!directory_) return;
    for (std::size_t b = 0; b < layout_.blocks; ++b) {
      for (std::size_t t = 0; t < layout_.tile_columns; ++t) {
        const auto* tile = directory_[b * layout_.tile_columns + t];
        if (!tile) continue;
        for (std::size_t d = 0; d < kTileWidth; ++d) {
          if (!(tile->used & (std::uint64_t{1} << d))) continue;
          visitor(b, t * kTileWidth + d, tile->values[d].lower, tile->values[d].upper);
        }
      }
    }
  }

 private:
  static constexpr std::size_t kTileWidth = 64;
  static constexpr int kInf = std::numeric_limits<int>::max() / 4;
  struct Interval { int lower = 0, upper = kInf; };
  struct Tile {
    std::array<Interval, kTileWidth> values;
    std::uint64_t used = 0;
  };
  struct Frame {
    std::size_t b = 0, j = 0, p = 0;
    Interval* interval = nullptr;
    int threshold = 0, mismatch = 0, next = 0;
  };
  enum class Verdict { True, False, Unknown, Expand };

  const std::string& x_;
  const std::string& y_;
  std::size_t n_, m_, block_;
  Layout layout_;
  std::uint64_t max_work_, memory_limit_;
  Stats stats_;
  bool ready_ = false, stopped_ = false;
  std::unique_ptr<Tile*[]> directory_;
  std::unique_ptr<std::uint64_t[]> source_, full_grams_, tail_grams_;
  std::unique_ptr<Frame[]> stack_;

  static std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
      throw std::length_error("sparse bound size addition overflow");
    return a + b;
  }
  static std::uint64_t multiply(std::uint64_t a, std::uint64_t b) {
    if (a && b > std::numeric_limits<std::uint64_t>::max() / a)
      throw std::length_error("sparse bound size multiplication overflow");
    return a * b;
  }

  int length_lower(std::size_t b, std::size_t j) const {
    const auto rx = n_ - std::min(n_, b * block_), ry = m_ - j;
    return static_cast<int>(rx > ry ? rx - ry : ry - rx);
  }
  int cheap_upper(std::size_t b, std::size_t j) const {
    return std::max(static_cast<int>(layout_.blocks - b), length_lower(b, j));
  }
  bool exact_base(std::size_t b, std::size_t j) const {
    return b == layout_.blocks || j + std::min(block_, n_ - b * block_) > m_;
  }
  Interval* find(std::size_t b, std::size_t j) const {
    if (!directory_ || b == layout_.blocks) return nullptr;
    auto* tile = directory_[b * layout_.tile_columns + j / kTileWidth];
    if (!tile || !(tile->used & (std::uint64_t{1} << (j % kTileWidth)))) return nullptr;
    return &tile->values[j % kTileWidth];
  }
  int lower(std::size_t b, std::size_t j) const {
    const auto* interval = find(b, j);
    return std::max(length_lower(b, j), interval ? interval->lower : 0);
  }
  Verdict known(std::size_t b, std::size_t j, int threshold) {
    const int lb = length_lower(b, j);
    if (lb >= threshold) return Verdict::True;
    if (exact_base(b, j)) return Verdict::False;
    if (cheap_upper(b, j) < threshold) {
      ++stats_.guard_rejections;
      return Verdict::False;
    }
    if (const auto* interval = find(b, j)) {
      if (interval->lower >= threshold) {
        ++stats_.cache_hits;
        return Verdict::True;
      }
      if (interval->upper < threshold) {
        ++stats_.cache_hits;
        return Verdict::False;
      }
    }
    return Verdict::Expand;
  }

  void stop(const std::string& reason) {
    stopped_ = true;
    if (stats_.fallback.empty()) stats_.fallback = reason;
  }

  static std::unique_ptr<std::uint64_t[]> pack_grams(const std::string& text,
                                                   std::size_t q,
                                                   std::size_t count) {
    if (!count) return {};
    auto result = std::unique_ptr<std::uint64_t[]>(new std::uint64_t[count]);
    const auto mask = q == 8 ? std::numeric_limits<std::uint64_t>::max()
                            : (std::uint64_t{1} << (8 * q)) - 1;
    std::uint64_t code = 0;
    for (std::size_t j = 0; j < text.size(); ++j) {
      code = ((code << 8) | static_cast<unsigned char>(text[j])) & mask;
      if (j + 1 >= q) result[j + 1 - q] = code;
    }
    return result;
  }

  bool ensure_index() {
    if (ready_) return true;
    if (stopped_) return false;
    const auto start = std::chrono::steady_clock::now();
    if (layout_.initial_bytes > memory_limit_) {
      stop("shared bound-memory limit omitted sparse index");
    } else {
      try {
        // Temporaries make failed allocations leave no partially active index.
        auto directory = std::unique_ptr<Tile*[]>(
            new Tile*[static_cast<std::size_t>(layout_.directory_slots)]());
        auto source = std::unique_ptr<std::uint64_t[]>(new std::uint64_t[layout_.blocks]);
        auto frames = std::unique_ptr<Frame[]>(new Frame[layout_.blocks]);
        auto full = pack_grams(y_, block_, layout_.full_grams);
        auto tail = pack_grams(y_, n_ % block_, layout_.tail_grams);
        for (std::size_t b = 0; b < layout_.blocks; ++b) {
          const auto begin = b * block_, len = std::min(block_, n_ - begin);
          std::uint64_t code = 0;
          for (std::size_t d = 0; d < len; ++d)
            code = (code << 8) | static_cast<unsigned char>(x_[begin + d]);
          source[b] = code;
        }
        directory_ = std::move(directory);
        source_ = std::move(source);
        stack_ = std::move(frames);
        full_grams_ = std::move(full);
        tail_grams_ = std::move(tail);
        stats_.bytes = layout_.initial_bytes;
        ready_ = true;
      } catch (const std::bad_alloc&) {
        stop("allocation failed; omitted sparse index");
      }
    }
    stats_.preprocess_ms += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
    return ready_;
  }

  Interval* acquire(std::size_t b, std::size_t j) {
    auto*& tile = directory_[b * layout_.tile_columns + j / kTileWidth];
    if (!tile) {
      if (stats_.bytes > memory_limit_ || sizeof(Tile) > memory_limit_ - stats_.bytes) {
        stop("shared bound-memory limit stopped sparse proof");
        return nullptr;
      }
      try { tile = new Tile; }
      catch (const std::bad_alloc&) {
        stop("allocation failed; stopped sparse proof");
        return nullptr;
      }
      stats_.bytes += sizeof(Tile);
      ++stats_.tiles;
    }
    const auto mask = std::uint64_t{1} << (j % kTileWidth);
    if (!(tile->used & mask)) {
      tile->used |= mask;
      ++stats_.cached_states;
    }
    return &tile->values[j % kTileWidth];
  }

  Verdict prove(std::size_t b, std::size_t j, int threshold) {
    std::size_t depth = 0;
    while (true) {
      auto verdict = known(b, j, threshold);
      if (verdict == Verdict::Expand) {
        if (stopped_) return Verdict::Unknown;
        if (stats_.expansions >= max_work_) {
          stop("sparse-max-work limit stopped sparse proof");
          return Verdict::Unknown;
        }
        if (!ensure_index()) return Verdict::Unknown;
        auto* interval = acquire(b, j);
        if (!interval) return Verdict::Unknown;
        ++stats_.expansions;
        const auto len = std::min(block_, n_ - b * block_);
        const auto* grams = len == block_ ? full_grams_.get() : tail_grams_.get();
        auto& frame = stack_[depth++];
        frame = {b, j, j + len, interval, threshold, source_[b] != grams[j], 0};
        stats_.peak_stack = std::max(stats_.peak_stack, static_cast<std::uint64_t>(depth));
        b = frame.b + 1;
        j = frame.p;
        threshold = frame.threshold - frame.mismatch;
        continue;
      }
      // A proven child false disproves the minimum threshold. A true child
      // advances to the next candidate; only all true children prove the node.
      while (depth) {
        auto& frame = stack_[depth - 1];
        if (verdict == Verdict::False) {
          frame.interval->upper = std::min(frame.interval->upper, frame.threshold - 1);
          --depth;
          continue;
        }
        ++frame.next;
        if (frame.next < (frame.p < m_ ? 3 : 2)) {
          b = frame.b + 1;
          j = frame.next == 1 ? frame.p - 1 : frame.p + 1;
          threshold = frame.threshold - 1;
          break;
        }
        frame.interval->lower = std::max(frame.interval->lower, frame.threshold);
        --depth;
      }
      if (!depth) return verdict;
    }
  }
};

}  // namespace cba
