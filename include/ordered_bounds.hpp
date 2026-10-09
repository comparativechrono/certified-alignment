#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace cba {

// Dense ordered block relaxation. Every possible target block placement is
// implicit in the rows; no selected seed match is required or forced.
// A block pays its target-length displacement, or one for an unequal
// equal-length chunk. See ORDERED_BOUND.md for the relaxation and proof.
class OrderedBounds {
 public:
  struct Layout {
    std::size_t blocks = 0, columns = 0;
    std::size_t full_grams = 0, tail_grams = 0;
    std::uint64_t cells = 0, bytes = 0;
  };

  // bytes counts this object and requested array payload, excluding allocator
  // bookkeeping. The packed target gram arrays are retained and included.
  static Layout checked_size(std::size_t n, std::size_t m, int block) {
    if (block < 1 || block > 8)
      throw std::invalid_argument("order-block must be between 1 and 8");
    constexpr auto score_limit = static_cast<std::size_t>(
        std::numeric_limits<int>::max() / 4);
    if (n >= score_limit || m >= score_limit || n + m >= score_limit)
      throw std::length_error("ordered bound exceeds the integer score limit");
    Layout size;
    const auto k = static_cast<std::size_t>(block);
    size.blocks = n / k + (n % k != 0);
    size.columns = m + 1;
    size.cells = multiply(size.blocks + 1, size.columns);
    if (n >= k && m >= k) size.full_grams = m - k + 1;
    const auto tail = n % k;
    if (tail && m >= tail) size.tail_grams = m - tail + 1;
    const auto score_bytes = multiply(size.cells, sizeof(int));
    const auto gram_bytes = multiply(add(size.full_grams, size.tail_grams),
                                     sizeof(std::uint64_t));
    size.bytes = add(add(sizeof(OrderedBounds), score_bytes), gram_bytes);
    if (size.bytes > std::numeric_limits<std::size_t>::max())
      throw std::length_error("ordered bound allocation exceeds address space");
    return size;
  }

  OrderedBounds(const std::string& x, const std::string& y, int block)
      : n_(x.size()), m_(y.size()), block_(static_cast<std::size_t>(block)),
        layout_(checked_size(n_, m_, block)),
        scores_(new int[static_cast<std::size_t>(layout_.cells)]),
        full_grams_(pack_grams(y, block_, layout_.full_grams)),
        tail_grams_(pack_grams(y, n_ % block_, layout_.tail_grams)) {
    int* terminal = scores_.get() + layout_.blocks * layout_.columns;
    for (std::size_t j = 0; j <= m_; ++j)
      terminal[j] = static_cast<int>(m_ - j);
    for (std::size_t b = layout_.blocks; b-- > 0;) {
      const auto start = b * block_;
      const auto len = std::min(block_, n_ - start);
      std::uint64_t source = 0;
      for (std::size_t d = 0; d < len; ++d)
        source = (source << 8) | static_cast<unsigned char>(x[start + d]);
      const auto* grams = len == block_ ? full_grams_.get() : tail_grams_.get();
      int* row = scores_.get() + b * layout_.columns;
      const int* next = row + layout_.columns;
      for (std::size_t j = 0; j <= m_; ++j) {
        const auto p = j + len;
        if (p > m_) {
          row[j] = static_cast<int>(p - m_) + next[m_];
        } else {
          // The full minimization over target chunk lengths reduces to the
          // exact equal-length choice and its two neighbors: next is
          // 1-Lipschitz. len is positive, so p-1 is always defined.
          int value = std::min(next[p] + (source != grams[j]), 1 + next[p - 1]);
          if (p < m_) value = std::min(value, 1 + next[p + 1]);
          row[j] = value;
        }
      }
    }
  }

  int query(std::size_t i, std::size_t j) const {
    if (i > n_ || j > m_) throw std::out_of_range("ordered suffix position");
    const auto b = i / block_ + (i % block_ != 0);
    const auto boundary = std::min(n_, b * block_);
    const auto target = j + boundary - i;
    // Relax the content of a partial source block, preserving only length.
    return scores_[b * layout_.columns + std::min(m_, target)] +
           static_cast<int>(target > m_ ? target - m_ : 0);
  }

 private:
  std::size_t n_, m_, block_;
  Layout layout_;
  std::unique_ptr<int[]> scores_;
  std::unique_ptr<std::uint64_t[]> full_grams_, tail_grams_;

  static std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    if (b > std::numeric_limits<std::uint64_t>::max() - a)
      throw std::length_error("ordered bound size addition overflow");
    return a + b;
  }

  static std::uint64_t multiply(std::uint64_t a, std::uint64_t b) {
    if (a && b > std::numeric_limits<std::uint64_t>::max() / a)
      throw std::length_error("ordered bound size multiplication overflow");
    return a * b;
  }

  static std::unique_ptr<std::uint64_t[]> pack_grams(const std::string& text,
                                                   std::size_t q,
                                                   std::size_t count) {
    if (!count) return {};
    auto codes = std::unique_ptr<std::uint64_t[]>(new std::uint64_t[count]);
    const auto mask = q == 8 ? std::numeric_limits<std::uint64_t>::max()
                            : (std::uint64_t{1} << (8 * q)) - 1;
    std::uint64_t code = 0;
    for (std::size_t j = 0; j < text.size(); ++j) {
      code = ((code << 8) | static_cast<unsigned char>(text[j])) & mask;
      if (j + 1 >= q) codes[j + 1 - q] = code;
    }
    return codes;
  }
};

}  // namespace cba
