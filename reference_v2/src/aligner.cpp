#include "engine.hpp"
#include "lower_bounds.hpp"
#include "predictor.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cba {
namespace {

using Clock = std::chrono::steady_clock;
constexpr int kInf = std::numeric_limits<int>::max() / 4;

double elapsed_ms(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

struct Row {
  std::size_t lo = 0, hi = 0, offset = 0;
};

struct Band {
  std::vector<Row> rows;
  std::vector<int> scores;

  bool contains(std::size_t i, std::size_t j) const {
    return i < rows.size() && j >= rows[i].lo && j <= rows[i].hi;
  }

  int get(std::size_t i, std::size_t j) const {
    if (!contains(i, j)) return kInf;
    const auto& row = rows[i];
    return scores[row.offset + j - row.lo];
  }

  int& at(std::size_t i, std::size_t j) {
    const auto& row = rows[i];
    return scores[row.offset + j - row.lo];
  }
};

std::string traceback(const Band& band, const std::string& x,
                      const std::string& y) {
  std::size_t i = x.size(), j = y.size();
  std::vector<char> reverse_ops;
  reverse_ops.reserve(i + j);
  while (i || j) {
    const int current = band.get(i, j);
    // Prefer a diagonal, then deletion, then insertion on equal scores.
    if (i && j && band.get(i - 1, j - 1) < kInf &&
        current == band.get(i - 1, j - 1) + (x[i - 1] != y[j - 1])) {
      reverse_ops.push_back(x[i - 1] == y[j - 1] ? '=' : 'X');
      --i;
      --j;
    } else if (i && band.get(i - 1, j) < kInf &&
               current == band.get(i - 1, j) + 1) {
      reverse_ops.push_back('D');
      --i;
    } else if (j && band.get(i, j - 1) < kInf &&
               current == band.get(i, j - 1) + 1) {
      reverse_ops.push_back('I');
      --j;
    } else {
      throw std::logic_error("traceback has no valid predecessor");
    }
  }
  std::reverse(reverse_ops.begin(), reverse_ops.end());
  std::string cigar;
  for (std::size_t begin = 0; begin < reverse_ops.size();) {
    std::size_t end = begin + 1;
    while (end < reverse_ops.size() && reverse_ops[end] == reverse_ops[begin]) ++end;
    cigar += std::to_string(end - begin);
    cigar.push_back(reverse_ops[begin]);
    begin = end;
  }
  return cigar;
}

}  // namespace

Result align(const std::string& x, const std::string& y,
             const Prediction& prediction, const Options& options) {
  Result result;
  const std::size_t n = x.size(), m = y.size();
  if (n >= static_cast<std::size_t>(kInf) || m >= static_cast<std::size_t>(kInf) ||
      n + m >= static_cast<std::size_t>(kInf)) {
    result.error = "sequence lengths exceed the integer score limit";
    return result;
  }
  if (prediction.centers.size() != n + 1) {
    result.error = "predictor must supply n+1 row centers";
    return result;
  }
  for (double center : prediction.centers) {
    if (!std::isfinite(center)) {
      result.error = "predictor supplied a non-finite row center";
      return result;
    }
  }
  if (!options.epsilon_den) {
    result.error = "epsilon denominator must be positive";
    return result;
  }
  if (options.bound != "length" && options.bound != "composition" && options.bound != "qgram") {
    result.error = "unknown suffix bound " + options.bound;
    return result;
  }
  SuffixBounds suffix_bounds(x, y, options, result);

  const auto max_dim = static_cast<std::uint64_t>(std::max(n, m));
  std::uint64_t radius = options.mode == "full" ? max_dim : std::min(options.width, max_dim);
  while (true) {
    const auto dp_start = Clock::now();
    Band band;
    if (result.dp_cells > options.max_cells ||
        n + 1 > options.max_cells - result.dp_cells) {
      result.dp_ms += elapsed_ms(dp_start);
      result.error = "max-cells budget exceeded before the next DP round";
      return result;
    }
    band.rows.resize(n + 1);
    std::uint64_t cells = 0;
    const bool full = options.mode == "full" || radius >= max_dim;
    for (std::size_t i = 0; i <= n; ++i) {
      auto& row = band.rows[i];
      if (full) {
        row.lo = 0;
        row.hi = m;
      } else {
        const double center = std::clamp(prediction.centers[i], 0.0, static_cast<double>(m));
        const double lower = std::max(0.0, std::floor(center - static_cast<double>(radius)));
        const double upper = std::min(static_cast<double>(m), std::ceil(center + static_cast<double>(radius)));
        row.lo = static_cast<std::size_t>(lower);
        row.hi = static_cast<std::size_t>(upper);
      }
      // Endpoints are always included, including the n=0 single-row case.
      if (i == 0) row.lo = 0;
      if (i == n) row.hi = m;
      row.offset = static_cast<std::size_t>(cells);
      const auto count = static_cast<std::uint64_t>(row.hi - row.lo + 1);
      if (result.dp_cells > options.max_cells || cells > options.max_cells - result.dp_cells ||
          count > options.max_cells - result.dp_cells - cells) {
        result.dp_ms += elapsed_ms(dp_start);
        result.error = "max-cells budget exceeded before the next DP round";
        result.certified_exact = false;
        result.certified_approx = false;
        return result;
      }
      cells += count;
    }
    if (cells > std::numeric_limits<std::size_t>::max() / sizeof(int)) {
      result.dp_ms += elapsed_ms(dp_start);
      result.error = "DP allocation would overflow the address space";
      return result;
    }
    band.scores.assign(static_cast<std::size_t>(cells), kInf);
    for (std::size_t i = 0; i <= n; ++i) {
      const auto& row = band.rows[i];
      for (std::size_t j = row.lo; j <= row.hi; ++j) {
        if (!i && !j) {
          band.at(i, j) = 0;
          continue;
        }
        int value = kInf;
        if (i && j) value = std::min(value, band.get(i - 1, j - 1) + (x[i - 1] != y[j - 1]));
        if (i) value = std::min(value, band.get(i - 1, j) + 1);
        if (j) value = std::min(value, band.get(i, j - 1) + 1);
        band.at(i, j) = value;
      }
    }
    result.dp_ms += elapsed_ms(dp_start);
    result.dp_cells += cells;
    result.peak_cells = std::max(result.peak_cells, cells);
    ++result.rounds;

    const auto certificate_start = Clock::now();
    const double previous_preprocess_ms = result.bound_preprocess_ms;
    // U is the restricted optimum, not merely the cost of an arbitrary path.
    const int upper = band.get(n, m);
    int lower = upper;
    const bool strong_bound = options.bound != "length" && upper < kInf;
    if (strong_bound) suffix_bounds.reset_round();
    auto inspect_exit = [&](std::size_t i, std::size_t j, int prefix, int weight,
                            ExitStream stream) {
      if (band.contains(i, j)) return;
      ++result.certificate_edges;
      const auto remain_x = static_cast<long long>(n - i);
      const auto remain_y = static_cast<long long>(m - j);
      const int suffix_bound = strong_bound
          ? suffix_bounds.query(i, j, lower - prefix - weight, stream)
          : static_cast<int>(std::llabs(remain_x - remain_y));
      lower = std::min(lower, prefix + weight + suffix_bound);
    };
    for (std::size_t i = 0; i <= n; ++i) {
      const auto& row = band.rows[i];
      // A horizontal edge exits only at the right end of an interval.
      if (row.hi < m && band.get(i, row.hi) < kInf)
        inspect_exit(i, row.hi + 1, band.get(i, row.hi), 1, ExitStream::Horizontal);
      if (i == n) continue;

      const auto& next = band.rows[i + 1];
      // Enumerate interval differences, rather than scanning interior cells.
      // Signed endpoints make empty ranges at columns 0 and m unambiguous.
      auto exit_range = [&](long long first, long long last, bool diagonal, ExitStream stream) {
        first = std::max(first, static_cast<long long>(row.lo));
        last = std::min(last, static_cast<long long>(row.hi));
        if (diagonal) last = std::min(last, static_cast<long long>(m) - 1);
        for (long long column = first; column <= last; ++column) {
          const auto j = static_cast<std::size_t>(column);
          const int prefix = band.get(i, j);
          if (prefix >= kInf) continue;
          if (diagonal) inspect_exit(i + 1, j + 1, prefix, x[i] != y[j], stream);
          else inspect_exit(i + 1, j, prefix, 1, stream);
        }
      };
      // Deletion leaves j fixed; a diagonal advances j by one.
      exit_range(row.lo, static_cast<long long>(next.lo) - 1, false, ExitStream::VerticalLeft);
      exit_range(static_cast<long long>(next.hi) + 1, row.hi, false, ExitStream::VerticalRight);
      exit_range(row.lo, static_cast<long long>(next.lo) - 2, true, ExitStream::DiagonalLeft);
      exit_range(next.hi, row.hi, true, ExitStream::DiagonalRight);
    }
    result.certificate_ms += std::max(0.0, elapsed_ms(certificate_start) -
        (result.bound_preprocess_ms - previous_preprocess_ms));
    result.distance = upper < kInf ? upper : -1;
    result.lower_bound = lower < kInf ? lower : -1;
    result.trace.push_back({radius, cells, result.distance, result.lower_bound});

    const bool finite = upper < kInf && lower < kInf;
    const bool exact = finite && lower >= upper;
    const bool approximate = finite &&
        static_cast<unsigned __int128>(upper) * options.epsilon_den <=
        static_cast<unsigned __int128>(lower) *
            (static_cast<unsigned __int128>(options.epsilon_den) + options.epsilon_num);
    if (exact || approximate) {
      const auto traceback_start = Clock::now();
      result.cigar = traceback(band, x, y);
      result.dp_ms += elapsed_ms(traceback_start);
      result.certified_exact = exact;
      result.certified_approx = approximate;
      return result;
    }
    if (full) {
      result.error = "internal error: full DP failed the exact certificate";
      return result;
    }
    radius = radius == 0 ? 1 : radius > max_dim / 2 ? max_dim : radius * 2;
  }
}

}  // namespace cba

namespace {

std::string json_string(const std::string& input) {
  std::ostringstream out;
  out << '"';
  for (unsigned char c : input) {
    switch (c) {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\b': out << "\\b"; break;
      case '\f': out << "\\f"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (c < 0x20) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<int>(c) << std::dec;
        } else {
          out << static_cast<char>(c);
        }
    }
  }
  out << '"';
  return out.str();
}

std::uint64_t parse_unsigned(const std::string& text, const std::string& name) {
  if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
    throw std::invalid_argument(name + " must be a nonnegative integer");
  std::size_t used = 0;
  const auto value = std::stoull(text, &used);
  if (used != text.size()) throw std::invalid_argument("invalid " + name);
  return value;
}

std::pair<std::uint64_t, std::uint64_t> parse_epsilon(std::string text) {
  constexpr std::uint64_t limit = 1000000000000000000ULL;
  if (text.empty() || text[0] == '-')
    throw std::invalid_argument("epsilon must be a nonnegative decimal");
  if (text[0] == '+') text.erase(0, 1);
  int exponent = 0;
  const auto ep = text.find_first_of("eE");
  if (ep != std::string::npos) {
    const std::string exp_text = text.substr(ep + 1);
    std::size_t used = 0;
    exponent = std::stoi(exp_text, &used);
    if (used != exp_text.size() || exponent < -18 || exponent > 18)
      throw std::invalid_argument("epsilon exponent must be between -18 and 18");
    text.resize(ep);
  }
  if (text.empty() || text.find_first_of("0123456789") == std::string::npos ||
      text.find_first_not_of("0123456789.") != std::string::npos ||
      (text.find('.') != std::string::npos && text.find('.', text.find('.') + 1) != std::string::npos))
    throw std::invalid_argument("epsilon must be a nonnegative decimal");
  const auto dot = text.find('.');
  int decimal_places = 0;
  if (dot != std::string::npos) {
    decimal_places = static_cast<int>(text.size() - dot - 1);
    text.erase(dot, 1);
  }
  // Trim zeroes before conversion, retaining at least one digit.
  while (text.size() > 1 && text.front() == '0') text.erase(0, 1);
  while (decimal_places > 0 && !text.empty() && text.back() == '0') {
    text.pop_back();
    --decimal_places;
  }
  if (text.empty()) text = "0";
  std::uint64_t numerator = parse_unsigned(text, "epsilon");
  if (numerator > limit) throw std::invalid_argument("epsilon has too many significant digits");
  if (!numerator) return {0, 1};
  const int power = decimal_places - exponent;
  if (power < -18 || power > 18)
    throw std::invalid_argument("epsilon requires more than 18 decimal places");
  std::uint64_t denominator = 1;
  for (int i = 0; i < power; ++i) denominator *= 10;
  for (int i = 0; i < -power; ++i) {
    if (numerator > limit / 10) throw std::invalid_argument("epsilon is too large");
    numerator *= 10;
  }
  const auto divisor = std::gcd(numerator, denominator);
  return {numerator / divisor, denominator / divisor};
}

void score_json(std::ostream& out, int score) {
  if (score < 0) out << "null";
  else out << score;
}

void print_result(const std::string& id, const std::string& x, const std::string& y,
                  const cba::Options& options, const cba::Prediction& prediction,
                  const cba::Result& result, double predictor_ms, double total_ms) {
  auto& out = std::cout;
  out << std::fixed << std::setprecision(6)
      << "{\"id\":" << json_string(id) << ",\"mode\":" << json_string(options.mode)
      << ",\"bound\":" << json_string(options.bound)
      << ",\"n\":" << x.size() << ",\"m\":" << y.size() << ",\"distance\":";
  score_json(out, result.distance);
  out << ",\"lower_bound\":";
  score_json(out, result.lower_bound);
  out << ",\"certified_exact\":" << (result.certified_exact ? "true" : "false")
      << ",\"certified_approx\":" << (result.certified_approx ? "true" : "false")
      << ",\"cigar\":" << json_string(result.cigar)
      << ",\"rounds\":" << result.rounds << ",\"dp_cells\":" << result.dp_cells
      << ",\"certificate_edges\":" << result.certificate_edges
      << ",\"probe_cells\":" << prediction.probe_cells
      << ",\"predictor_ms\":" << predictor_ms << ",\"dp_ms\":" << result.dp_ms
      << ",\"certificate_ms\":" << result.certificate_ms << ",\"total_ms\":" << total_ms
      << ",\"bound_preprocess_ms\":" << result.bound_preprocess_ms
      << ",\"bound_bytes\":" << result.bound_bytes
      << ",\"bound_memory_limit\":" << options.bound_memory
      << ",\"bound_calls\":" << result.bound_calls
      << ",\"composition_calls\":" << result.composition_calls
      << ",\"q2_calls\":" << result.q2_calls << ",\"q4_calls\":" << result.q4_calls
      << ",\"bound_cursor_updates\":" << result.bound_cursor_updates
      << ",\"bound_active\":" << json_string(result.bound_active)
      << ",\"bound_fallback\":" << json_string(result.bound_fallback)
      << ",\"peak_cells\":" << result.peak_cells << ",\"trace\":[";
  for (std::size_t i = 0; i < result.trace.size(); ++i) {
    if (i) out << ',';
    const auto& round = result.trace[i];
    out << "{\"width\":" << round.width << ",\"cells\":" << round.cells << ",\"upper\":";
    score_json(out, round.upper);
    out << ",\"lower\":";
    score_json(out, round.lower);
    out << '}';
  }
  out << ']';
  if (!result.error.empty()) out << ",\"error\":" << json_string(result.error);
  out << "}\n";
}

}  // namespace

int main(int argc, char** argv) {
  cba::Options options;
  try {
    for (int i = 1; i < argc; ++i) {
      const std::string option = argv[i];
      if (option == "--help") {
        std::cout << "certalign [--mode full|linear|cost|anchor-linear|anchor-cost] "
                     "[--width 8] [--block 64] [--seed-k 11] [--epsilon 0] "
                     "[--bound length|composition|qgram] "
                     "[--max-cells 50000000]\n"
                     "Read id<TAB>x<TAB>y records from stdin; emit JSONL. "
                     "max-cells limits cumulative DP work per record.\n";
        return 0;
      }
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + option);
      const std::string value = argv[++i];
      if (option == "--mode") options.mode = value;
      else if (option == "--bound") options.bound = value;
      else if (option == "--width") options.width = parse_unsigned(value, "width");
      else if (option == "--max-cells") options.max_cells = parse_unsigned(value, "max-cells");
      else if (option == "--block" || option == "--seed-k") {
        const auto parsed = parse_unsigned(value, option);
        if (!parsed || parsed > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
          throw std::invalid_argument(option + " must be a positive int");
        if (option == "--block") options.block = static_cast<int>(parsed);
        else options.seed_k = static_cast<int>(parsed);
      } else if (option == "--epsilon") {
        const auto rational = parse_epsilon(value);
        options.epsilon_num = rational.first;
        options.epsilon_den = rational.second;
      } else throw std::invalid_argument("unknown option " + option);
    }
    if (options.mode != "full" && options.mode != "linear" && options.mode != "cost" &&
        options.mode != "anchor-linear" && options.mode != "anchor-cost")
      throw std::invalid_argument("unknown mode " + options.mode);
    if (options.bound != "length" && options.bound != "composition" && options.bound != "qgram")
      throw std::invalid_argument("unknown suffix bound " + options.bound);
  } catch (const std::exception& error) {
    std::cerr << "certalign: " << error.what() << '\n';
    return 2;
  }

  std::string line;
  while (std::getline(std::cin, line)) {
    const auto total_start = std::chrono::steady_clock::now();
    if (!line.empty() && line.back() == '\r') line.pop_back();
    std::string id, x, y;
    cba::Prediction prediction;
    cba::Result result;
    double predictor_ms = 0;
    try {
      const auto first = line.find('\t');
      const auto second = first == std::string::npos ? std::string::npos : line.find('\t', first + 1);
      id = line.substr(0, first);
      if (first == std::string::npos || second == std::string::npos ||
          line.find('\t', second + 1) != std::string::npos)
        throw std::invalid_argument("expected exactly three TSV fields: id, x, y");
      x = line.substr(first + 1, second - first - 1);
      y = line.substr(second + 1);
      const auto predictor_start = std::chrono::steady_clock::now();
      prediction = cba::predict(x, y, options.mode, options.block, options.seed_k);
      predictor_ms = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - predictor_start).count();
      result = cba::align(x, y, prediction, options);
    } catch (const std::exception& error) {
      result.error = error.what();
      result.certified_exact = false;
      result.certified_approx = false;
    }
    const double total_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - total_start).count();
    print_result(id, x, y, options, prediction, result, predictor_ms, total_ms);
  }
  return std::cin.bad() ? 1 : 0;
}
