// Native, exact global edit-distance baselines. No score bounds or oracle input.
#include "edlib.h"
#include "engine.hpp"
#include "predictor.hpp"
extern "C" {
#include "wavefront/wavefront_align.h"
#include "wavefront/wfa.h"
}

#include <chrono>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Clock = std::chrono::steady_clock;
static double elapsed(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
struct Pair { std::string id, x, y; };
struct Result {
  int distance = -1;
  int lower_bound = -1;
  std::string cigar;
  double total_ms = 0;
  std::uint64_t rounds = 0, dp_cells = 0, sparse_expansions = 0, bound_bytes = 0;
  double certificate_ms = 0, bound_preprocess_ms = 0;
  bool certified_exact = true;
  std::string error, sparse_fallback;
  std::vector<cba::RoundTrace> trace;
};

static std::string json_string(const std::string& value) {
  const char hex[] = "0123456789abcdef";
  std::string out = "\"";
  for (unsigned char c : value) {
    if (c == '"' || c == '\\') { out += '\\'; out += char(c); }
    else if (c < 32 || c >= 127) {
      out += "\\u00"; out += hex[c >> 4]; out += hex[c & 15];
    } else out += char(c);
  }
  return out + '"';
}

template<class Operation> static std::string compact_cigar(int length, Operation op) {
  std::string out;
  for (int begin = 0; begin < length;) {
    const char c = op(begin);
    int end = begin + 1;
    while (end < length && op(end) == c) ++end;
    out += std::to_string(end - begin);
    out += c;
    begin = end;
  }
  return out;
}

static Result align_edlib(const Pair& pair) {
  const auto start = Clock::now();
  const auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
  auto raw = edlibAlign(pair.x.data(), int(pair.x.size()),
                        pair.y.data(), int(pair.y.size()), config);
  if (raw.status != EDLIB_STATUS_OK || raw.editDistance < 0) {
    edlibFreeAlignResult(raw);
    throw std::runtime_error("Edlib did not return an exact alignment");
  }
  // Edlib INSERT consumes query (x), hence D in our x->y convention.
  const int distance = raw.editDistance;
  auto cigar = compact_cigar(raw.alignmentLength, [&](int pos) {
    switch (raw.alignment[pos]) {
      case EDLIB_EDOP_MATCH: return '=';
      case EDLIB_EDOP_MISMATCH: return 'X';
      case EDLIB_EDOP_INSERT: return 'D';
      case EDLIB_EDOP_DELETE: return 'I';
      default: throw std::runtime_error("invalid Edlib path operation");
    }
  });
  // Edlib returns no alignment operations on its empty-input fast path.
  // The unique global path here consists only of insertions or deletions.
  if (pair.x.empty() && !pair.y.empty()) cigar = std::to_string(pair.y.size()) + 'I';
  if (pair.y.empty() && !pair.x.empty()) cigar = std::to_string(pair.x.size()) + 'D';
  edlibFreeAlignResult(raw);
  Result result;
  result.distance = result.lower_bound = distance;
  result.cigar = std::move(cigar);
  result.total_ms = elapsed(start);
  return result;
}

static Result align_wfa(const Pair& pair, wavefront_aligner_t* workspace) {
  const auto start = Clock::now();
  const int status = wavefront_align(workspace, pair.x.data(), int(pair.x.size()),
                                     pair.y.data(), int(pair.y.size()));
  if (status != WF_STATUS_ALG_COMPLETED)
    throw std::runtime_error("WFA2 failed with status " + std::to_string(status));
  const auto* raw = workspace->cigar;
  auto cigar = compact_cigar(raw->end_offset - raw->begin_offset, [&](int pos) {
    const char op = raw->operations[raw->begin_offset + pos];
    if (op == 'M' || op == '=') return '=';
    if (op == 'X' || op == 'I' || op == 'D') return op;
    throw std::runtime_error("invalid WFA2 path operation");
  });
  // Use the public edit scorer, as in WFA2's benchmark utilities. In ultralow
  // mode its empty-input fast path can leave cigar->score unset (INT_MIN).
  Result result;
  result.distance = result.lower_bound = cigar_score_edit(raw);
  result.cigar = std::move(cigar);
  result.total_ms = elapsed(start);
  return result;
}

static Result align_prototype(const Pair& pair, const cba::Options& options) {
  const auto start = Clock::now();
  Result result;
  {
    // These are the frozen v5 engine and predictor, compiled without source changes.
    auto prediction = cba::predict(pair.x, pair.y, options.mode, options.block, options.seed_k);
    auto raw = cba::align(pair.x, pair.y, prediction, options);
    result.distance = raw.distance;
    result.lower_bound = raw.lower_bound;
    result.cigar = std::move(raw.cigar);
    result.rounds = raw.rounds;
    result.dp_cells = raw.dp_cells;
    result.sparse_expansions = raw.sparse_expansions;
    result.bound_bytes = raw.bound_bytes;
    result.certificate_ms = raw.certificate_ms;
    result.bound_preprocess_ms = raw.bound_preprocess_ms;
    result.certified_exact = raw.certified_exact;
    result.error = std::move(raw.error);
    result.sparse_fallback = std::move(raw.sparse_fallback);
    result.trace = std::move(raw.trace);
  }  // Include destruction of per-alignment work structures in the measurement.
  result.total_ms = elapsed(start);
  return result;
}

static std::vector<Pair> read_pairs(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open input " + path);
  std::vector<Pair> pairs;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto first = line.find('\t');
    const auto second = first == std::string::npos ? first : line.find('\t', first + 1);
    if (first == std::string::npos || second == std::string::npos ||
        line.find('\t', second + 1) != std::string::npos)
      throw std::runtime_error("expected id<TAB>x<TAB>y on every input line");
    Pair pair{line.substr(0, first), line.substr(first + 1, second - first - 1),
              line.substr(second + 1)};
    if (pair.x.size() >= INT_MAX / 4 || pair.y.size() >= INT_MAX / 4)
      throw std::runtime_error("sequence exceeds native integer limits");
    pairs.push_back(std::move(pair));
  }
  if (pairs.empty()) throw std::runtime_error("empty batch");
  return pairs;
}

int main(int argc, char** argv) {
  try {
    std::string method, input;
    int repeats = 1, warmups = 0;
    for (int i = 1; i < argc; ++i) {
      const std::string key = argv[i];
      if (i + 1 == argc) throw std::runtime_error("missing option value");
      const std::string value = argv[++i];
      if (key == "--method") method = value;
      else if (key == "--input") input = value;
      else if (key == "--repeats") repeats = std::stoi(value);
      else if (key == "--warmups") warmups = std::stoi(value);
      else throw std::runtime_error("unknown option " + key);
    }
    const bool is_wfa = method == "wfa2-high" || method == "wfa2-ultralow";
    const bool is_prototype = method == "prototype-length" || method == "prototype-qgram" ||
                              method == "prototype-sparse";
    if (method != "edlib" && !is_wfa && !is_prototype)
      throw std::runtime_error("unknown --method");
    if (repeats < 1 || warmups < 0 || repeats > 1000 || warmups > 1000)
      throw std::runtime_error("invalid repetition count");
    auto pairs = read_pairs(input);
    std::vector<Result> results;
    results.reserve(pairs.size() * std::size_t(repeats));
    wavefront_aligner_t* workspace = nullptr;
    const auto wall_start = Clock::now();
    const auto setup_start = Clock::now();
    cba::Options prototype_options;
    if (is_prototype) prototype_options.bound = method.substr(10);
    if (is_wfa) {
      auto attributes = wavefront_aligner_attr_default;
      attributes.distance_metric = edit;
      attributes.alignment_scope = compute_alignment;
      attributes.alignment_form.span = alignment_end2end;
      attributes.alignment_form.extension = false;
      attributes.heuristic.strategy = wf_heuristic_none;
      attributes.memory_mode = method == "wfa2-high" ?
          wavefront_memory_high : wavefront_memory_ultralow;
      attributes.system.max_num_threads = 1;
      workspace = wavefront_aligner_new(&attributes);
      if (!workspace) throw std::runtime_error("cannot create WFA2 workspace");
    }
    const double setup_ms = elapsed(setup_start);
    double warmup_ms = 0, alignment_sum_ms = 0;
    for (int repeat = -warmups; repeat < repeats; ++repeat) {
      for (const auto& pair : pairs) {
        auto result = method == "edlib" ? align_edlib(pair) :
            (is_wfa ? align_wfa(pair, workspace) : align_prototype(pair, prototype_options));
        if (repeat >= 0) {
          alignment_sum_ms += result.total_ms;
          results.push_back(std::move(result));
        }
        else warmup_ms += result.total_ms;
      }
    }
    const auto teardown_start = Clock::now();
    if (workspace) wavefront_aligner_delete(workspace);
    const double teardown_ms = elapsed(teardown_start);
    const double wallcompute_ms = elapsed(wall_start);
    std::cout << std::setprecision(17);
    std::cout << "{\"type\":\"batch\",\"method\":" << json_string(method)
              << ",\"pairs\":" << pairs.size() << ",\"repeats\":" << repeats
              << ",\"warmups\":" << warmups << ",\"batch_setup_ms\":" << setup_ms
              << ",\"batch_teardown_ms\":" << teardown_ms
              << ",\"alignment_sum_ms\":" << alignment_sum_ms
              << ",\"wallcompute_ms\":" << wallcompute_ms
              << ",\"warmup_alignment_ms\":" << warmup_ms << "}\n";
    for (int repeat = 0; repeat < repeats; ++repeat) {
      for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto& pair = pairs[i];
        const auto& result = results[std::size_t(repeat) * pairs.size() + i];
        std::cout << "{\"type\":\"alignment\",\"method\":" << json_string(method)
                  << ",\"id\":" << json_string(pair.id) << ",\"repeat\":" << repeat
                  << ",\"n\":" << pair.x.size() << ",\"m\":" << pair.y.size()
                  << ",\"distance\":" << result.distance
                  << ",\"lower_bound\":" << result.lower_bound
                  << ",\"cigar\":" << json_string(result.cigar)
                  << ",\"total_ms\":" << result.total_ms
                  << ",\"certified_exact\":" << (result.certified_exact ? "true" : "false")
                  << ",\"rounds\":" << result.rounds << ",\"dp_cells\":" << result.dp_cells
                  << ",\"sparse_expansions\":" << result.sparse_expansions
                  << ",\"bound_bytes\":" << result.bound_bytes
                  << ",\"certificate_ms\":" << result.certificate_ms
                  << ",\"bound_preprocess_ms\":" << result.bound_preprocess_ms
                  << ",\"sparse_fallback\":" << json_string(result.sparse_fallback)
                  << ",\"error\":" << json_string(result.error) << ",\"trace\":[";
        for (std::size_t t = 0; t < result.trace.size(); ++t) {
          if (t) std::cout << ',';
          const auto& trace = result.trace[t];
          std::cout << "{\"width\":" << trace.width << ",\"cells\":" << trace.cells
                    << ",\"upper\":" << trace.upper << ",\"lower\":" << trace.lower << '}';
        }
        std::cout << "]}\n";
      }
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "external_adapter: " << error.what() << '\n';
    return 2;
  }
}
