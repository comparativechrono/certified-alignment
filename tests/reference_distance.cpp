// Independent unrestricted unit-cost Levenshtein oracle.
// Deliberately includes no alignment-engine code, bounds, or predictors.
// Protocol: id<TAB>x<TAB>y<NEWLINE>; sequences are literal bytes other than
// TAB and NEWLINE. Output is one JSON object per input record.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

static std::string json_string(const std::string& value) {
    static constexpr char hex[] = "0123456789abcdef";
    std::string result = "\"";
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') {
            result += '\\';
            result += static_cast<char>(c);
        } else if (c < 32 || c >= 127) {
            result += "\\u00";
            result += hex[c >> 4];
            result += hex[c & 15];
        } else {
            result += static_cast<char>(c);
        }
    }
    return result + '"';
}

static std::uint32_t distance(const std::string& x, const std::string& y) {
    // Keep the shorter input along columns; this only reduces oracle memory.
    const std::string& rows = x.size() >= y.size() ? x : y;
    const std::string& columns = x.size() >= y.size() ? y : x;
    if (rows.size() >= std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("input too long for 32-bit oracle scores");
    std::vector<std::uint32_t> previous(columns.size() + 1);
    std::vector<std::uint32_t> current(columns.size() + 1);
    for (std::size_t j = 0; j <= columns.size(); ++j)
        previous[j] = static_cast<std::uint32_t>(j);
    for (std::size_t i = 1; i <= rows.size(); ++i) {
        current[0] = static_cast<std::uint32_t>(i);
        for (std::size_t j = 1; j <= columns.size(); ++j) {
            const std::uint32_t deletion = previous[j] + 1;
            const std::uint32_t insertion = current[j - 1] + 1;
            const std::uint32_t substitution = previous[j - 1]
                + static_cast<std::uint32_t>(rows[i - 1] != columns[j - 1]);
            current[j] = std::min({deletion, insertion, substitution});
        }
        previous.swap(current);
    }
    return previous.back();
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    try {
        std::string line;
        std::size_t record = 0;
        while (std::getline(std::cin, line)) {
            ++record;
            const auto first = line.find('\t');
            const auto second = first == std::string::npos
                ? std::string::npos : line.find('\t', first + 1);
            if (first == std::string::npos || second == std::string::npos
                || line.find('\t', second + 1) != std::string::npos)
                throw std::runtime_error("malformed TSV at record "
                                         + std::to_string(record));
            const auto id = line.substr(0, first);
            const auto x = line.substr(first + 1, second - first - 1);
            const auto y = line.substr(second + 1);
            std::cout << "{\"id\":" << json_string(id)
                      << ",\"distance\":" << distance(x, y) << "}\n";
        }
        if (!std::cin.eof()) throw std::runtime_error("input read failure");
    } catch (const std::exception& error) {
        std::cerr << "reference_distance: " << error.what() << '\n';
        return 2;
    }
}
