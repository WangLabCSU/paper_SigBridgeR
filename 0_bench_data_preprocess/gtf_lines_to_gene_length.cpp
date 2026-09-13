#include <Rcpp.h>
#include <zlib.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]

struct Interval {
  int start;
  int end;
};

using GroupMap = unordered_map<string, vector<Interval>>;

// ============================================================
// Line reader
// ============================================================

class LineReader {
private:
  bool is_gz_;
  gzFile gz_fp_;
  ifstream file_;
  vector<char> buffer_;

public:
  explicit LineReader(const string &path)
      : is_gz_(false), gz_fp_(nullptr), buffer_(1 << 16) {

    const bool gz =
        path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;

    if (gz) {
      is_gz_ = true;
      gz_fp_ = gzopen(path.c_str(), "rb");

      if (gz_fp_ == nullptr) {
        stop("Cannot open gzip file: " + path);
      }
    } else {
      file_.open(path);

      if (!file_) {
        stop("Cannot open file: " + path);
      }
    }
  }

  ~LineReader() {
    if (gz_fp_ != nullptr) {
      gzclose(gz_fp_);
      gz_fp_ = nullptr;
    }
  }

  bool next(string &line) {
    line.clear();

    if (is_gz_) {
      char tmp[1 << 16];

      if (gzgets(gz_fp_, tmp, sizeof(tmp)) == nullptr) {
        return false;
      }

      line.assign(tmp);

      // 处理超过缓冲区长度的长行
      while (!line.empty() && line.back() != '\n' && !gzeof(gz_fp_)) {
        if (gzgets(gz_fp_, tmp, sizeof(tmp)) == nullptr) {
          break;
        }
        line += tmp;
      }
    } else {
      if (!std::getline(file_, line)) {
        return false;
      }
    }

    if (!line.empty() && line.back() == '\n') {
      line.pop_back();
    }

    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }

    return true;
  }
};

// ============================================================
// Helpers
// ============================================================

static bool is_canonical_chr(const string &chr) {
  string x = chr;

  if (x.rfind("chr", 0) == 0) {
    x = x.substr(3);
  }

  if (x == "X" || x == "Y" || x == "M" || x == "MT") {
    return true;
  }

  if (x.empty()) {
    return false;
  }

  for (char c : x) {
    if (c < '0' || c > '9') {
      return false;
    }
  }

  int n = atoi(x.c_str());
  return n >= 1 && n <= 22;
}

static bool get_column(const char *line, size_t n, int target, string &out) {
  int column = 0;
  size_t begin = 0;

  for (size_t i = 0; i <= n; ++i) {
    if (i == n || line[i] == '\t') {
      if (column == target) {
        out.assign(line + begin, i - begin);
        return true;
      }

      ++column;
      begin = i + 1;
    }
  }

  return false;
}

static bool parse_gene_id(const string &attributes, string &gene_id) {
  const string key = "gene_id";
  size_t pos = attributes.find(key);

  while (pos != string::npos) {
    bool valid_left = pos == 0 || attributes[pos - 1] == ' ' ||
                      attributes[pos - 1] == ';' || attributes[pos - 1] == '\t';

    size_t p = pos + key.size();

    while (p < attributes.size() && attributes[p] == ' ') {
      ++p;
    }

    if (valid_left && p < attributes.size() && attributes[p] == '"') {
      ++p;
      size_t q = attributes.find('"', p);

      if (q != string::npos) {
        gene_id = attributes.substr(p, q - p);
        return !gene_id.empty();
      }
    }

    pos = attributes.find(key, pos + 1);
  }

  return false;
}

// ============================================================
// Add one GTF line
// ============================================================

static void add_gtf_line_fast(const char *line, size_t n, GroupMap &groups,
                              bool canonical_chr_only) {
  if (n == 0 || line[0] == '#') {
    return;
  }

  string seqname;
  string feature;
  string start_str;
  string end_str;
  string strand;
  string attributes;

  if (!get_column(line, n, 0, seqname) || !get_column(line, n, 2, feature) ||
      !get_column(line, n, 3, start_str) || !get_column(line, n, 4, end_str) ||
      !get_column(line, n, 6, strand) || !get_column(line, n, 8, attributes)) {
    return;
  }

  if (feature != "exon") {
    return;
  }

  if (canonical_chr_only && !is_canonical_chr(seqname)) {
    return;
  }

  string gene_id;
  if (!parse_gene_id(attributes, gene_id)) {
    return;
  }

  int start;
  int end;

  try {
    start = stoi(start_str);
    end = stoi(end_str);
  } catch (...) {
    return;
  }

  if (start > end) {
    swap(start, end);
  }

  // gene_id + chromosome + strand 分组
  string key = gene_id + '\t' + seqname + '\t' + strand;
  groups[key].push_back({start, end});
}

// ============================================================
// Read file
// ============================================================

static void read_gtf_lines(const string &path, GroupMap &groups,
                           bool canonical_chr_only) {
  LineReader reader(path);
  string line;
  size_t count = 0;

  while (reader.next(line)) {
    add_gtf_line_fast(line.data(), line.size(), groups, canonical_chr_only);

    if ((count++ & 131071) == 0) {
      checkUserInterrupt();
    }
  }
}

// ============================================================
// Merge intervals and calculate gene lengths
// ============================================================

static NumericVector finalize_gene_lengths(GroupMap &groups) {
  unordered_map<string, long long> gene_lengths;
  gene_lengths.reserve(groups.size());

  for (auto &item : groups) {
    vector<Interval> &intervals = item.second;

    if (intervals.empty()) {
      continue;
    }

    sort(intervals.begin(), intervals.end(),
         [](const Interval &a, const Interval &b) {
           if (a.start != b.start) {
             return a.start < b.start;
           }
           return a.end < b.end;
         });

    int current_start = intervals[0].start;
    int current_end = intervals[0].end;
    long long total = 0;

    for (size_t i = 1; i < intervals.size(); ++i) {
      const Interval &current = intervals[i];

      if (current.start <= current_end + 1) {
        if (current.end > current_end) {
          current_end = current.end;
        }
      } else {
        total += static_cast<long long>(current_end - current_start + 1);

        current_start = current.start;
        current_end = current.end;
      }
    }

    total += static_cast<long long>(current_end - current_start + 1);

    const size_t tab = item.first.find('\t');

    if (tab != string::npos) {
      const string gene_id = item.first.substr(0, tab);
      gene_lengths[gene_id] += total;
    }
  }

  NumericVector result(gene_lengths.size());
  CharacterVector names(gene_lengths.size());

  R_xlen_t i = 0;

  for (const auto &item : gene_lengths) {
    names[i] = item.first;
    result[i] = static_cast<double>(item.second);
    ++i;
  }

  result.attr("names") = names;
  return result;
}

// ============================================================
// R character vector
// ============================================================

// [[Rcpp::export]]
NumericVector gene_length_from_r_lines(const CharacterVector &lines,
                                       bool canonical_chr_only = false) {

  GroupMap groups;
  groups.reserve(65536);

  for (R_xlen_t i = 0; i < lines.size(); ++i) {
    SEXP element = STRING_ELT(lines, i);

    if (element == NA_STRING) {
      continue;
    }

    const char *line = CHAR(element);

    add_gtf_line_fast(line, strlen(line), groups, canonical_chr_only);

    if ((i & 131071) == 0) {
      checkUserInterrupt();
    }
  }

  return finalize_gene_lengths(groups);
}

// ============================================================
// R file interface
// ============================================================

// [[Rcpp::export]]
NumericVector gtf_file_to_gene_length(std::string path,
                                      bool canonical_chr_only = false) {

  GroupMap groups;
  groups.reserve(65536);

  read_gtf_lines(path, groups, canonical_chr_only);

  return finalize_gene_lengths(groups);
}