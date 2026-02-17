#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>

using namespace Rcpp;

// ---- helpers: hash a row by exact double bit patterns ----
static inline uint64_t dbl_bits(double x) {
  uint64_t u;
  std::memcpy(&u, &x, sizeof(double));
  return u;
}

struct RowKey {
  std::vector<uint64_t> bits;
  bool operator==(RowKey const& other) const {
    return bits == other.bits;
  }
};

struct RowKeyHash {
  std::size_t operator()(RowKey const& k) const {
    // simple hash combine
    std::size_t h = 1469598103934665603ULL;
    for (auto b : k.bits) {
      h ^= (std::size_t)b + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    }
    return h;
  }
};

static RowKey make_key(const NumericMatrix& X, int i) {
  int d = X.ncol();
  RowKey k;
  k.bits.reserve(d);
  for (int col = 0; col < d; ++col) k.bits.push_back(dbl_bits(X(i, col)));
  return k;
}

// Fisher–Yates shuffle using R RNG (respects set.seed)
static void shuffle_inplace(std::vector<int>& a) {
  for (int i = (int)a.size() - 1; i > 0; --i) {
    int j = (int)std::floor(R::runif(0.0, 1.0) * (i + 1));
    std::swap(a[i], a[j]);
  }
}

// squared Euclidean distance between rows i and j
static inline double dist2_row(const NumericMatrix& X, int i, int j) {
  int d = X.ncol();
  double s = 0.0;
  for (int k = 0; k < d; ++k) {
    double diff = X(i, k) - X(j, k);
    s += diff * diff;
  }
  return s;
}

// [[Rcpp::export]]
IntegerMatrix knn_tie_duplicate_aware_cpp(NumericMatrix X,
                                          IntegerMatrix nn_idx,
                                          NumericMatrix nn_dists,
                                          int Knn,
                                          double tol = 1e-12) {
  RNGScope scope; // use R RNG

  const int n = X.nrow();
  if (Knn < 1) stop("Knn must be >= 1");
  if (Knn > n - 1) Knn = n - 1;

  // ---- build exact-duplicate groups by hashing rows ----
  std::unordered_map<RowKey, std::vector<int>, RowKeyHash> groups;
  groups.reserve((size_t)n);

  for (int i = 0; i < n; ++i) {
    groups[make_key(X, i)].push_back(i); // store 0-based indices
  }

  IntegerMatrix out(n, Knn);

  // nn_idx and nn_dists come from RANN::nn2; nn_idx is 1-based indices in R
  // We'll treat them as given and convert to 0-based when needed.

  for (int i = 0; i < n; ++i) {
    std::vector<int> chosen;
    chosen.reserve(Knn);

    // ---- Step A: handle duplicates at distance 0 exactly ----
    auto &dup_group = groups[make_key(X, i)];
    std::vector<int> dup_candidates;
    dup_candidates.reserve(dup_group.size());

    for (int id : dup_group) {
      if (id != i) dup_candidates.push_back(id);
    }

    if (!dup_candidates.empty()) {
      shuffle_inplace(dup_candidates);
      int take = std::min((int)dup_candidates.size(), Knn);
      for (int t = 0; t < take; ++t) chosen.push_back(dup_candidates[t]);
    }

    if ((int)chosen.size() == Knn) {
      for (int j = 0; j < Knn; ++j) out(i, j) = chosen[j] + 1; // 1-based
      continue;
    }

    // ---- Step B: fill remaining using RANN ordering unless boundary tie forces full enumeration ----
    int need = Knn - (int)chosen.size();

    // Build a small candidate list from RANN (skipping self and already chosen)
    // Note: nn_idx(i,0) is usually self.
    int k_cols = nn_idx.ncol();

    std::vector<int> rann_list;
    rann_list.reserve(k_cols);

    // mark chosen for fast skip
    std::vector<char> is_chosen(n, 0);
    for (int id : chosen) is_chosen[id] = 1;

    for (int col = 0; col < k_cols; ++col) {
      int id0 = nn_idx(i, col) - 1; // to 0-based
      if (id0 < 0 || id0 >= n) continue;
      if (id0 == i) continue;
      if (is_chosen[id0]) continue;
      rann_list.push_back(id0);
    }

    // If RANN list is too short (can happen if k_cols small), we'll still handle via full scan below.
    // Propose a cutoff distance r based on RANN if possible.
    bool do_full = false;
    double r = NA_REAL;

    if ((int)rann_list.size() >= need) {
      // r is the distance of the "need-th" element in the RANN-proposed list.
      // We need to locate its column distance in nn_dists.
      // Easiest: read nn_dists row and walk columns collecting usable candidates with their dists.
      std::vector<double> cand_d;
      cand_d.reserve(k_cols);

      for (int col = 0; col < k_cols; ++col) {
        int id0 = nn_idx(i, col) - 1;
        if (id0 < 0 || id0 >= n) continue;
        if (id0 == i) continue;
        if (is_chosen[id0]) continue;
        cand_d.push_back(nn_dists(i, col));
      }

      // distance at boundary among remaining picks:
      r = cand_d[need - 1];

      // tie detect: if we can see the next candidate distance and it's equal (within tol), tie exists.
      if ((int)cand_d.size() > need) {
        if (std::fabs(cand_d[need] - r) <= tol) do_full = true;
      } else {
        // no lookahead info -> be safe if r is 0 or extremely small; otherwise assume no tie
        // (ties are usually caught with lookahead given you call nn2 with Knn+2)
      }
    } else {
      do_full = true; // insufficient candidates -> full enumeration
    }

    if (!do_full) {
      // take first 'need' from RANN order
      for (int t = 0; t < need; ++t) {
        chosen.push_back(rann_list[t]);
      }
      for (int j = 0; j < Knn; ++j) out(i, j) = chosen[j] + 1;
      continue;
    }

    // ---- Step C: full enumeration to sample uniformly from the full tie set at the cutoff radius ----
    // If r is NA (e.g., not enough RANN candidates), we can define r as the smallest radius that yields enough points:
    // simplest: compute all distances, sort, and pick r as the need-th smallest among remaining points.
    // We'll do: compute dist2 to all candidates not in chosen and not self, then determine r2 accordingly.

    std::vector<std::pair<double,int>> all; // (dist2, id)
    all.reserve(n - 1);

    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      if (is_chosen[j]) continue;
      double d2 = dist2_row(X, i, j);
      all.push_back({d2, j});
    }

    std::sort(all.begin(), all.end(),
              [](auto const& a, auto const& b){ return a.first < b.first; });

    if ((int)all.size() < need) {
      stop("Not enough points to fill neighbor set (unexpected).");
    }

    double r2 = all[need - 1].first;

    // collect strictly closer and tied-at-boundary
    std::vector<int> closer;
    std::vector<int> tied;
    closer.reserve(need);
    tied.reserve(need);

    for (auto const& pr : all) {
      double d2 = pr.first;
      int id = pr.second;
      if (d2 < r2 - tol) closer.push_back(id);
      else if (std::fabs(d2 - r2) <= tol) tied.push_back(id);
      else break; // sorted
    }

    // fill from closer first
    for (int id : closer) {
      if ((int)chosen.size() < Knn) chosen.push_back(id);
    }

    int remaining = Knn - (int)chosen.size();
    if (remaining > 0) {
      shuffle_inplace(tied);
      if ((int)tied.size() < remaining) {
        // should not happen, but guard:
        remaining = (int)tied.size();
      }
      for (int t = 0; t < remaining; ++t) chosen.push_back(tied[t]);
    }

    if ((int)chosen.size() != Knn) stop("Failed to construct Knn neighbors.");

    for (int j = 0; j < Knn; ++j) out(i, j) = chosen[j] + 1; // 1-based
  }

  return out;
}
