#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix handle_duplicates_and_ties(NumericMatrix X,
                                             IntegerMatrix nn_idx,
                                             NumericMatrix nn_dists,
                                             int Knn) {
  int n = X.nrow();
  int d = X.ncol();

  // Initialize output matrix with the initial KNN (excluding self)
  IntegerMatrix nn_index_X(n, Knn);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Knn; j++) {
      nn_index_X(i, j) = nn_idx(i, j + 1); // Skip first column (self)
    }
  }

  // Find duplicate data points (distance to nearest neighbor is 0)
  std::vector<int> repeat_data;
  for (int i = 0; i < n; i++) {
    if (nn_dists(i, 1) == 0) {  // nn_dists(i,1) corresponds to nn_dists[i,2] in R
      repeat_data.push_back(i);
    }
  }

  // Handle duplicate points
  if (repeat_data.size() > 0) {
    // Group duplicate points
    std::unordered_map<int, std::vector<int>> groups;
    for (int idx : repeat_data) {
      int group_id = nn_idx(idx, 0) - 1;  // Convert to 0-based
      groups[group_id].push_back(idx);
    }

    // Process each duplicate
    for (int i : repeat_data) {
      int group_id = nn_idx(i, 0) - 1;
      std::vector<int>& group_indices = groups[group_id];

      if (group_indices.size() > Knn) {
        // More duplicates than Knn - sample from group excluding self
        std::vector<int> candidates;
        for (int idx : group_indices) {
          if (idx != i) candidates.push_back(idx + 1);  // Convert to 1-based
        }

        // Random shuffle and take first Knn
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(candidates.begin(), candidates.end(), gen);

        for (int j = 0; j < Knn; j++) {
          nn_index_X(i, j) = candidates[j];
        }
      }
      else {
        // Fewer duplicates than Knn
        if (nn_dists(i, Knn) < nn_dists(i, Knn + 1)) {
          // No tie at Knn - use original result excluding self
          std::vector<int> neighbors;
          for (int j = 0; j <= Knn; j++) {
            if (nn_idx(i, j) - 1 != i) {
              neighbors.push_back(nn_idx(i, j));
            }
          }
          for (int j = 0; j < Knn; j++) {
            nn_index_X(i, j) = neighbors[j];
          }
        }
        else {
          // Handle ties at Knn boundary
          double tie_dist = nn_dists(i, Knn);
          std::vector<int> closer_indices;
          std::vector<int> tied_indices;

          // Calculate distances to all other points
          for (int j = 0; j < n; j++) {
            if (j != i) {
              double dist = 0;
              for (int k = 0; k < d; k++) {
                double diff = X(i, k) - X(j, k);
                dist += diff * diff;
              }
              dist = sqrt(dist);

              if (dist < tie_dist) {
                closer_indices.push_back(j + 1);  // 1-based
              } else if (fabs(dist - tie_dist) < 1e-10) {
                tied_indices.push_back(j + 1);  // 1-based
              }
            }
          }

          // Fill with closer points first
          int pos = 0;
          for (int idx : closer_indices) {
            nn_index_X(i, pos++) = idx;
          }

          // Randomly sample from tied points
          int needed = Knn - closer_indices.size();
          std::random_device rd;
          std::mt19937 gen(rd());
          std::shuffle(tied_indices.begin(), tied_indices.end(), gen);

          for (int j = 0; j < needed; j++) {
            nn_index_X(i, pos++) = tied_indices[j];
          }
        }
      }
    }
  }

  // Find points with ties at Knn boundary (excluding duplicates)
  std::vector<int> ties;
  for (int i = 0; i < n; i++) {
    if (nn_dists(i, Knn) == nn_dists(i, Knn + 1)) {
      // Check if not in repeat_data
      if (std::find(repeat_data.begin(), repeat_data.end(), i) == repeat_data.end()) {
        ties.push_back(i);
      }
    }
  }

  // Handle ties
  for (int i : ties) {
    double tie_dist = nn_dists(i, Knn);
    std::vector<int> closer_indices;
    std::vector<int> tied_indices;

    // Calculate distances to all other points
    for (int j = 0; j < n; j++) {
      if (j != i) {
        double dist = 0;
        for (int k = 0; k < d; k++) {
          double diff = X(i, k) - X(j, k);
          dist += diff * diff;
        }
        dist = sqrt(dist);

        if (dist < tie_dist) {
          closer_indices.push_back(j + 1);  // 1-based
        } else if (fabs(dist - tie_dist) < 1e-10) {
          tied_indices.push_back(j + 1);  // 1-based
        }
      }
    }

    // Fill with closer points first
    int pos = 0;
    for (int idx : closer_indices) {
      if (pos < Knn) nn_index_X(i, pos++) = idx;
    }

    // Randomly sample from tied points
    int needed = Knn - pos;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(tied_indices.begin(), tied_indices.end(), gen);

    for (int j = 0; j < needed && j < tied_indices.size(); j++) {
      nn_index_X(i, pos++) = tied_indices[j];
    }
  }

  return nn_index_X;
}
