# KernelIR

R package implementing **kernel integrated R²** (`D(Y, X)`), a nonparametric measure of statistical dependence that extends [integrated R²](https://arxiv.org/abs/2505.18146) to responses `Y` taking values in general spaces equipped with a characteristic kernel (multivariate, functional, or structured data such as SO(3)). Full details are in the paper: [*Kernel Integrated R²: A Measure of Dependence*](https://arxiv.org/pdf/2602.22985) (arXiv:2602.22985).

`D(Y, X) ∈ [0, 1]`, with `D(Y, X) = 0` iff `X` and `Y` are independent, and `D(Y, X) = 1` iff `Y` is almost surely a measurable function of `X`. Unlike symmetric measures (HSIC, distance covariance), `D` is directional, like [Chatterjee's ξ](https://doi.org/10.1080/01621459.2021.1934784) and [Azadkia–Chatterjee's ν](https://arxiv.org/abs/2505.18146).

Two estimators are provided:

- **`kir_graph()`** — K-nearest-neighbour graph-based estimator. Only requires `X` to lie in a metric space. Consistent, with a convergence rate that adapts to the intrinsic dimension of `X` (see Theorem 2 in the paper).
- **`kir_rkhs()`** — RKHS-based estimator built on conditional mean embeddings. Requires a kernel on `X` rather than a metric.

## Installation

```r
# install.packages("devtools")
devtools::install_github("PouyaRoudaki/KernelIR")
```

The package links against Rcpp/RcppArmadillo, so a working C++ toolchain (Rtools on Windows) is required to build from source.

## Usage

```r
library(KernelIR)

set.seed(1)
n <- 200
X <- matrix(rnorm(n * 2), n, 2)

# Strong dependence: Y is a noisy function of X
Y_dep <- X[, 1, drop = FALSE] + 0.3 * rnorm(n)
kir_graph(X, Y_dep, Knn = 5)   # close to 1
kir_rkhs(X, Y_dep)             # close to 1

# Near independence
Y_ind <- matrix(rnorm(n), n, 1)
kir_graph(X, Y_ind, Knn = 5)   # close to 0
kir_rkhs(X, Y_ind)             # close to 0
```

Both estimators accept a custom kernel on `Y` (and, for `kir_rkhs()`, on `X`) via `kernlab` kernel objects or a `function(A, B = NULL)` returning a Gram matrix; see `?kir_graph` and `?kir_rkhs` for full argument details.

## Repository structure

```
R/            Estimator implementations (kir_graph.R, kir_rkhs.R, get_knn.R, helper.R)
src/          Rcpp/RcppArmadillo backends for the KNN graph and centering routines
man/          Function documentation
inst/         Scripts reproducing the paper's experiments (see below)
```

### Reproducing the paper's experiments (`inst/`)

| Folder | Contents |
|---|---|
| `7.1.Power_Comparision` | Power comparison study (Section 7.1): heteroscedastic Euclidean example (`power_heteroscedastic.R`), non-Euclidean SO(3) example (`power_SO3.R`), and the Figure 1 plotting script (`plot_Figure1.R`) |
| `7.1.Runtime_Comparison` | Runtime comparison study (Section 7.1): `runtime_comparison.R` |
| `7.2.Million_Song_Data_Study` | Power comparison on the [Million Song Dataset](https://archive.ics.uci.edu/dataset/203/yearpredictionmsd) (Section 7.2): `Million_Song_Data_study.R` plus the Figure 2 plotting script (`plot_Figure2.R`) |
| `7.3.Convergence_Rate` | Intrinsic-dimension effect in Theorem 2 (Appendix D): the Haar-distributed orthogonal matrices example (`convergence_rate.R`) and the Figure 3 plotting script (`tikz_convergence_rate_exp2_subset.R`) |

## Citation

If you use this package, please cite:

```bibtex
@InProceedings{roudaki2026kernelir,
  title     = {Kernel Integrated $R^2$: A Measure of Dependence},
  author    = {Roudaki, Pouya and Gavioli-Akilagun, Shakeel and Kalinke, Florian and Azadkia, Mona and Szab\'o, Zolt\'an},
  booktitle = {Proceedings of the 42nd Conference on Uncertainty in Artificial Intelligence (UAI)},
  year      = {2026}
}
```

## License

MIT — see [LICENSE](LICENSE).
