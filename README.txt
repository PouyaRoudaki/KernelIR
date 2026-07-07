KernelIR
========
R package implementing the KNN-based and RKHS-based estimators of Kernel Integrated R^2.

- KNN-based estimator:  R/kir_graph.R
- RKHS-based estimator: R/kir_rkhs.R

Contents of inst/
------------------

7.1.Power_Comparision -- Power comparison study (Section 7.1)
  - power_heteroscedastic.R : Euclidean example with a heteroscedastic alternative hypothesis
  - power_SO3.R             : Non-Euclidean example with an SO(3) alternative hypothesis
  - plot_Figure1.R          : Plotting script for Figure 1

7.1.Runtime_Comparison -- Runtime comparison study (Section 7.1)
  - runtime_comparison.R

7.2.Million_Song_Data_Study -- Power comparison study for the Million Song Dataset (Section 7.2)
  - Dataset: https://archive.ics.uci.edu/dataset/203/yearpredictionmsd
  - Million_Song_Data_study.R : R implementation of the power study
  - plot_Figure2.R            : Plotting script for Figure 2

7.3.Convergence_Rate -- Intrinsic-dimension effect in Theorem 2 (Appendix D additional experiment)
  - convergence_rate.R                  : Example of independent Haar-distributed orthogonal matrices (convergence rate)
  - tikz_convergence_rate_exp2_subset.R : Plotting script for Figure 3
