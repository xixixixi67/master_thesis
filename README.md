# master_thesis

Penalized reduced-rank multi-outcome Cox models (survRRR) optimized via an accelerated proximal gradient descent approach (**B = αΓᵀ**) — including the model solvers (R + C++/Eigen), a simulation study, and a real-data application.

## Solvers

Each penalty variant is solved by accelerated proximal gradient descent. They share the same layout: an R entry point (`*_survrrr.R`), a C++/Eigen backend (`*_survrrr.cpp`), and type aliases (`*_types.h`).

- `lasso_gd/` — Lasso penalty on (**α and Γ**). Entry point: `solve_RR_Lasso()`.
- `ridge_gd/` — Ridge penalty on (**α and Γ**). Entry point: `solve_RR_Ridge()`.
- `group_lasso_gd/` — Group-lasso penalty on (**α**). Entry point: `solve_RR_GrpLasso()`.

Other solvers:

- `lasso_cc/survRRR_functions.R` — Lasso reduced-rank method based on criss-cross optimization, built on `glmnet`.
- `mrcox/` — Multi-response Cox (a full-rank method). Entry point: `solve_aligned()`.

## Simulation study

- `simdata_function_and_performance_functions.R` — Functions for the simulation study to compare the estimation performance (simulated-data generators, performance metrics, plotting functions, and so on).
- `simulation_estimated_Beta.Rmd` — Simulation study on coefficient (**B**) estimation, comparing all methods.
- `simulation_prediction.R` — Simulation study on out-of-sample prediction, including the functions.

## Application

- `NEO_analysis.R` — Real-data application to the NEO cohort.
