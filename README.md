# HPV vaccination model with risk heterogeneity

MATLAB code and processed data for the HPV vaccination study in South Korea.
Tested with MATLAB R2025b. Refitting requires Optimization Toolbox; PRCC and
bootstrapping require Statistics and Machine Learning Toolbox.

Open this folder in MATLAB and run `setup_paths` once. Then run the analyses needed:

```matlab
setup_paths
DataFitting                         % refit the baseline model
DataFittingBootstrap                % plot using saved bootstrap estimates
SensitivityAnalysis('R')             % saved results: 'R', 'Cf', 'Cm', or 'Cmsm'
HIT_Contour                         % baseline vaccination coverage contours
HIT_Contour(1.35)                    % lower within-MSM contact scenario
AvertedCancers                      % vaccination and spillover simulations
```

To generate new bootstrap estimates, use `DataFittingBootstrap(true,1000)`.
To repeat sensitivity sampling and screening, use
`SensitivityAnalysis('Cf',true)` (replace `'Cf'` with the desired response).
For the additional PRCC runs, use `SensitivityAnalysis('all',false,seed)` with
seeds 2, 3 and 4; the saved OAT selections and ranges are retained.

Inputs are in `data/` and saved estimates and samples are in `results/`:
- `DataFittingVariables.mat`: baseline fit, initial state and fit residuals.
- `DataFittingBootstrapResult.mat`: 1,000 fitted parameter pairs and their intervals.
- `Sensitivity_*.mat`: LHS samples, responses, OAT curves and PRCCs for each outcome.
- `PRCC_stability.csv`: additional seeds, PRCCs, p-values and ranks (seeds 2--4).

New figures and results are written to `outputs/`. Each analysis uses the
bundled baseline in `results/`; a new fit in `outputs/` is not used automatically.

Data sources are listed in `data/README.md`. The adapted LHS/PRCC code is
acknowledged in `THIRD_PARTY_NOTICES.md`.
