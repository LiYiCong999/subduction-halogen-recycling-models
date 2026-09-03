# Fig. 5 code: usage and interpretation notes

Script:

    code/fig5_scenario1_fluid_mantle_mixing.py

## Purpose

This script simulates slab-derived fluid halogen concentrations from bulk-rock compositions, water loss, residual mineral modes and mineral–fluid partition coefficients. The modeled fluid is then mixed with depleted MORB mantle (DMM).

## Required directory

The script defines:

    BASE_DIR = Path("data source")

The path is resolved relative to the terminal working directory, not automatically relative to the script. The recommended procedure is to enter the code directory before running and create:

    code/data source/

## Required input files

### AOC-BS-HP-UHP halogens and H2O_data.csv

Required columns:

- stage
- H2O_min_wt%
- H2O_max_wt%
- F_ppm
- Cl_ppm
- Br_ppb
- I_ppb

### mineral_modes.csv

Required columns:

- stage
- mineral
- min_vol%
- max_vol%

### partition_coefficients.csv

Required column:

- mineral

For each mineral and element, use either:

- D_F_min and D_F_max, or D_F;
- D_Cl_min and D_Cl_max, or D_Cl;
- D_Br_min and D_Br_max, or D_Br;
- D_I_min and D_I_max, or D_I.

### stages.csv

Required columns:

- step
- initial_stage
- residual_stage

The stage labels must match exactly across all input tables.

## Monte Carlo settings

- Number of iterations: 100,000.
- Random seed: 42.
- The initial bulk-rock composition is sampled from one observed row for each initial stage.
- Initial and residual water contents are sampled independently from uniform minimum–maximum ranges.
- Each mineral abundance is sampled independently from a uniform range and then all sampled modes are normalized to sum to one.
- Each partition coefficient is sampled independently from a uniform minimum–maximum range when both limits exist.

These independent uniform distributions are model assumptions. Correlations among water content, mineral modes, partition coefficients and bulk compositions are not represented.

## Critical water-content convention

The code calculates:

\[
f_{\mathrm{fluid}}=W_{\mathrm{initial}}-W_{\mathrm{residual}}
\]

and uses this value directly in:

\[
C_{i,\mathrm{fluid}}=
\frac{C_{i,0}}
{f_{\mathrm{fluid}}+(1-f_{\mathrm{fluid}})\overline{D}_i}.
\]

Although the columns are labelled wt%, the script performs no division by 100. Confirm whether 5 wt% is stored as 0.05 or 5. If whole percentage numbers are stored, the equation will not use a valid fractional fluid mass.

This unit convention must be checked before publication and stated explicitly in the README and Methods.

## Partition coefficients

The bulk partition coefficient is:

\[
\overline{D}_i=\sum_j X_jD_{i,j}.
\]

Important behaviors:

- If a mineral is absent from partition_coefficients.csv, every halogen partition coefficient for that mineral is set to zero.
- If the element-specific columns are absent, that coefficient is also set to zero.
- The code does not issue a warning when zero is assigned because of missing information.

Before release, audit the input table so that zero means an intentional scientific assumption rather than an omitted value.

## Fluid-composition exclusions

The calculation returns a missing value when:

- the initial concentration is missing;
- the released-fluid fraction is no greater than 1 × 10⁻⁹; or
- the mass-balance denominator is non-positive.

The number of retained simulations may therefore differ among elements and steps. Report retained simulation counts.

## DMM mixing

Fixed DMM concentrations:

- F = 12.0 ppm
- Cl = 5.0 ppm
- Br = 0.013 ppm
- I = 0.0003 ppm

The mixing equation is:

\[
C_{i,\mathrm{mix}}=(1-r)C_{i,\mathrm{fluid}}+rC_{i,\mathrm{DMM}}.
\]

The code assumes concentration-linear binary mixing. It does not model mineral reaction, fluid consumption, precipitation, kinetic fractionation or element-specific partitioning during mantle interaction.

The DMM composition is treated as fixed and carries no uncertainty. Cite and justify the source of all four DMM values.

## Visualization

### Fixed 50% mantle mixing

- The code plots F/Cl versus I/Cl and Br/Cl versus I/Cl for Steps 1–3.
- Ratios are transformed to log10 space before bivariate Gaussian KDE.
- KDEs are evaluated on a 100 × 100 grid.
- Filled levels correspond to 20%, 40%, 60%, 80% and 100% of the maximum density.
- These are not contours enclosing the same percentages of probability mass.
- Marginal histograms use 20 logarithmic bins.
- Marginal KDEs are fitted in log10 space and evaluated at 200 positions.

### Mantle fractions from 10% to 90%

- Fractions are evaluated at 10% intervals.
- Values below the 2.5th percentile and above the 97.5th percentile are removed.
- The median of the retained central 95% is plotted.
- Error bars show ±1 sample standard deviation of the retained distribution.
- The error bars are not 95% confidence intervals.

## Expected outputs

    fluid_concentration_all_iterations.csv
    fluid_concentration_summary.csv
    density_50mantle_steps.png
    halogen_ratios_mantle_mixing_10_90.png

The figures are saved at 300 dpi as PNG files. For publication-quality line graphics, consider whether an SVG or PDF export is needed; if the code is changed, release a new version.

The full-iteration CSV may be large. It can be regenerated from the fixed seed and inputs. If it exceeds practical GitHub limits, store it in Zenodo or another data repository rather than committing it repeatedly.

## Additional implementation notes

- RUN_MONTE_CARLO=True reruns the complete model.
- RUN_MONTE_CARLO=False reads fluid_concentration_all_iterations.csv.
- SHOW_FIGURE=True opens interactive figures.
- The validation error messages refer to bulk_data.csv even though the configured bulk filename is different; this affects messages, not calculations.
- The results depend on the order and exact spelling of stage and mineral names.

## Items to report in the manuscript

- Geological meaning of Steps 1–3.
- Sources of water-content, mineral-mode and partition-coefficient ranges.
- Distribution assigned to each uncertain parameter.
- Number of observations in every initial-stage population.
- Number of Monte Carlo iterations and seed.
- Water-content unit convention.
- Missing-coefficient equals zero assumption.
- Fluid mass-balance equation.
- DMM composition and reference.
- Mixing equation and excluded processes.
- KDE construction and density-level definition.
- Central-95% trimming and exact meaning of ±1σ.
