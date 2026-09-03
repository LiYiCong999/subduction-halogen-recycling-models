# Fig. 4 code: usage and interpretation notes

Script:

    code/fig4_aoc_release_contributions.py

## Purpose

This script estimates stage-specific release contributions of F, Cl, Br and I during the AOC–blueschist–HP-eclogite–UHP-eclogite sequence using independent bootstrap resampling of the four lithological groups.

## Required input

Filename:

    AOC-BS-HP-UHP Halogens data.csv

Required column names:

- Rock_Type
- F_ppm
- Cl_ppm
- Br_ppm or Br_ppb
- I_ppm or I_ppb

Required Rock_Type values:

- AOC
- Blueschist
- HP eclogite
- UHP eclogite

The group strings are matched exactly and are case-sensitive.

## Unit handling

- F and Cl must be supplied in ppm.
- If Br_ppb exists, the script creates Br_ppm by dividing by 1,000 and deletes Br_ppb.
- If I_ppb exists, the script creates I_ppm by dividing by 1,000 and deletes I_ppb.
- If both ppb and ppm columns are present, the ppb-derived values overwrite the existing ppm column.
- Rows are retained only when all four halogen concentrations are positive.

## Monte Carlo model

- Number of iterations: 10,000.
- Each lithological group is resampled independently with replacement.
- The number of bootstrap draws equals the original sample number in each group.
- Group means are recalculated in every iteration.
- The three sequential contributions are normalized to the initial AOC inventory.

For a given element:

\[
R_1=100\frac{C_{\mathrm{AOC}}-C_{\mathrm{BS}}}{C_{\mathrm{AOC}}},
\]

\[
R_2=100\left(\frac{C_{\mathrm{BS}}-C_{\mathrm{HP}}}{C_{\mathrm{BS}}}\right)
\left(1-\frac{C_{\mathrm{AOC}}-C_{\mathrm{BS}}}{C_{\mathrm{AOC}}}\right),
\]

and

\[
R_3=100\left(\frac{C_{\mathrm{HP}}-C_{\mathrm{UHP}}}{C_{\mathrm{HP}}}\right)
\left(1-e_1\right)\left(1-e_2\right).
\]

Algebraically, the stage-2 and stage-3 contributions correspond to concentration differences normalized to AOC when the sequential means are internally consistent.

Positive values represent apparent release; negative values represent apparent uptake or enrichment.

## Important reproducibility issue

The current script does not set a random seed. Consequently, repeated runs will not produce identical bootstrap distributions, medians or CSV files.

Before the archival release, choose one of two approaches:

1. add and report a fixed seed; or
2. retain the stochastic implementation and explicitly state that numerical summaries vary slightly among runs.

For Nature Communications, a fixed seed is strongly preferable when the plotted numerical results are expected to be reproduced exactly.

## Robust statistics

- The displayed annotated median is calculated after excluding values outside 1.5 times the interquartile range.
- The helper function also calculates the 2.5th and 97.5th percentiles of the IQR-filtered distribution.
- Those percentile bounds are not plotted in the current figure.
- The printed pandas summary is based on the unfiltered simulation distribution.
- Do not describe the annotated values as means.
- Do not describe the unused percentile bounds as displayed 95% confidence intervals.

## Visualization

- Layout: 2 × 2 panels for F, Cl, Br and I.
- Each stage is represented by a boxplot, normalized horizontal histogram and Gaussian KDE.
- Boxplot outliers are hidden from display but remain in the underlying simulated distribution.
- F, Cl and Br use 100 histogram bins; I uses 200.
- KDEs are evaluated at 200 positions.
- The y-axis is fixed at −400% to 400% for F and −40% to 100% for Cl, Br and I.

Check whether simulated values fall outside these fixed ranges. Clipped values remain in the calculations but will not appear in the displayed panel.

## Output behavior

The script exports:

    F_ppm_stage_contributions.csv
    Cl_ppm_stage_contributions.csv
    Br_ppm_stage_contributions.csv
    I_ppm_stage_contributions.csv

The figure is shown interactively using plt.show() but is not saved automatically. Closing the interactive figure window may be required before CSV export is reached.

Before archiving the code, either:

- save the final Fig. 4 separately and explain that the code displays rather than saves it; or
- add an explicit savefig command and create a new documented software version.

## Failure conditions

- Any missing lithological group produces an empty array and causes bootstrap sampling to fail.
- Constant simulated values can cause gaussian_kde to fail because of a singular covariance matrix.
- Missing or non-numeric columns cause KeyError or conversion-related failure.
- Zero and negative concentrations are removed before group separation.

## Items to report in the manuscript

- Original sample number for each lithological group.
- Data sources and inclusion criteria.
- Exact iteration number and random-seed policy.
- Equations for the three normalized stage contributions.
- Meaning of negative release values.
- Outlier-display and median-annotation rules.
- Histogram bin counts and KDE implementation if the distribution shape supports a main conclusion.
- The fact that simulation counts are not independent geological sample numbers.
