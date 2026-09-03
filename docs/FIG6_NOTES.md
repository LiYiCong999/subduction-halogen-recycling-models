# Fig. 6 code: usage and interpretation notes

Script:

    code/fig6_scenario2_melange_mantle_mixing.py

## Purpose

This script creates three simulated mélange populations using fixed proportions of serpentinite, metasediment and a metamorphic endmember. The alternative metamorphic endmembers are blueschist, HP eclogite and UHP eclogite. Each mélange population is subsequently mixed with DMM.

## Required input

Expected filename:

    Forearc serpentinite-Metasediment-BS-HP-UHP halogen_data.csv

The path is resolved relative to the terminal working directory. Run the script from the code directory or update the path deliberately and document the change.

Recommended columns:

- Rock_type
- F_ppm
- Cl_ppm
- Br_ppb or Br_ppm
- I_ppb or I_ppm

Recommended exact group names:

- Serpentinite
- Metasediment
- Blue schist
- HP eclogite
- UHP eclogite

Use standardized exact names to minimize ambiguous keyword matching.

## Column and group detection

- The group column is detected from a list of common English names unless GROUP_COLUMN is set explicitly.
- Element columns are detected from predefined aliases.
- Br and I are multiplied by 0.001 only when the normalized detected column name ends in ppb.
- F and Cl are assumed to be in ppm.
- Rows with missing, infinite, zero or negative values for any halogen are removed.

The alias "Meta" is broad. If exact matching fails, substring matching may incorrectly classify a non-metasedimentary group containing the same character sequence. Inspect the printed group values and use exact standardized labels.

## Three-endmember mélange model

Each simulation uses fixed proportions:

- serpentinite = 0.50;
- metasediment = 0.05; and
- metamorphic endmember = 0.45.

For element \(i\):

\[
C_{i,\mathrm{m\acute{e}lange}}
=0.50C_{i,\mathrm{serpentinite}}
+0.05C_{i,\mathrm{metasediment}}
+0.45C_{i,\mathrm{metamorphic}}.
\]

For each of 100,000 realizations, one measured row is independently sampled with replacement from each source group. The three metamorphic alternatives generate three separate mélange populations.

Random seed:

    42

The model propagates within-group concentration variability but does not propagate uncertainty in the 50:5:45 proportions. If those proportions are interpretive rather than directly constrained, provide a geological justification and sensitivity analysis.

Serpentinite and metasediment are resampled separately for each of the three mélange populations. The three scenario outputs are therefore not paired realization by realization.

## DMM mixing

Fixed DMM concentrations:

- F = 12.0 ppm
- Cl = 5.0 ppm
- Br = 0.013 ppm
- I = 0.0003 ppm

Mixing is linear in concentration:

\[
C_{i,\mathrm{mix}}=(1-r)C_{i,\mathrm{m\acute{e}lange}}+rC_{i,\mathrm{DMM}}.
\]

The code evaluates:

- a fixed DMM fraction of 0.50 for the density figure; and
- DMM fractions from 0.10 to 0.90 at intervals of 0.10 for the ratio figure.

The model does not include reaction, partitioning, precipitation, selective fluid loss or isotope fractionation during mélange–mantle interaction.

## Density visualization

- Layout: 3 × 2.
- Rows: blueschist-, HP-eclogite- and UHP-eclogite-bearing mélange.
- Columns: F/Cl versus I/Cl and Br/Cl versus I/Cl.
- Bivariate KDE is fitted in log10 ratio space.
- The grid contains 100 × 100 positions and extends 5% beyond the modeled log-space range.
- Filled levels are 20%, 40%, 60%, 80% and 100% of maximum density.
- The levels are relative peak densities, not probability-mass contours.
- Marginal histograms use 20 logarithmic bins.
- Marginal univariate KDEs are evaluated at 240 positions and scaled to histogram counts.

Fixed plotting ranges:

- I/Cl: 10⁻⁶ to 10⁻¹;
- F/Cl: 10⁻² to 10¹; and
- Br/Cl: 10⁻³ to 10⁻¹.

Simulated values outside these ranges are calculated but not visible in the plotted panels. Quantify the fraction outside each range before publication.

## Mantle-fraction visualization

- Layout: 1 × 3 for F/Cl, Br/Cl and I/Cl.
- DMM fractions: 10%–90% in 10% increments.
- Values outside the 2.5th–97.5th percentile interval are removed.
- Points represent medians of the retained distributions.
- Error bars represent ±1 sample standard deviation of the retained values.
- The three scenarios are offset horizontally to avoid overlap.

Do not describe the error bars as 95% confidence intervals. Do not use the number of simulated realizations as the geological sample size.

## Expected outputs

    1_melange_and_mantle_mixing_three_by_two.svg
    halogen_ratios_melange_50-45-5_and_mantle_mixing_10-90.svg
    mantle_mixing_10-90_central95_median_1sigma.csv

The figures are exported as editable SVG files, with text retained as text objects.

## Failure conditions

- Missing source groups cause group-extraction failure.
- Fewer than five requested simulations are rejected.
- Mixture proportions must sum to one and cannot be negative.
- Fewer than five valid points prevent bivariate KDE calculation.
- Constant or collinear log-ratio data may produce a singular KDE covariance matrix; the function returns no density field.

## Items to report in the manuscript

- Original sample number and source for every endmember group.
- Exact group-classification rules.
- Fixed 50:5:45 proportions and their justification.
- Iteration number and random seed.
- Whether endmember rows are sampled independently.
- DMM composition and reference.
- Binary concentration-mixing equation.
- Processes excluded from the mixing model.
- KDE space, grid resolution and relative-density definition.
- Fixed plot limits and any excluded visual range.
- Central-95% trimming and meaning of ±1σ.
