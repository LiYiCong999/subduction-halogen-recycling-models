# Input-data requirements

The input data are not included in this template. Before public release, add the shareable source tables or provide permanent repository identifiers and access statements.

## Fig. 4 input

Filename:

    AOC-BS-HP-UHP Halogens data.csv

Required columns:

- Rock_Type
- F_ppm
- Cl_ppm
- either Br_ppm or Br_ppb
- either I_ppm or I_ppb

Rock_Type values are matched exactly and are case-sensitive:

- AOC
- Blueschist
- HP eclogite
- UHP eclogite

All four halogen concentrations must be positive. If both ppm and ppb columns are present for Br or I, the current code prioritizes conversion from the ppb column and overwrites the corresponding ppm column.

## Fig. 5 inputs

Directory expected by the current script:

    code/data source/

### AOC-BS-HP-UHP halogens and H2O_data.csv

Required columns:

- stage
- H2O_min_wt%
- H2O_max_wt%
- F_ppm
- Cl_ppm
- Br_ppb
- I_ppb

Stage names must exactly match those used in stages.csv and mineral_modes.csv.

Important: the water-content values are used directly as fractions in the mass-balance equation. Verify whether 5 wt% is stored as 0.05 or 5. The final repository must document the chosen convention.

### mineral_modes.csv

Required columns:

- stage
- mineral
- min_vol%
- max_vol%

Mineral modes are sampled independently from uniform ranges and then normalized to sum to one.

### partition_coefficients.csv

Required column:

- mineral

For each element, supply either a fixed coefficient such as D_F or a range such as D_F_min and D_F_max. The same naming pattern applies to Cl, Br and I. If the required partition-coefficient columns are absent, the current code silently assigns zero.

### stages.csv

Required columns:

- step
- initial_stage
- residual_stage

The final manuscript and README should explicitly map Steps 1–3 to their geological transitions.

## Fig. 6 input

Filename:

    Forearc serpentinite-Metasediment-BS-HP-UHP halogen_data.csv

Recommended group column:

- Rock_type

Recommended elemental columns:

- F_ppm
- Cl_ppm
- Br_ppb or Br_ppm
- I_ppb or I_ppm

Recommended exact group values:

- Serpentinite
- Metasediment
- Blue schist
- HP eclogite
- UHP eclogite

Use exact standardized group names. Broad keyword aliases such as "Meta" may cause ambiguous substring matching when non-standard labels contain the same text.

## Data provenance

For every compiled or published datum, provide:

- original reference or dataset DOI;
- sample identifier;
- lithological classification;
- original reported unit;
- any conversion applied;
- inclusion or exclusion criteria; and
- redistribution license or access conditions.
