# Halogen recycling models for subduction-zone metamorphism

This repository contains the Python code used to generate Figs. 4–6 of the associated Nature Communications manuscript. The scripts quantify stage-specific halogen release from altered oceanic crust (AOC), simulate fluid compositions during progressive metamorphic dehydration, and evaluate two alternative mechanisms of mantle modification.

## Repository contents

| Manuscript figure | Script | Purpose |
|---|---|---|
| Fig. 4 | code/fig4_aoc_release_contributions.py | Bootstrap Monte Carlo estimates of F, Cl, Br and I release contributions during the AOC–blueschist–HP-eclogite–UHP-eclogite sequence |
| Fig. 5 | code/fig5_scenario1_fluid_mantle_mixing.py | Mineral-mode and partition-coefficient constrained simulation of slab-derived fluid compositions, followed by fluid–DMM mixing |
| Fig. 6 | code/fig6_scenario2_melange_mantle_mixing.py | Three-endmember mélange construction followed by mélange–DMM mixing |

Detailed documentation:

- [GitHub–Zenodo publication guide (Chinese)](docs/GITHUB_ZENODO_UPLOAD_GUIDE_CN.md)
- [Fig. 4 usage and interpretation notes](docs/FIG4_NOTES.md)
- [Fig. 5 usage and interpretation notes](docs/FIG5_NOTES.md)
- [Fig. 6 usage and interpretation notes](docs/FIG6_NOTES.md)
- [Input-data requirements](data/README.md)

## Software requirements

- Python 3.9 or later
- NumPy
- pandas
- Matplotlib
- SciPy

Install the required packages from the repository root:

    python -m venv .venv

Windows PowerShell:

    .venv\Scripts\Activate.ps1
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt

macOS or Linux:

    source .venv/bin/activate
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt

Before the archived release is created, the package versions used to generate the final manuscript figures should be recorded in requirements-lock.txt or environment.yml.

## Input data

The scripts currently use relative paths. Run each script from the code directory and place its input files at the exact locations described below.

### Fig. 4

Place the following file directly in the code directory:

    AOC-BS-HP-UHP Halogens data.csv

### Fig. 5

Create the following directory:

    code/data source/

Place these four files inside it:

    AOC-BS-HP-UHP halogens and H2O_data.csv
    mineral_modes.csv
    partition_coefficients.csv
    stages.csv

### Fig. 6

Place the following file directly in the code directory:

    Forearc serpentinite-Metasediment-BS-HP-UHP halogen_data.csv

Do not rename input files unless the corresponding path constants in the scripts are updated.

## Running the analyses

From the repository root:

    cd code

Run Fig. 4:

    python fig4_aoc_release_contributions.py

Run Fig. 5:

    python fig5_scenario1_fluid_mantle_mixing.py

Run Fig. 6:

    python fig6_scenario2_melange_mantle_mixing.py

The scripts may open interactive Matplotlib windows. Close the Fig. 4 window to allow execution to continue to CSV export. For non-interactive execution of Figs. 5 and 6, set SHOW_FIGURE or SHOW_FIGURES to False before running.

## Expected outputs

### Fig. 4

- F_ppm_stage_contributions.csv
- Cl_ppm_stage_contributions.csv
- Br_ppm_stage_contributions.csv
- I_ppm_stage_contributions.csv

The current Fig. 4 script displays the figure interactively but does not save it automatically.

### Fig. 5

- fluid_concentration_all_iterations.csv
- fluid_concentration_summary.csv
- density_50mantle_steps.png
- halogen_ratios_mantle_mixing_10_90.png

### Fig. 6

- 1_melange_and_mantle_mixing_three_by_two.svg
- halogen_ratios_melange_50-45-5_and_mantle_mixing_10-90.svg
- mantle_mixing_10-90_central95_median_1sigma.csv

Generated outputs should normally not be committed repeatedly. For the archived manuscript release, include either the exact final figure source data and outputs or clear instructions that reproduce them.

## Reproducibility

- Figs. 5 and 6 use a fixed random seed of 42.
- Fig. 4 currently has no fixed random seed; its bootstrap results will vary slightly among runs.
- The final release must state the exact input-data version, software versions and parameter sources.
- The numerical source data underlying the published figures should be provided as source data or in a permanent data repository.
- The meanings of KDE density levels, percentile trimming and error bars must be defined in the manuscript and figure legends.

## Code and data availability

Replace this section after the Zenodo release has been created:

The custom Python code used for the bootstrap Monte Carlo calculations, fluid-composition simulations, mélange construction, mantle-mixing calculations and visualization of halogen-ratio distributions is archived on Zenodo at https://doi.org/10.5281/zenodo.XXXXXXX. The version corresponding to the published analysis is v1.0.0. Input datasets and source data are available at [location and identifier].

## Citation

Before creating the GitHub release, complete CITATION.cff.template, rename it to CITATION.cff and replace all placeholders with the final authors, software title, version, release date, repository URL and DOI.

## License

No reuse license is granted merely by making a repository public. Select a license approved by all code owners before release. An MIT license template is included for consideration; replace the placeholders before renaming it to LICENSE.

## Contact

Correspondence and requests concerning the code should be directed to:

[Corresponding author name]  
[Institution]  
[Email address]
