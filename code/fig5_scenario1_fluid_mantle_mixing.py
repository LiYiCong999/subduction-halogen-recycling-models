from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde


# =============================================================================
# 1. Configuration
# =============================================================================

# Change this directory to the folder containing the four input CSV files.
BASE_DIR = Path(r"data source")
OUTPUT_DIR = BASE_DIR

BULK_DATA_FILE = BASE_DIR / "AOC-BS-HP-UHP halogens and H2O_data.csv"
MINERAL_MODES_FILE = BASE_DIR / "mineral_modes.csv"
PARTITION_COEFFICIENTS_FILE = BASE_DIR / "partition_coefficients.csv"
STAGES_FILE = BASE_DIR / "stages.csv"

ALL_RESULTS_FILE = OUTPUT_DIR / "fluid_concentration_all_iterations.csv"
SUMMARY_FILE = OUTPUT_DIR / "fluid_concentration_summary.csv"
FIGURE_FILE = OUTPUT_DIR / "density_50mantle_steps.png"
MIXING_RATIO_FIGURE_FILE = OUTPUT_DIR / "halogen_ratios_mantle_mixing_10_90.png"

# True: run Monte Carlo simulation and then plot.
# False: skip simulation and read ALL_RESULTS_FILE directly for plotting.
RUN_MONTE_CARLO = True
N_ITER = 100000
RANDOM_SEED = 42
SHOW_FIGURE = True

# Halogen column names and unit conversions to ppm.
ELEMENT_CONFIG = {
    "F": {"col": "F_ppm", "factor": 1.0},
    "Cl": {"col": "Cl_ppm", "factor": 1.0},
    "Br": {"col": "Br_ppb", "factor": 0.001},
    "I": {"col": "I_ppb", "factor": 0.001},
}
ELEMENTS = list(ELEMENT_CONFIG)

# DMM composition in ppm: F, Cl, Br, I.
MANTLE_CONC = np.array([12.0, 5.0, 0.013, 0.0003])
MANTLE_RATIO = 0.5
MANTLE_FRACTIONS = np.arange(0.1, 1.0, 0.1)
STEPS = [1, 2, 3]

# Plot settings.
FIG_SIZE = (18, 18)
GRID_RES = 100
DENSITY_LEVELS_PCT = [0.2, 0.4, 0.6, 0.8, 1.0]
X_LIM = (1e-6, 0.1)
Y_LIM_F = (0.01, 100)
Y_LIM_BR = {
    1: (0.001, 0.1),
    2: (0.001, 0.01),
    3: (0.001, 0.01),
}
HIST_BINS = 20
HIST_SIZE_RATIO = 0.4
INNER_HSPACE = 0.06
INNER_WSPACE = 0.06

STEP_COLORS = {
    1: (62 / 255, 158 / 255, 179 / 255),
    2: (241 / 255, 163 / 255, 97 / 255),
    3: (136 / 255, 176 / 255, 124 / 255),
}


# =============================================================================
# 2. Input and validation
# =============================================================================

def require_columns(df: pd.DataFrame, required, table_name: str) -> None:
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"{table_name} is missing columns: {sorted(missing)}")


def load_model_inputs():
    bulk_data = pd.read_csv(BULK_DATA_FILE)
    mineral_modes = pd.read_csv(MINERAL_MODES_FILE)
    partition_coeff = pd.read_csv(PARTITION_COEFFICIENTS_FILE)
    stages_def = pd.read_csv(STAGES_FILE)

    bulk_required = ["stage", "H2O_min_wt%", "H2O_max_wt%"] + [
        config["col"] for config in ELEMENT_CONFIG.values()
    ]
    require_columns(bulk_data, bulk_required, "bulk_data.csv")
    require_columns(
        mineral_modes,
        ["stage", "mineral", "min_vol%", "max_vol%"],
        "mineral_modes.csv",
    )
    require_columns(partition_coeff, ["mineral"], "partition_coefficients.csv")
    require_columns(
        stages_def,
        ["step", "initial_stage", "residual_stage"],
        "stages.csv",
    )
    return bulk_data, mineral_modes, partition_coeff, stages_def


def load_simulation_results(filepath: Path) -> pd.DataFrame:
    df = pd.read_csv(filepath)
    required = ["F", "Cl", "Br", "I", "step"]
    require_columns(df, required, filepath.name)
    df = df.dropna(subset=required).copy()
    df = df[
        (df["F"] > 0)
        & (df["Cl"] > 0)
        & (df["Br"] > 0)
        & (df["I"] > 0)
    ]
    return df


# =============================================================================
# 3. Monte Carlo sampling and fluid-composition calculation
# =============================================================================

def sample_water_range(df: pd.DataFrame, stage: str) -> float:
    selected = df[df["stage"] == stage]
    if selected.empty:
        raise ValueError(f"Stage not found in bulk_data.csv: {stage}")
    row = selected.iloc[0]
    return np.random.uniform(row["H2O_min_wt%"], row["H2O_max_wt%"])


def sample_halogen_concentration(
    df: pd.DataFrame, stage: str
) -> Dict[str, float]:
    selected = df[df["stage"] == stage]
    if selected.empty:
        raise ValueError(f"Stage not found in bulk_data.csv: {stage}")
    sample = selected.sample(1).iloc[0]
    result = {}
    for element, config in ELEMENT_CONFIG.items():
        value = sample[config["col"]]
        result[element] = (
            np.nan if pd.isna(value) else float(value) * config["factor"]
        )
    return result


def sample_mineral_mode(
    df: pd.DataFrame, stage: str
) -> Dict[str, float]:
    modes = df[df["stage"] == stage]
    if modes.empty:
        raise ValueError(f"Stage not found in mineral_modes.csv: {stage}")

    sampled = {}
    total = 0.0
    for _, row in modes.iterrows():
        value = np.random.uniform(row["min_vol%"], row["max_vol%"])
        sampled[row["mineral"]] = value
        total += value

    if total <= 0:
        raise ValueError(f"The sum of the randomly sampled mineral abundances for stage {stage} is not greater than 0")
    return {mineral: value / total for mineral, value in sampled.items()}


def sample_D_for_minerals(
    mineral_fractions: Dict[str, float],
    partition_coeff: pd.DataFrame,
) -> Dict[str, Dict[str, float]]:
    sampled_D = {}
    for mineral in mineral_fractions:
        selected = partition_coeff[partition_coeff["mineral"] == mineral]
        if selected.empty:
            sampled_D[mineral] = {element: 0.0 for element in ELEMENTS}
            continue

        row = selected.iloc[0]
        sampled_D[mineral] = {}
        for element in ELEMENTS:
            min_col = f"D_{element}_min"
            max_col = f"D_{element}_max"
            single_col = f"D_{element}"

            if min_col in row.index and max_col in row.index:
                d_value = np.random.uniform(row[min_col], row[max_col])
            elif single_col in row.index:
                d_value = row[single_col]
            else:
                d_value = 0.0
            sampled_D[mineral][element] = float(d_value)
    return sampled_D


def get_bulk_D(
    mineral_fractions: Dict[str, float],
    sampled_D: Dict[str, Dict[str, float]],
    element: str,
) -> float:
    return sum(
        fraction * sampled_D[mineral][element]
        for mineral, fraction in mineral_fractions.items()
    )


def calc_fluid_concentration(
    initial_concentration: float,
    initial_water: float,
    residual_water: float,
    mineral_fractions: Dict[str, float],
    sampled_D: Dict[str, Dict[str, float]],
    element: str,
) -> float:
    if np.isnan(initial_concentration):
        return np.nan

    fluid_fraction = initial_water - residual_water
    if fluid_fraction <= 1e-9:
        return np.nan

    bulk_D = get_bulk_D(mineral_fractions, sampled_D, element)
    denominator = fluid_fraction + (1.0 - fluid_fraction) * bulk_D
    if denominator <= 0:
        return np.nan
    return initial_concentration / denominator


def monte_carlo_simulation(
    bulk_data: pd.DataFrame,
    mineral_modes: pd.DataFrame,
    partition_coeff: pd.DataFrame,
    stages_def: pd.DataFrame,
    n_iter: int = N_ITER,
    seed: int = RANDOM_SEED,
) -> pd.DataFrame:
    np.random.seed(seed)
    step_values = sorted(stages_def["step"].astype(int).unique())
    results = {
        step: {element: [] for element in ELEMENTS}
        for step in step_values
    }

    for _ in range(n_iter):
        for _, stage_row in stages_def.iterrows():
            step = int(stage_row["step"])
            initial_stage = stage_row["initial_stage"]
            residual_stage = stage_row["residual_stage"]

            initial = sample_halogen_concentration(bulk_data, initial_stage)
            initial_water = sample_water_range(bulk_data, initial_stage)
            residual_water = sample_water_range(bulk_data, residual_stage)
            mineral_fractions = sample_mineral_mode(
                mineral_modes, residual_stage
            )
            sampled_D = sample_D_for_minerals(
                mineral_fractions, partition_coeff
            )

            for element in ELEMENTS:
                concentration = calc_fluid_concentration(
                    initial[element],
                    initial_water,
                    residual_water,
                    mineral_fractions,
                    sampled_D,
                    element,
                )
                results[step][element].append(concentration)

    frames = []
    for step, element_data in results.items():
        frame = pd.DataFrame(element_data)
        frame["step"] = step
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


# =============================================================================
# 4. Output simulation tables
# =============================================================================

def build_summary(sim_results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for step in sorted(sim_results["step"].unique()):
        subset = sim_results[sim_results["step"] == step]
        for element in ELEMENTS:
            data = subset[element].dropna()
            if data.empty:
                rows.append(
                    {
                        "step": step,
                        "element": element,
                        "mean": np.nan,
                        "median": np.nan,
                        "std": np.nan,
                        "p2.5": np.nan,
                        "p97.5": np.nan,
                        "p5": np.nan,
                        "p95": np.nan,
                        "count": 0,
                    }
                )
            else:
                rows.append(
                    {
                        "step": step,
                        "element": element,
                        "mean": data.mean(),
                        "median": data.median(),
                        "std": data.std(),
                        "p2.5": data.quantile(0.025),
                        "p97.5": data.quantile(0.975),
                        "p5": data.quantile(0.05),
                        "p95": data.quantile(0.95),
                        "count": len(data),
                    }
                )
    return pd.DataFrame(rows)


def print_summary(sim_results: pd.DataFrame) -> None:
    print("\n===== Fluid halogen concentration distribution statistics (ppm) =====")
    for step in sorted(sim_results["step"].unique()):
        print(f"\n--- Step {int(step)} ---")
        subset = sim_results[sim_results["step"] == step]
        for element in ELEMENTS:
            data = subset[element].dropna()
            if data.empty:
                print(f"{element}: no valid data")
            else:
                print(
                    f"{element}: mean={data.mean():.4f}, "
                    f"median={data.median():.4f}, "
                    f"std={data.std():.4f}, "
                    f"5%={data.quantile(0.05):.4f}, "
                    f"95%={data.quantile(0.95):.4f}"
                )


def save_simulation_outputs(sim_results: pd.DataFrame) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = build_summary(sim_results)
    summary.to_csv(SUMMARY_FILE, index=False)
    sim_results.to_csv(ALL_RESULTS_FILE, index=False)
    print(f"Statistical summary saved to: {SUMMARY_FILE}")
    print(f"Complete simulation results saved to: {ALL_RESULTS_FILE}")
    print_summary(sim_results)


# =============================================================================
# 5. KDE density plot and marginal histograms
# =============================================================================

def compute_kde_log(x, y, grid_res=GRID_RES):
    log_x = np.log10(x)
    log_y = np.log10(y)
    try:
        kde = gaussian_kde(np.vstack([log_x, log_y]))
    except np.linalg.LinAlgError:
        return None, None, None

    margin = 0.05
    x_range = log_x.max() - log_x.min() or 1.0
    y_range = log_y.max() - log_y.min() or 1.0
    xmin = log_x.min() - margin * x_range
    xmax = log_x.max() + margin * x_range
    ymin = log_y.min() - margin * y_range
    ymax = log_y.max() + margin * y_range

    X_log, Y_log = np.mgrid[
        xmin:xmax:complex(0, grid_res),
        ymin:ymax:complex(0, grid_res),
    ]
    positions = np.vstack([X_log.ravel(), Y_log.ravel()])
    density = np.reshape(kde(positions).T, X_log.shape)
    return 10**X_log, 10**Y_log, density


def plot_density_contours(ax, x, y, color, levels_pct):
    valid = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if len(x) < 5:
        return None

    X_orig, Y_orig, density = compute_kde_log(x, y)
    if density is None or density.max() == 0:
        return None

    levels = [percentage * density.max() for percentage in levels_pct]
    base_rgba = to_rgba(color)
    colors = [
        (*base_rgba[:3], levels_pct[index + 1])
        for index in range(len(levels) - 1)
    ]
    return ax.contourf(
        X_orig,
        Y_orig,
        density,
        levels=levels,
        colors=colors,
        antialiased=True,
    )


def add_marginal_histogram(
    ax_hist,
    data,
    bin_edges,
    color,
    orientation="vertical",
    add_kde=True,
):
    if len(data) == 0:
        return

    counts, _ = np.histogram(data, bins=bin_edges)
    if orientation == "vertical":
        ax_hist.bar(
            bin_edges[:-1],
            counts,
            width=np.diff(bin_edges),
            align="edge",
            color=color,
            alpha=0.5,
            edgecolor="none",
        )
        ax_hist.set_xlim(bin_edges[0], bin_edges[-1])
        ax_hist.set_ylabel("Count")
    else:
        ax_hist.barh(
            bin_edges[:-1],
            counts,
            height=np.diff(bin_edges),
            align="edge",
            color=color,
            alpha=0.5,
            edgecolor="none",
        )
        ax_hist.set_ylim(bin_edges[0], bin_edges[-1])
        ax_hist.set_xlabel("Count")

    if not add_kde or len(data) <= 2:
        return

    log_data = np.log10(data)
    try:
        kde = gaussian_kde(log_data)
    except np.linalg.LinAlgError:
        return

    log_start = np.log10(bin_edges[0])
    log_end = np.log10(bin_edges[-1])
    smooth_log = np.linspace(log_start, log_end, 200)
    density_log = kde(smooth_log)
    bin_width_log = (log_end - log_start) / (len(bin_edges) - 1)
    density_scaled = density_log * len(data) * bin_width_log

    if orientation == "vertical":
        ax_hist.plot(
            10**smooth_log, density_scaled, color="darkred", linewidth=1.5
        )
    else:
        ax_hist.plot(
            density_scaled, 10**smooth_log, color="darkred", linewidth=1.5
        )


def mix_with_mantle(subset: pd.DataFrame):
    original = subset[ELEMENTS].to_numpy(dtype=float)
    mixed = (1.0 - MANTLE_RATIO) * original + MANTLE_RATIO * MANTLE_CONC

    valid_cl = mixed[:, 1] > 0
    mixed = mixed[valid_cl]
    f_cl = mixed[:, 0] / mixed[:, 1]
    br_cl = mixed[:, 2] / mixed[:, 1]
    i_cl = mixed[:, 3] / mixed[:, 1]
    valid = (
        (f_cl > 0)
        & (br_cl > 0)
        & (i_cl > 0)
        & np.isfinite(f_cl)
        & np.isfinite(br_cl)
        & np.isfinite(i_cl)
    )
    return i_cl[valid], f_cl[valid], br_cl[valid]


def plot_results(sim_results: pd.DataFrame) -> None:
    sim_results = sim_results.dropna(
        subset=["F", "Cl", "Br", "I", "step"]
    ).copy()
    positive = (
        (sim_results["F"] > 0)
        & (sim_results["Cl"] > 0)
        & (sim_results["Br"] > 0)
        & (sim_results["I"] > 0)
    )
    sim_results = sim_results[positive]

    top_bins = np.logspace(
        np.log10(X_LIM[0]), np.log10(X_LIM[1]), HIST_BINS + 1
    )

    fig = plt.figure(figsize=FIG_SIZE)
    outer_grid = GridSpec(
        3,
        2,
        figure=fig,
        hspace=0.4,
        wspace=0.35,
        left=0.08,
        right=0.95,
        top=0.95,
        bottom=0.05,
    )

    for row, step in enumerate(STEPS):
        subset = sim_results[sim_results["step"] == step]
        if subset.empty:
            print(f"Warning: no valid data for Step {step}; plotting skipped.")
            continue

        x, y_f, y_br = mix_with_mantle(subset)
        step_color = STEP_COLORS[step]
        ylim_f = Y_LIM_F
        ylim_br = Y_LIM_BR[step]

        right_bins_f = np.logspace(
            np.log10(ylim_f[0]), np.log10(ylim_f[1]), HIST_BINS + 1
        )
        right_bins_br = np.logspace(
            np.log10(ylim_br[0]), np.log10(ylim_br[1]), HIST_BINS + 1
        )
        column_configs = [
            (0, y_f, "F/Cl", ylim_f, right_bins_f),
            (1, y_br, "Br/Cl", ylim_br, right_bins_br),
        ]

        for column, y, y_label, y_lim, right_bins in column_configs:
            inner_grid = outer_grid[row, column].subgridspec(
                2,
                2,
                width_ratios=[4, HIST_SIZE_RATIO],
                height_ratios=[HIST_SIZE_RATIO, 4],
                hspace=INNER_HSPACE,
                wspace=INNER_WSPACE,
            )

            ax_main = fig.add_subplot(inner_grid[1, 0])
            ax_main.set_xscale("log")
            ax_main.set_yscale("log")
            ax_main.set_xlim(X_LIM)
            ax_main.set_ylim(y_lim)
            ax_main.grid(True, linestyle="--", alpha=0.5)
            ax_main.set_box_aspect(1)
            ax_main.set_xlabel("I/Cl")
            ax_main.set_ylabel(y_label)
            ax_main.set_title(
                f"{y_label} vs I/Cl (Step {step}, 50% Mantle)",
                fontsize=9,
            )

            ax_top = fig.add_subplot(inner_grid[0, 0], sharex=ax_main)
            ax_top.set_xscale("log")
            ax_top.set_xlim(X_LIM)
            ax_top.tick_params(
                axis="x", labelbottom=False, bottom=False
            )
            add_marginal_histogram(
                ax_top,
                x,
                top_bins,
                step_color,
                orientation="vertical",
                add_kde=True,
            )

            ax_right = fig.add_subplot(inner_grid[1, 1], sharey=ax_main)
            ax_right.set_yscale("log")
            ax_right.set_ylim(y_lim)
            ax_right.tick_params(
                axis="y", labelleft=False, left=False
            )
            add_marginal_histogram(
                ax_right,
                y,
                right_bins,
                step_color,
                orientation="horizontal",
                add_kde=True,
            )

            ax_blank = fig.add_subplot(inner_grid[0, 1])
            ax_blank.axis("off")
            plot_density_contours(
                ax_main, x, y, step_color, DENSITY_LEVELS_PCT
            )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_FILE, dpi=300, bbox_inches="tight")
    print(f"Density plot saved to: {FIGURE_FILE}")
    if not SHOW_FIGURE:
        plt.close(fig)


# =============================================================================
# 6. Halogen-ratio evolution from 10% to 90% mantle mixing
# =============================================================================

def summarize_central_95(values):
    """
    Retain the central 95% of a simulated distribution and calculate its
    median and sample standard deviation. The median is used as the plotted
    central value and ±1σ is used as the error-bar range.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return np.nan, np.nan

    lower, upper = np.quantile(values, [0.025, 0.975])
    central = values[(values >= lower) & (values <= upper)]
    if len(central) < 2:
        return np.nan, np.nan
    return float(np.median(central)), float(np.std(central, ddof=1))


def calculate_mixing_ratio_statistics(sim_results: pd.DataFrame) -> pd.DataFrame:
    """Calculate median and ±1σ for each step and mantle fraction."""
    rows = []
    ratio_definitions = {
        "F/Cl": (0, 1),
        "Br/Cl": (2, 1),
        "I/Cl": (3, 1),
    }

    for step in STEPS:
        subset = sim_results[sim_results["step"] == step]
        if subset.empty:
            continue

        fluid = subset[ELEMENTS].to_numpy(dtype=float)
        valid_fluid = np.all(np.isfinite(fluid), axis=1) & np.all(fluid > 0, axis=1)
        fluid = fluid[valid_fluid]

        for mantle_fraction in MANTLE_FRACTIONS:
            mixed = (
                (1.0 - mantle_fraction) * fluid
                + mantle_fraction * MANTLE_CONC
            )

            for ratio_name, (numerator, denominator) in ratio_definitions.items():
                valid = mixed[:, denominator] > 0
                ratios = mixed[valid, numerator] / mixed[valid, denominator]
                ratios = ratios[(ratios > 0) & np.isfinite(ratios)]
                median, std = summarize_central_95(ratios)
                rows.append(
                    {
                        "step": step,
                        "mantle_fraction": mantle_fraction,
                        "ratio": ratio_name,
                        "median": median,
                        "std": std,
                    }
                )

    return pd.DataFrame(rows)


def plot_mantle_fraction_ratios(sim_results: pd.DataFrame) -> None:
    """
    Plot F/Cl, Br/Cl and I/Cl against mantle fraction in a 1 × 3 layout.
    Values are medians of the central 95% of each distribution; vertical
    error bars show ±1σ calculated from that retained distribution.
    """
    statistics = calculate_mixing_ratio_statistics(sim_results)
    if statistics.empty:
        print("Warning: no valid data are available for the mantle-mixing fraction plot.")
        return

    ratio_names = ["F/Cl", "Br/Cl", "I/Cl"]
    x_centers = np.arange(1, len(MANTLE_FRACTIONS) + 1)
    step_offsets = {1: -0.2, 2: 0.0, 3: 0.2}

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, ratio_name in zip(axes, ratio_names):
        for step in STEPS:
            selected = statistics[
                (statistics["step"] == step)
                & (statistics["ratio"] == ratio_name)
            ].sort_values("mantle_fraction")
            if selected.empty:
                continue

            median = selected["median"].to_numpy(dtype=float)
            std = selected["std"].to_numpy(dtype=float)
            x = x_centers + step_offsets[step]
            color = STEP_COLORS[step]

            ax.errorbar(
                x,
                median,
                yerr=std,
                fmt="o-",
                color=color,
                ecolor=color,
                linewidth=1.8,
                elinewidth=1.8,
                capsize=4,
                capthick=1.0,
                markersize=8,
                markeredgecolor="#333333",
                markeredgewidth=0.6,
                label=f"Step {step} ±1σ",
                zorder=3,
            )

            for x_value, median_value in zip(x, median):
                if not np.isfinite(median_value):
                    continue
                ax.annotate(
                    f"{median_value:.3g}",
                    xy=(x_value, median_value),
                    xytext=(0, 12),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color=color,
                    bbox={
                        "facecolor": "white",
                        "edgecolor": "none",
                        "alpha": 0.7,
                        "pad": 1.0,
                    },
                    zorder=4,
                )

        ax.set_title(f"{ratio_name} ratio", fontsize=15)
        ax.set_xlabel("Mantle fraction", fontsize=12)
        ax.set_ylabel(ratio_name, fontsize=12)
        ax.set_xticks(x_centers)
        ax.set_xticklabels(
            [f"{fraction:.0%}" for fraction in MANTLE_FRACTIONS],
            rotation=45,
            ha="right",
        )
        ax.set_xlim(0.4, len(MANTLE_FRACTIONS) + 0.6)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.margins(y=0.08)
        ax.legend(loc="upper right", fontsize=11)

    fig.suptitle(
        "Halogen ratios with mantle mixing (10%–90%)",
        fontsize=16,
        y=1.01,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(MIXING_RATIO_FIGURE_FILE, dpi=300, bbox_inches="tight")
    print(f"Mantle-mixing fraction plot saved to: {MIXING_RATIO_FIGURE_FILE}")
    if not SHOW_FIGURE:
        plt.close(fig)


# =============================================================================
# 7. One-click workflow
# =============================================================================

def main():
    if RUN_MONTE_CARLO:
        print(f"Running Monte Carlo simulation with {N_ITER} iterations ...")
        inputs = load_model_inputs()
        sim_results = monte_carlo_simulation(
            *inputs,
            n_iter=N_ITER,
            seed=RANDOM_SEED,
        )
        save_simulation_outputs(sim_results)
    else:
        print(f"Skipping simulation and reading existing results from: {ALL_RESULTS_FILE}")
        sim_results = load_simulation_results(ALL_RESULTS_FILE)

    print("\nStarting 50% mantle mixing and plotting ...")
    plot_results(sim_results)

    print("\nStarting the 10%–90% mantle-mixing fraction plot ...")
    plot_mantle_fraction_ratios(sim_results)

    if SHOW_FIGURE:
        plt.show()
    print("\nAll calculations and plots have been completed.")


if __name__ == "__main__":
    main()
