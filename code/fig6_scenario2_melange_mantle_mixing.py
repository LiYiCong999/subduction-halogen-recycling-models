
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm, to_rgba
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import gaussian_kde


# =============================================================================
# 1. User configuration
# =============================================================================

# Read only this original dataset; do not read any mixed_samples_*.csv files.
HALOGEN_DATA_FILE = Path(
    r"Forearc serpentinite-Metasediment-BS-HP-UHP halogen_data.csv"
)

# Lithological group column in the original table. Set to None for automatic
# detection. If automatic detection fails, enter the actual column name,
# such as "Rock_type", "Lithology" or "Group".
GROUP_COLUMN = None

# Three-endmember Monte Carlo mixing settings.
N_MONTE_CARLO = 100000
RANDOM_SEED = 42
MELANGE_PROPORTIONS = {
    "Serpentinite": 0.50,
    "Metasediment": 0.05,
    "Metamorphic_endmember": 0.45,
}

# Common names for each lithology in the original table. The program first
# performs exact matching and then keyword matching. If halogen_data.csv uses
# another name, add it to the corresponding tuple.
GROUP_ALIASES = {
    "Serpentinite": (
        "Serpentinite", "Serpentine", "Serp", "Forearc serpentinite",
        "Serpentinite rock", "Forearc serpentinite rock",
    ),
    "Metasediment": (
        "Metasediment", "Metasedimentary rock", "Meta", "Sediment",
        "Metasedimentary material", "Metamorphosed sediment", "Sedimentary material",
    ),
    "Blue schist": (
        "Blue schist", "Blueschist", "Blue", "Blueschist rock",
    ),
    "HP eclogite": (
        "HP eclogite", "High-pressure eclogite", "Eclogite", "Ecl",
        "High pressure eclogite", "Eclogite rock",
    ),
    "UHP eclogite": (
        "UHP eclogite", "Ultrahigh-pressure eclogite", "UHP Ecl",
        "Ultrahigh pressure eclogite",
    ),
}

# Save results to the directory containing this script by default.
OUTPUT_DIR = Path(__file__).resolve().parent

# Export both final visualizations as SVG vector graphics.
DENSITY_FIGURE_FILE = OUTPUT_DIR / "1_melange_and_mantle_mixing_three_by_two.svg"
RATIO_FIGURE_FILE = OUTPUT_DIR / "halogen_ratios_melange_50-45-5_and_mantle_mixing_10-90.svg"

# The statistical table is separate from the figures; set to False if unneeded.
EXPORT_STATS_CSV = True
STATS_FILE = OUTPUT_DIR / "mantle_mixing_10-90_central95_median_1sigma.csv"

# Whether to display the figure windows after completion.
SHOW_FIGURES = True


# =============================================================================
# 2. Shared constants and plotting parameters
# =============================================================================

ELEMENTS = ["F", "Cl", "Br", "I"]
RATIO_NAMES = ["F/Cl", "Br/Cl", "I/Cl"]
ROCK_ORDER = ["Blue schist", "HP eclogite", "UHP eclogite"]

# DMM endmember concentrations (ppm): F, Cl, Br and I.
MANTLE_CONC = np.array([12.0, 5.0, 0.013, 0.0003], dtype=float)

# Figure 1 uses fixed 50% mantle mixing; Figure 2 uses 10%–90% in 10% steps.
DENSITY_MANTLE_FRACTION = 0.50
MANTLE_FRACTIONS = np.arange(0.10, 1.00, 0.10)

# Shared colors for the three sample groups, retaining the original palette.
COLORS = {
    "Blue schist": (62 / 255, 158 / 255, 179 / 255),
    "HP eclogite": (241 / 255, 163 / 255, 97 / 255),
    "UHP eclogite": (136 / 255, 176 / 255, 124 / 255),
}

ROW_TITLES = {
    "Blue schist": "Mantle domain (1)",
    "HP eclogite": "Mantle domain (2)",
    "UHP eclogite": "Mantle domain (3)",
}

# Figure 1 axis ranges, consistent with the reference six-panel figure.
X_LIM = (1e-6, 1e-1)
Y_LIM_F_CL = (1e-2, 1e1)
Y_LIM_BR_CL = (1e-3, 1e-1)

# Figure 1 KDE and marginal-histogram settings.
GRID_RES = 100
DENSITY_LEVELS = np.array([0.20, 0.40, 0.60, 0.80, 1.00])
HIST_BINS = 20
FIGURE1_SIZE = (16, 23)

# Preserve text objects in SVG files for editing in Adobe Illustrator, Inkscape, etc.
plt.rcParams["svg.fonttype"] = "none"

# Figure 2 settings.
FIGURE2_SIZE = (18, 6)
ROCK_X_OFFSETS = {
    "Blue schist": -0.20,
    "HP eclogite": 0.00,
    "UHP eclogite": 0.20,
}


# =============================================================================
# 3. Data loading, validation and mixing calculations
# =============================================================================

def normalize_text(value) -> str:
    """Normalize column and group names, ignoring case, spaces, underscores and hyphens."""
    return re.sub(r"[^0-9a-zA-Z\u4e00-\u9fff]+", "", str(value).strip().lower())


def detect_group_column(df: pd.DataFrame) -> str:
    """Automatically detect the lithology or group column in halogen_data.csv."""
    if GROUP_COLUMN is not None:
        if GROUP_COLUMN not in df.columns:
            raise ValueError(
                f"GROUP_COLUMN={GROUP_COLUMN!r} is not present in the original table; "
                f"available columns: {list(df.columns)}"
            )
        return GROUP_COLUMN

    candidates = (
        "Rock_type", "Rock type", "Rock", "Lithology", "Type", "Group",
        "Facies", "Stage", "Category", "Endmember", "Sample_type",
        "Lithological_unit", "Data_type", "Sample_group", "Metamorphic_facies", "End_member",
    )
    normalized_columns = {normalize_text(col): col for col in df.columns}
    for candidate in candidates:
        matched = normalized_columns.get(normalize_text(candidate))
        if matched is not None:
            return matched

    raise ValueError(
        "The lithology or group column could not be detected automatically. "
        "Set GROUP_COLUMN at the beginning of the script to the actual column name."
        f"\nColumns in the original table: {list(df.columns)}"
    )


def detect_element_columns(df: pd.DataFrame) -> dict[str, tuple[str, float]]:
    """
    Identify the F, Cl, Br and I data columns.

    Return {element: (original column name, conversion factor to ppm)}.
    Multiply Br and I columns reported in ppb by 0.001.
    """
    alias_priority = {
        "F": ("F_ppm", "F (ppm)", "Fppm", "F"),
        "Cl": ("Cl_ppm", "Cl (ppm)", "Clppm", "Cl"),
        "Br": ("Br_ppm", "Br (ppm)", "Brppm", "Br_ppb", "Br (ppb)", "Brppb", "Br"),
        "I": ("I_ppm", "I (ppm)", "Ippm", "I_ppb", "I (ppb)", "Ippb", "I"),
    }
    normalized_columns = {normalize_text(col): col for col in df.columns}
    result = {}
    for element, aliases in alias_priority.items():
        matched = None
        for alias in aliases:
            matched = normalized_columns.get(normalize_text(alias))
            if matched is not None:
                break
        if matched is None:
            raise ValueError(
                f"The {element} data column could not be identified. "
                f"Columns in the original table: {list(df.columns)}"
            )
        normalized_name = normalize_text(matched)
        factor = 0.001 if normalized_name.endswith("ppb") else 1.0
        result[element] = (matched, factor)
    return result


def load_halogen_data(filepath: Path) -> tuple[pd.DataFrame, str]:
    """Read the original halogen_data.csv and convert all four elements to ppm."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Original data file not found: {filepath}")

    raw = pd.read_csv(filepath)
    if raw.empty:
        raise ValueError(f"{filepath.name} contains no data.")

    group_column = detect_group_column(raw)
    element_columns = detect_element_columns(raw)
    data = pd.DataFrame({group_column: raw[group_column].astype(str).str.strip()})
    for element in ELEMENTS:
        source_column, ppm_factor = element_columns[element]
        data[element] = pd.to_numeric(raw[source_column], errors="coerce") * ppm_factor

    data = data.replace([np.inf, -np.inf], np.nan).dropna(subset=ELEMENTS)
    data = data[(data[ELEMENTS] > 0).all(axis=1)].copy()
    if data.empty:
        raise ValueError(
            "The original table contains no valid records with positive "
            "F, Cl, Br and I concentrations."
        )

    print(f"Original file loaded: {filepath}")
    print(f"Identified group column: {group_column}")
    print("All four elements have been converted to ppm.")
    return data, group_column


def find_group_rows(
    data: pd.DataFrame,
    group_column: str,
    group_key: str,
) -> pd.DataFrame:
    """Extract one lithological endmember using GROUP_ALIASES."""
    values = data[group_column].astype(str)
    normalized_values = values.map(normalize_text)
    normalized_aliases = [normalize_text(v) for v in GROUP_ALIASES[group_key]]

    # Prioritize exact matching to prevent Eclogite from matching UHP eclogite.
    mask = normalized_values.isin(normalized_aliases)
    if not mask.any():
        aliases_by_length = sorted(normalized_aliases, key=len, reverse=True)
        mask = normalized_values.map(
            lambda value: any(alias and alias in value for alias in aliases_by_length)
        )
        if group_key == "HP eclogite":
            uhp_aliases = [normalize_text(v) for v in GROUP_ALIASES["UHP eclogite"]]
            mask &= ~normalized_values.map(
                lambda value: any(alias and alias in value for alias in uhp_aliases)
            )

    selected = data.loc[mask, ELEMENTS].copy()
    if selected.empty:
        available = sorted(values.dropna().unique().tolist())
        raise ValueError(
            f"{group_key!r} could not be found in {group_column!r}."
            f"\nAvailable group values: {available}"
            f"\nAdd the actual name to GROUP_ALIASES[{group_key!r}] "
            "at the beginning of the script."
        )
    return selected


def monte_carlo_three_endmember_mixing(
    raw_data: pd.DataFrame,
    group_column: str,
    n_samples: int = N_MONTE_CARLO,
    seed: int = RANDOM_SEED,
) -> dict[str, pd.DataFrame]:
    """
    Generate three 50%–5%–45% mélange compositions by random sampling
    directly from halogen_data.csv.

    In each iteration, sample one row with replacement from the serpentinite,
    metasediment and target metamorphic endmembers and linearly mix the F, Cl,
    Br and I concentrations. The workflow does not read or depend on
    premixed CSV files.
    """
    if n_samples < 5:
        raise ValueError("N_MONTE_CARLO must be at least 5.")

    proportions = np.array(
        [
            MELANGE_PROPORTIONS["Serpentinite"],
            MELANGE_PROPORTIONS["Metasediment"],
            MELANGE_PROPORTIONS["Metamorphic_endmember"],
        ],
        dtype=float,
    )
    if not np.isclose(proportions.sum(), 1.0):
        raise ValueError(
            f"MELANGE_PROPORTIONS must sum to 1; the current sum is "
            f"{proportions.sum():.6g}."
        )
    if (proportions < 0).any():
        raise ValueError("MELANGE_PROPORTIONS cannot contain negative values.")

    serpentinite = find_group_rows(
        raw_data, group_column, "Serpentinite"
    )[ELEMENTS].to_numpy(dtype=float)
    metasediment = find_group_rows(
        raw_data, group_column, "Metasediment"
    )[ELEMENTS].to_numpy(dtype=float)

    rng = np.random.default_rng(seed)
    output = {}
    for rock_name in ROCK_ORDER:
        metamorphic = find_group_rows(
            raw_data, group_column, rock_name
        )[ELEMENTS].to_numpy(dtype=float)

        serpentinite_sample = serpentinite[
            rng.integers(0, len(serpentinite), size=n_samples)
        ]
        metasediment_sample = metasediment[
            rng.integers(0, len(metasediment), size=n_samples)
        ]
        metamorphic_sample = metamorphic[
            rng.integers(0, len(metamorphic), size=n_samples)
        ]

        mixed = (
            proportions[0] * serpentinite_sample
            + proportions[1] * metasediment_sample
            + proportions[2] * metamorphic_sample
        )
        output[rock_name] = pd.DataFrame(mixed, columns=ELEMENTS)
        print(
            f"{rock_name}: generated {n_samples:,} "
            "50%–5%–45% mélange samples from the original data."
        )
    return output


def load_all_data() -> dict[str, pd.DataFrame]:
    """Read halogen_data.csv and generate three mélange datasets in memory."""
    raw_data, group_column = load_halogen_data(HALOGEN_DATA_FILE)
    return monte_carlo_three_endmember_mixing(raw_data, group_column)


def mix_with_mantle(values: np.ndarray, mantle_fraction: float) -> np.ndarray:
    """Linearly mix concentrations: fluid/rock endmember (1-r) and DMM (r)."""
    r = float(mantle_fraction)
    if not 0.0 <= r <= 1.0:
        raise ValueError("The mantle-mixing fraction must be between 0 and 1.")
    return (1.0 - r) * values + r * MANTLE_CONC.reshape(1, -1)


def calculate_ratios(mixed_values: np.ndarray) -> dict[str, np.ndarray]:
    """Calculate F/Cl, Br/Cl and I/Cl from F, Cl, Br and I concentrations."""
    valid = np.isfinite(mixed_values).all(axis=1) & (mixed_values[:, 1] > 0)
    mixed = mixed_values[valid]
    ratios = {
        "F/Cl": mixed[:, 0] / mixed[:, 1],
        "Br/Cl": mixed[:, 2] / mixed[:, 1],
        "I/Cl": mixed[:, 3] / mixed[:, 1],
    }
    for name, values in ratios.items():
        values = np.asarray(values, dtype=float)
        ratios[name] = values[np.isfinite(values) & (values > 0)]
    return ratios


# =============================================================================
# 4. Figure 1: three-by-two density plots
# =============================================================================

def compute_kde_log2d(x: np.ndarray, y: np.ndarray):
    """Calculate a two-dimensional KDE in log10(x)–log10(y) space."""
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    log_x = np.log10(x[mask])
    log_y = np.log10(y[mask])
    if len(log_x) < 5:
        return None

    try:
        kde = gaussian_kde(np.vstack([log_x, log_y]))
    except (np.linalg.LinAlgError, ValueError):
        return None

    x_span = max(float(np.ptp(log_x)), 1.0)
    y_span = max(float(np.ptp(log_y)), 1.0)
    x_min = log_x.min() - 0.05 * x_span
    x_max = log_x.max() + 0.05 * x_span
    y_min = log_y.min() - 0.05 * y_span
    y_max = log_y.max() + 0.05 * y_span

    grid_x, grid_y = np.mgrid[
        x_min:x_max:complex(0, GRID_RES),
        y_min:y_max:complex(0, GRID_RES),
    ]
    positions = np.vstack([grid_x.ravel(), grid_y.ravel()])
    density = kde(positions).reshape(grid_x.shape)
    return 10**grid_x, 10**grid_y, density


def plot_density_contours(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    base_color,
):
    """Plot continuous KDE regions using the original colors and 20%–100% relative-density levels."""
    result = compute_kde_log2d(x, y)
    if result is None:
        return None
    grid_x, grid_y, density = result
    density_max = float(np.nanmax(density))
    if not np.isfinite(density_max) or density_max <= 0:
        return None

    absolute_levels = DENSITY_LEVELS * density_max
    rgba = to_rgba(base_color)
    fill_colors = [
        (*rgba[:3], float(alpha)) for alpha in DENSITY_LEVELS[1:]
    ]
    return ax.contourf(
        grid_x,
        grid_y,
        density,
        levels=absolute_levels,
        colors=fill_colors,
        antialiased=True,
        zorder=1,
    )


def add_marginal_histogram(
    ax: plt.Axes,
    data: np.ndarray,
    bin_edges: np.ndarray,
    base_color,
    orientation: str,
) -> None:
    """Plot a log-binned histogram with a frequency-scaled one-dimensional log-KDE curve."""
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data) & (data > 0)]
    if len(data) == 0:
        return

    counts, _ = np.histogram(data, bins=bin_edges)
    if orientation == "vertical":
        ax.bar(
            bin_edges[:-1],
            counts,
            width=np.diff(bin_edges),
            align="edge",
            color=base_color,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.6,
        )
    elif orientation == "horizontal":
        ax.barh(
            bin_edges[:-1],
            counts,
            height=np.diff(bin_edges),
            align="edge",
            color=base_color,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.6,
        )
    else:
        raise ValueError("orientation must be either 'vertical' or 'horizontal'.")

    if len(data) <= 2 or np.ptp(np.log10(data)) == 0:
        return
    try:
        kde = gaussian_kde(np.log10(data))
    except (np.linalg.LinAlgError, ValueError):
        return

    z = np.linspace(np.log10(bin_edges[0]), np.log10(bin_edges[-1]), 240)
    delta_z = (z[-1] - z[0]) / (len(bin_edges) - 1)
    scaled_density = kde(z) * len(data) * delta_z
    if orientation == "vertical":
        ax.plot(10**z, scaled_density, color="darkred", linewidth=1.2)
    else:
        ax.plot(scaled_density, 10**z, color="darkred", linewidth=1.2)


def style_main_axis(
    ax: plt.Axes,
    y_label: str,
    y_limits: tuple[float, float],
    panel_letter: str,
    row_title: str,
) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(X_LIM)
    ax.set_ylim(y_limits)
    ax.set_xlabel("I/Cl", fontsize=13)
    ax.set_ylabel(y_label, fontsize=13)
    ax.set_box_aspect(1)
    ax.grid(False)
    ax.tick_params(which="both", direction="in", length=5, width=1.0)
    ax.tick_params(which="major", labelsize=11)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    ax.text(
        0.055,
        0.965,
        panel_letter,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=14,
        fontweight="bold",
        zorder=10,
    )
    ax.text(
        0.50,
        0.965,
        row_title,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
        zorder=10,
    )


def style_top_histogram(ax: plt.Axes) -> None:
    ax.set_xscale("log")
    ax.set_xlim(X_LIM)
    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Count", fontsize=8, labelpad=2)
    ax.tick_params(axis="y", labelsize=8, pad=1)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)


def style_right_histogram(ax: plt.Axes, y_limits: tuple[float, float]) -> None:
    ax.set_yscale("log")
    ax.set_ylim(y_limits)
    ax.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("Count", fontsize=8, labelpad=2)
    ax.tick_params(axis="x", labelsize=8, pad=1)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)


def add_density_scale(ax: plt.Axes, base_color) -> None:
    """Add a horizontal relative-density scale to the lower left of each left panel."""
    rgba = to_rgba(base_color)
    colors = [(*rgba[:3], float(alpha)) for alpha in DENSITY_LEVELS]
    cmap = ListedColormap(colors)
    boundaries = np.linspace(0.0, 1.0, len(colors) + 1)
    norm = BoundaryNorm(boundaries, cmap.N)

    cax = inset_axes(
        ax,
        width="47%",
        height="2.3%",
        loc="lower left",
        bbox_to_anchor=(0.035, 0.075, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = ax.figure.colorbar(
        scalar,
        cax=cax,
        orientation="horizontal",
        boundaries=boundaries,
        ticks=np.linspace(0, 1, 6),
    )
    cbar.set_ticklabels(["0", "20%", "40%", "60%", "80%", "100%"])
    cbar.ax.tick_params(labelsize=6, length=2, pad=1)
    cbar.ax.set_title("Density (fraction of max)", fontsize=6, pad=1)
    cbar.outline.set_linewidth(0.6)


def create_density_figure(data_by_rock: dict[str, pd.DataFrame]) -> plt.Figure:
    """Create a three-by-two density figure matching the reference layout."""
    fig = plt.figure(figsize=FIGURE1_SIZE)
    outer = GridSpec(
        3,
        2,
        figure=fig,
        left=0.065,
        right=0.965,
        bottom=0.045,
        top=0.985,
        wspace=0.25,
        hspace=0.18,
    )

    x_bins = np.logspace(np.log10(X_LIM[0]), np.log10(X_LIM[1]), HIST_BINS + 1)
    y_bins_f = np.logspace(
        np.log10(Y_LIM_F_CL[0]), np.log10(Y_LIM_F_CL[1]), HIST_BINS + 1
    )
    y_bins_br = np.logspace(
        np.log10(Y_LIM_BR_CL[0]), np.log10(Y_LIM_BR_CL[1]), HIST_BINS + 1
    )
    panel_letters = iter("ABCDEF")

    for row, rock_name in enumerate(ROCK_ORDER):
        original = data_by_rock[rock_name][ELEMENTS].to_numpy(dtype=float)
        mixed = mix_with_mantle(original, DENSITY_MANTLE_FRACTION)
        ratios = calculate_ratios(mixed)
        x = ratios["I/Cl"]
        row_color = COLORS[rock_name]

        panel_config = [
            ("F/Cl", ratios["F/Cl"], Y_LIM_F_CL, y_bins_f),
            ("Br/Cl", ratios["Br/Cl"], Y_LIM_BR_CL, y_bins_br),
        ]

        for col, (y_name, y, y_limits, y_bins) in enumerate(panel_config):
            # All ratios originate from the same mixed array and should have
            # equal lengths; apply a joint filter again for additional robustness.
            n = min(len(x), len(y))
            xx = x[:n]
            yy = y[:n]
            valid = np.isfinite(xx) & np.isfinite(yy) & (xx > 0) & (yy > 0)
            xx, yy = xx[valid], yy[valid]

            inner = outer[row, col].subgridspec(
                2,
                2,
                width_ratios=[4.0, 0.42],
                height_ratios=[0.42, 4.0],
                hspace=0.035,
                wspace=0.035,
            )
            ax_top = fig.add_subplot(inner[0, 0])
            ax_main = fig.add_subplot(inner[1, 0], sharex=ax_top)
            ax_right = fig.add_subplot(inner[1, 1], sharey=ax_main)
            ax_blank = fig.add_subplot(inner[0, 1])
            ax_blank.axis("off")

            style_main_axis(
                ax_main,
                y_name,
                y_limits,
                next(panel_letters),
                ROW_TITLES[rock_name],
            )
            plot_density_contours(ax_main, xx, yy, row_color)

            style_top_histogram(ax_top)
            add_marginal_histogram(
                ax_top, xx, x_bins, row_color, orientation="vertical"
            )

            style_right_histogram(ax_right, y_limits)
            add_marginal_histogram(
                ax_right, yy, y_bins, row_color, orientation="horizontal"
            )

            if col == 0:
                add_density_scale(ax_main, row_color)

    return fig


# =============================================================================
# 5. Figure 2: 10%–90% mantle mixing and three halogen ratios
# =============================================================================

def central_95_stats(values: np.ndarray) -> dict[str, float]:
    """
    Retain the central 95% of the data (2.5th–97.5th percentiles), then
    calculate the median and sample standard deviation.

    Figure points = medians; error bars = median ± 1σ.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n_total = len(values)
    if n_total == 0:
        return {
            "median": np.nan,
            "std": np.nan,
            "q025": np.nan,
            "q975": np.nan,
            "n_total": 0,
            "n_used": 0,
        }

    q025, q975 = np.percentile(values, [2.5, 97.5])
    retained = values[(values >= q025) & (values <= q975)]
    if len(retained) == 0:
        retained = values

    std = np.std(retained, ddof=1) if len(retained) > 1 else 0.0
    return {
        "median": float(np.median(retained)),
        "std": float(std),
        "q025": float(q025),
        "q975": float(q975),
        "n_total": int(n_total),
        "n_used": int(len(retained)),
    }


def compute_ratio_statistics(
    data_by_rock: dict[str, pd.DataFrame],
) -> tuple[dict, pd.DataFrame]:
    """Calculate central-95% statistics for three ratios, three groups and nine mantle fractions."""
    all_stats = {}
    summary_rows = []

    for rock_name in ROCK_ORDER:
        original = data_by_rock[rock_name][ELEMENTS].to_numpy(dtype=float)
        all_stats[rock_name] = {
            ratio_name: {"median": [], "std": []}
            for ratio_name in RATIO_NAMES
        }

        for fraction in MANTLE_FRACTIONS:
            ratios = calculate_ratios(mix_with_mantle(original, fraction))
            for ratio_name in RATIO_NAMES:
                stats = central_95_stats(ratios[ratio_name])
                all_stats[rock_name][ratio_name]["median"].append(stats["median"])
                all_stats[rock_name][ratio_name]["std"].append(stats["std"])
                summary_rows.append(
                    {
                        "Rock_type": rock_name,
                        "Mantle_fraction": fraction,
                        "Mantle_fraction_percent": int(round(fraction * 100)),
                        "Ratio": ratio_name,
                        "Median": stats["median"],
                        "Std_1sigma": stats["std"],
                        "Central95_lower": stats["q025"],
                        "Central95_upper": stats["q975"],
                        "N_total": stats["n_total"],
                        "N_used": stats["n_used"],
                    }
                )

        for ratio_name in RATIO_NAMES:
            for key in ("median", "std"):
                all_stats[rock_name][ratio_name][key] = np.asarray(
                    all_stats[rock_name][ratio_name][key], dtype=float
                )

    return all_stats, pd.DataFrame(summary_rows)


def create_ratio_figure(all_stats: dict) -> plt.Figure:
    """Create a one-by-three ratio–mantle-fraction figure matching the reference layout."""
    x = np.arange(len(MANTLE_FRACTIONS), dtype=float)
    x_labels = [f"{int(round(f * 100))}%" for f in MANTLE_FRACTIONS]

    fig, axes = plt.subplots(1, 3, figsize=FIGURE2_SIZE)
    for ratio_name, ax in zip(RATIO_NAMES, axes):
        ax.set_title(f"{ratio_name} ratio", fontsize=14)
        ax.set_xlabel("Mantle fraction")
        ax.set_ylabel(ratio_name)
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=45, ha="right")
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.set_box_aspect(1)

        lower_candidates = []
        upper_candidates = []
        for rock_name in ROCK_ORDER:
            median = all_stats[rock_name][ratio_name]["median"]
            std = all_stats[rock_name][ratio_name]["std"]
            color = COLORS[rock_name]
            shifted_x = x + ROCK_X_OFFSETS[rock_name]

            ax.errorbar(
                shifted_x,
                median,
                yerr=std,
                fmt="o-",
                color=color,
                capsize=4,
                capthick=1.5,
                elinewidth=1.5,
                linewidth=1.8,
                markersize=8,
                markeredgecolor="black",
                markeredgewidth=0.6,
                label=f"{rock_name} ±1σ",
                zorder=3,
            )
            lower_candidates.extend((median - std).tolist())
            upper_candidates.extend((median + std).tolist())

        finite_lower = np.asarray(lower_candidates, dtype=float)
        finite_upper = np.asarray(upper_candidates, dtype=float)
        finite_lower = finite_lower[np.isfinite(finite_lower)]
        finite_upper = finite_upper[np.isfinite(finite_upper)]
        if len(finite_lower) and len(finite_upper):
            y_low = float(finite_lower.min())
            y_high = float(finite_upper.max())
            y_span = y_high - y_low
        else:
            y_span = 1.0
        text_offset = 0.025 * (y_span if y_span > 0 else 1.0)

        # Annotate median values above their data points, matching the reference figure.
        for rock_name in ROCK_ORDER:
            median = all_stats[rock_name][ratio_name]["median"]
            shifted_x = x + ROCK_X_OFFSETS[rock_name]
            color = COLORS[rock_name]
            for xi, yi in zip(shifted_x, median):
                if np.isfinite(yi):
                    ax.text(
                        xi,
                        yi + text_offset,
                        f"{yi:.3g}",
                        ha="center",
                        va="bottom",
                        fontsize=6,
                        color=color,
                        bbox={
                            "facecolor": "white",
                            "alpha": 0.70,
                            "edgecolor": "none",
                            "pad": 1,
                        },
                        zorder=4,
                    )

        ax.legend(loc="best", fontsize=9)

    fig.suptitle(
        "Halogen ratios with mantle mixing "
        "(Serp 50%, Meta 5%, Endmember 45%; Mantle 10%–90%)",
        fontsize=16,
        y=1.03,
    )
    fig.tight_layout()
    return fig


# =============================================================================
# 6. Main program: generate two figures in one run
# =============================================================================

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data_by_rock = load_all_data()

    # Figure 1: 50% mantle mixing in six panels arranged in three rows and two columns.
    density_figure = create_density_figure(data_by_rock)
    density_figure.savefig(DENSITY_FIGURE_FILE, format="svg", bbox_inches="tight")

    # Figure 2: 10%–90% mantle mixing with central 95%, median and ±1σ.
    ratio_stats, summary = compute_ratio_statistics(data_by_rock)
    ratio_figure = create_ratio_figure(ratio_stats)
    ratio_figure.savefig(RATIO_FIGURE_FILE, format="svg", bbox_inches="tight")

    if EXPORT_STATS_CSV:
        summary.to_csv(STATS_FILE, index=False, float_format="%.8g")

    print("Run completed.")
    print(f"Figure 1 (three-by-two density plot): {DENSITY_FIGURE_FILE}")
    print(f"Figure 2 (10%–90% ratio plot): {RATIO_FIGURE_FILE}")
    if EXPORT_STATS_CSV:
        print(f"Statistical table: {STATS_FILE}")

    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close(density_figure)
        plt.close(ratio_figure)


if __name__ == "__main__":
    main()
