import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# ================== 1. Read data ==================
df = pd.read_csv('AOC-BS-HP-UHP Halogens data.csv')
print("Data preview:")
print(df.head())
print("Column names:", df.columns.tolist())
print("Unique rock types:", df['Rock_Type'].unique())

# Unit conversion: ppb → ppm
if 'Br_ppb' in df.columns:
    df['Br_ppm'] = df['Br_ppb'] / 1000
    df.drop('Br_ppb', axis=1, inplace=True)
if 'I_ppb' in df.columns:
    df['I_ppm'] = df['I_ppb'] / 1000
    df.drop('I_ppb', axis=1, inplace=True)

elements = ['F_ppm', 'Cl_ppm', 'Br_ppm', 'I_ppm']
element_labels = ['F', 'Cl', 'Br', 'I']

# Filter out rows in which any halogen concentration is ≤0 to avoid division by zero
df = df[(df[elements] > 0).all(axis=1)].copy()
print(f"Number of samples after retaining positive values: {len(df)}")

# Separate data by rock type
aoc_data = df[df['Rock_Type'] == 'AOC'][elements].values
bs_data = df[df['Rock_Type'] == 'Blueschist'][elements].values
hp_data = df[df['Rock_Type'] == 'HP eclogite'][elements].values
uhp_data = df[df['Rock_Type'] == 'UHP eclogite'][elements].values

# Print diagnostic sample counts and concentration ranges for each rock type
for rock, data in zip(['AOC','Blueschist','HP eclogite','UHP eclogite'], 
                      [aoc_data, bs_data, hp_data, uhp_data]):
    if len(data) > 0:
        print(f"{rock}: n={len(data)}, min={data.min(axis=0)}, max={data.max(axis=0)}")
    else:
        print(f"Warning: no valid data for {rock}!")

# ================== 2. Monte Carlo simulation (sequential stage-wise depletion) ==================
def mc_stage_efficiencies(aoc, bs, hp, uhp, n_iter=10000):
    n_aoc, n_bs = len(aoc), len(bs)
    n_hp, n_uhp = len(hp), len(uhp)
    
    results = {elem: {'stage1': [], 'stage2': [], 'stage3': []} for elem in elements}
    
    for _ in range(n_iter):
        aoc_boot = aoc[np.random.choice(n_aoc, size=n_aoc, replace=True)]
        bs_boot = bs[np.random.choice(n_bs, size=n_bs, replace=True)]
        hp_boot = hp[np.random.choice(n_hp, size=n_hp, replace=True)]
        uhp_boot = uhp[np.random.choice(n_uhp, size=n_uhp, replace=True)]
        
        mean_aoc = np.mean(aoc_boot, axis=0)
        mean_bs = np.mean(bs_boot, axis=0)
        mean_hp = np.mean(hp_boot, axis=0)
        mean_uhp = np.mean(uhp_boot, axis=0)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            eff1 = np.where(mean_aoc != 0, (mean_aoc - mean_bs) / mean_aoc, np.nan)
            eff2 = np.where(mean_bs != 0, (mean_bs - mean_hp) / mean_bs, np.nan)
            eff3 = np.where(mean_hp != 0, (mean_hp - mean_uhp) / mean_hp, np.nan)
        
        remaining1 = 1.0 - eff1
        remaining2 = remaining1 * (1.0 - eff2)
        
        contrib1 = eff1 * 100.0
        contrib2 = eff2 * remaining1 * 100.0
        contrib3 = eff3 * remaining2 * 100.0
        
        for i, elem in enumerate(elements):
            results[elem]['stage1'].append(contrib1[i] if not np.isnan(contrib1[i]) else np.nan)
            results[elem]['stage2'].append(contrib2[i] if not np.isnan(contrib2[i]) else np.nan)
            results[elem]['stage3'].append(contrib3[i] if not np.isnan(contrib3[i]) else np.nan)
    
    dfs = {}
    for elem in elements:
        df_elem = pd.DataFrame(results[elem])
        df_elem.columns = ['Blueschist (stage1)', 'HP Eclogite (stage2)', 'UHP Eclogite (stage3)']
        dfs[elem] = df_elem
    return dfs

# ================== 3. Run simulation ==================
dfs_efficiencies = mc_stage_efficiencies(aoc_data, bs_data, hp_data, uhp_data, n_iter=10000)

# Print statistical summary
for elem in elements:
    print(f"\n===== {elem} normalized stage-specific release efficiency (%) =====")
    print(dfs_efficiencies[elem].describe(percentiles=[0.025, 0.5, 0.975]))

# ================== Helper functions: IQR outlier filtering + 95% confidence interval ==================
def remove_outliers_iqr(data, factor=1.5):
    if len(data) < 4:
        return data
    q1, q3 = np.percentile(data, 25), np.percentile(data, 75)
    iqr = q3 - q1
    lower, upper = q1 - factor * iqr, q3 + factor * iqr
    return data[(data >= lower) & (data <= upper)]

def robust_median_ci(series):
    arr = series.dropna().values
    if len(arr) < 4:
        return np.nan, np.nan, np.nan
    arr_iqr = remove_outliers_iqr(arr)
    if len(arr_iqr) < 4:
        return np.nan, np.nan, np.nan
    median = np.median(arr_iqr)
    low_ci, high_ci = np.percentile(arr_iqr, [2.5, 97.5])
    return median, low_ci, high_ci

# ================== 4. Visualization (element-specific numbers of bins) ==================
ylim_dict = {
    'F': (-400, 400),
    'Cl': (-40, 100),
    'Br': (-40, 100),
    'I': (-40, 100)
}

# Number of histogram bins for each element
bin_count_dict = {
    'F': 100,
    'Cl': 100,
    'Br': 100,
    'I': 200
}

fig, axes = plt.subplots(2, 2, figsize=(14, 14))
axes = axes.flatten()

colors = [(62/255, 158/255, 179/255),   # stage1
          (241/255, 163/255, 97/255),   # stage2
          (136/255, 176/255, 124/255)]  # stage3
stage_labels = ['Blueschist', 'HP Eclogite', 'UHP Eclogite']

box_width = 0.25
hist_max_width = 0.35
box_offset = -0.2

for i, elem in enumerate(elements):
    ax = axes[i]
    ax.set_box_aspect(1)
    
    df_elem = dfs_efficiencies[elem]
    data_raw = [df_elem[col].dropna().values for col in df_elem.columns]
    
    x_pos = np.arange(1, 4)
    
    # Retrieve the number of bins for this element
    bin_count = bin_count_dict[element_labels[i]]
    
    # ---- Calculate the full data range of this panel to define consistent bins ----
    if len(data_raw) > 0 and any(len(d) > 0 for d in data_raw):
        all_data_sub = np.concatenate([d for d in data_raw if len(d) > 0])
        sub_min = all_data_sub.min()
        sub_max = all_data_sub.max()
        if sub_max > sub_min:
            sub_bin_width = (sub_max - sub_min) / bin_count
            sub_bins = np.arange(sub_min, sub_max + sub_bin_width, sub_bin_width)
        else:
            # Use one narrow bin interval when all values are identical
            sub_bins = np.array([sub_min - 0.5, sub_min + 0.5])
    else:
        sub_bins = np.array([0, 1])  # Fallback
    
    # Plot boxplots (left half)
    bp = ax.boxplot(data_raw, positions=x_pos + box_offset,
                    widths=box_width, patch_artist=True,
                    showfliers=False, zorder=5,
                    boxprops=dict(facecolor='lightgray', edgecolor='black'),
                    whiskerprops=dict(color='black'),
                    capprops=dict(color='black'),
                    medianprops=dict(color='red', linewidth=1.5))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    # Plot horizontal histograms and density curves (right half) using panel-wide consistent bins
    for j, (series, pos) in enumerate(zip(data_raw, x_pos)):
        if len(series) < 2:
            continue
        # Use bins defined from the full data range of the panel
        hist, bin_edges = np.histogram(series, bins=sub_bins)
        if hist.max() > 0:
            hist_norm = hist / hist.max() * hist_max_width
        else:
            hist_norm = np.zeros_like(hist)
        ax.barh(bin_edges[:-1], hist_norm, height=np.diff(bin_edges)[0], left=pos,
                color=colors[j], alpha=0.5, edgecolor='none')
        
        # Kernel density estimate (KDE)
        kde = gaussian_kde(series)
        y_vals = np.linspace(series.min(), series.max(), 200)
        density = kde(y_vals)
        density_norm = density / density.max() * hist_max_width
        ax.plot(pos + density_norm, y_vals, color=colors[j], linewidth=1.5, zorder=10)
    
    # Annotate robust median values
    medians = []
    for series in data_raw:
        if len(series) >= 4:
            med, _, _ = robust_median_ci(pd.Series(series))
            medians.append(med)
        else:
            medians.append(np.nan)
    
    y_low, y_high = ylim_dict[element_labels[i]]
    ax.set_ylim(y_low, y_high)
    ax.axhline(y=0, color='black', linewidth=0.8, linestyle='-')
    ax.grid(axis='y', linestyle='--', alpha=0.6)
    
    y_range = y_high - y_low
    text_offset = y_range * 0.08
    for xi, med, col in zip(x_pos, medians, colors):
        if not np.isnan(med):
            ax.text(xi + box_offset, med + text_offset, f'{med:.2f}%',
                    ha='center', va='bottom', fontsize=9, color=col,
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(stage_labels)
    ax.set_title(f'{element_labels[i]} release contribution (% of AOC)', fontsize=14)
    ax.set_ylabel('Release (+) / Uptake (–) (% of AOC)')

plt.tight_layout(pad=2.0)
plt.show()

# ================== 5. Export results ==================
for elem in elements:
    dfs_efficiencies[elem].to_csv(f'{elem}_stage_contributions.csv', index=False)
