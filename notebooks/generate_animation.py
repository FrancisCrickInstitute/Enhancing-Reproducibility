from utility_functions import generate_superplot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams['font.size'] = 26
plt.rcParams['axes.linewidth'] = 2

# Set seed for reproducibility
np.random.seed(42)

# Parameters
n_experiments = 3
treatments = {
    'Untreated Cells': {'mean': 10.0, 'std': 2.5},
    'Negative Control': {'mean': 10.0, 'std': 2.5},
    'Positive Control': {'mean': 11.0, 'std': 3.0},
    'Test Treatment': {'mean': 10.5, 'std': 3.0},
}
n_samples_per_group = 100  # observations per treatment per experiment

# Build the DataFrame
rows = []
for experiment in range(1, n_experiments + 1):
    for treatment_name, params in treatments.items():
        values = np.random.normal(
            loc=params['mean'],
            scale=params['std'],
            size=n_samples_per_group
        )
        for value in values:
            rows.append({
                'Experiment': experiment,
                'Treatment': treatment_name,
                'Metric of Interest': value
            })

df = pd.DataFrame(rows)

# Quick look
print(df.head())
print("\nShape:", df.shape)
print("\nSummary by treatment:")
print(df.groupby('Treatment')['Metric of Interest'].agg(['mean', 'std']))

# --- Precompute running averages (mean across experiments) for every sample size ---
treatment_names = list(treatments.keys())
sample_axis = np.arange(1, n_samples_per_group + 1)
running_averages = pd.DataFrame(index=sample_axis, columns=treatment_names, dtype=float)

for treatment_name in treatment_names:
    per_experiment_running_means = []
    for experiment in range(1, n_experiments + 1):
        values = df[(df['Treatment'] == treatment_name) &
                    (df['Experiment'] == experiment)]['Metric of Interest'].values
        # Running mean for this experiment at each sample size
        per_experiment_running_means.append(np.cumsum(values) / sample_axis)
    # Average across the three experimental replicates
    running_averages[treatment_name] = np.mean(per_experiment_running_means, axis=0)

# Distinct palette for the 4 treatments (separate from the experiment colors in the superplot)
treatment_palette = {
    'Untreated Cells':  '#4C72B0',
    'Negative Control': '#55A868',
    'Positive Control': '#C44E52',
    'Test Treatment':   '#8172B2',
}

point_size = 3
samples_sizes = range(n_samples_per_group)

for s in samples_sizes:
    print(f'Running sample size {s}')
    s_string = str(s).zfill(3)

    # --- Existing superplot frame ---
    generate_superplot(treatments.keys(), treatments.keys(),
                       df, 'Treatment', 'Experiment', 'Metric of Interest',
                       'Metric of Interest', sample_size=s, point_size=point_size,
                       filename=f'../plots/plot_{s_string}.png')

    # --- New: replicate-averaged metric vs sample size ---
    plt.figure(figsize=(10, 8))
    ax = plt.subplot(1, 1, 1)

    if s >= 1:
        x = sample_axis[:s]
        for treatment_name in treatment_names:
            y = running_averages[treatment_name].iloc[:s].values
            plt.plot(x, y, linewidth=3,
                     label=treatment_name,
                     color=treatment_palette[treatment_name])

    plt.xlim(0, n_samples_per_group)
    plt.ylim(7.5, 15.0)  # matches superplot y-range; tighten (e.g. 7–14) to emphasise separation
    plt.xlabel('Number of cells per treatment, per experiment')
    plt.ylabel('Mean Metric of Interest')
    plt.title(' ')
    plt.legend(loc='upper right', fontsize=16, frameon=False)

    plt.tight_layout()
    plt.savefig(f'../plots/avg_{s_string}.png')
    plt.close()