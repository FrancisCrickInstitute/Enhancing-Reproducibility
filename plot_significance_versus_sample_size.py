import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import scipy.stats as stats

plt.rcParams['font.size'] = 16

# Load your data
data = pd.read_csv('./plots/selected_wells_raw_data.csv')  # Replace with your file path

# Define the treatment column and variable of interest
treatment_col = 'Treatment'
variable_of_interest = 'Fascin_Ratio'

# Define sample sizes to test
sample_sizes = [*range(10, 500, 10)]  # Modify as needed

# Initialize dictionaries to store results
mean_values = {treatment: [] for treatment in data[treatment_col].unique()}
kruskal_p_values = []

num_iterations = 100

# Define specific pairs for Dunn's test
dunn_pairs = [('Untreated', 'DMSO'), ('DMSO', 'SN0212398523'), ('SN0212398523', 'Leptomycin b')]

# Modify the dictionary initialization to store multiple p-values per sample size
dunn_p_values = {pair: [[] for _ in range(len(sample_sizes))] for pair in dunn_pairs}

# Initialize dictionaries to store multiple mean values per sample size for each treatment
mean_values = {treatment: [[] for _ in range(len(sample_sizes))] for treatment in data[treatment_col].unique()}

for sample_size_index, sample_size in enumerate(sample_sizes):
    for _ in range(num_iterations):
        combined_data = pd.DataFrame()

        for treatment in data[treatment_col].unique():
            subsample = data[data[treatment_col] == treatment].sample(n=sample_size, replace=False)
            mean = subsample[variable_of_interest].mean()
            mean_values[treatment][sample_size_index].append(mean)
            combined_data = pd.concat([combined_data, subsample])

        # Perform Kruskal-Wallis test
        _, p_value = stats.kruskal(
            *(combined_data[combined_data[treatment_col] == t][variable_of_interest] for t in
              combined_data[treatment_col].unique()))

        # Perform Dunn's test if Kruskal-Wallis test is significant
        if p_value < 0.05:
            dunn_result = sp.posthoc_dunn(combined_data, val_col=variable_of_interest, group_col=treatment_col)
            for pair in dunn_pairs:
                dunn_p_values[pair][sample_size_index].append(dunn_result.loc[pair[0], pair[1]])
        else:
            for pair in dunn_pairs:
                dunn_p_values[pair][sample_size_index].append(np.nan)

# Calculate the minimum and maximum for the p-values
dunn_p_means = {pair: np.nanmean(dunn_p_values[pair], axis=1) for pair in dunn_pairs}
dunn_p_25th = {pair: np.nanpercentile(dunn_p_values[pair], 25, axis=1) for pair in dunn_pairs}
dunn_p_75th = {pair: np.nanpercentile(dunn_p_values[pair], 75, axis=1) for pair in dunn_pairs}

# Plotting the Dunn's test p-values with uncertainty ranges
plt.figure(figsize=(14, 10))
for pair in dunn_pairs:
    mean_p_values = dunn_p_means[pair]
    min_p_values = dunn_p_25th[pair]
    max_p_values = dunn_p_75th[pair]
    plt.plot(sample_sizes, mean_p_values, label=f'{pair[0]} vs {pair[1]}')
    plt.fill_between(sample_sizes, min_p_values, max_p_values, alpha=0.2)

plt.xlabel('Sample Size')
plt.ylabel('Dunn Test P-Value')
plt.axhline(y=0.05, color='red', linestyle='dotted', label='p = 0.05')
plt.legend()
plt.show()

# Calculate the mean, minimum, and maximum for the mean values
mean_values_mean = {treatment: np.nanmean(mean_values[treatment], axis=1) for treatment in data[treatment_col].unique()}
mean_values_25th = {treatment: np.nanpercentile(mean_values[treatment], 25, axis=1) for treatment in data[treatment_col].unique()}
mean_values_75th = {treatment: np.nanpercentile(mean_values[treatment], 75, axis=1) for treatment in data[treatment_col].unique()}

# Plotting the mean Fascin_Ratio for each treatment with uncertainty ranges
plt.figure(figsize=(14, 10))
for treatment in data[treatment_col].unique():
    plt.plot(sample_sizes, mean_values_mean[treatment], label=treatment)
    plt.fill_between(sample_sizes, mean_values_25th[treatment], mean_values_75th[treatment], alpha=0.2)

plt.xlabel('Sample Size')
plt.ylabel('Mean Fascin_Ratio')
plt.legend()
plt.show()
