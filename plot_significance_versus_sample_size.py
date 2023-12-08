import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import scikit_posthocs as sp
import numpy as np

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

# Define specific pairs for Dunn's test
dunn_pairs = [('Untreated', 'DMSO'), ('DMSO', 'SN0212398523'), ('SN0212398523', 'Leptomycin b')]
dunn_p_values = {pair: [] for pair in dunn_pairs}

for sample_size in sample_sizes:
    combined_data = pd.DataFrame()

    for treatment in data[treatment_col].unique():
        # Sampling without replacement for each treatment
        subsample = data[data[treatment_col] == treatment].sample(n=sample_size, replace=False)
        mean = subsample[variable_of_interest].mean()
        mean_values[treatment].append(mean)
        combined_data = pd.concat([combined_data, subsample])

    # Perform Kruskal-Wallis test
    _, p_value = stats.kruskal(
        *(combined_data[combined_data[treatment_col] == t][variable_of_interest] for t in combined_data[treatment_col].unique()))
    kruskal_p_values.append(p_value)

    # Perform Dunn's test if Kruskal-Wallis test is significant
    if p_value < 0.05:
        dunn_result = sp.posthoc_dunn(combined_data, val_col=variable_of_interest, group_col=treatment_col)
        for pair in dunn_pairs:
            dunn_p_values[pair].append(dunn_result.loc[pair[0], pair[1]])
    else:
        for pair in dunn_pairs:
            dunn_p_values[pair].append(np.nan)  # Append NaN if Kruskal-Wallis is not significant

# Plotting the mean Fascin_Ratio for each treatment
plt.figure(figsize=(14, 10))
for treatment, means in mean_values.items():
    plt.plot(sample_sizes, means, label=treatment)
plt.xlabel('Sample Size')
plt.ylabel('Mean_Nuclear_Fascin_Intensity /\n (Mean_Nuclear_Fascin_Intensity + Mean_Cytoplasmic_Fascin_Intensity)')
#plt.title('Mean Fascin_Ratio by Treatment and Sample Size')
plt.legend()
plt.show()

# Plotting the Kruskal-Wallis p-values
plt.figure(figsize=(14, 10))
plt.plot(sample_sizes, kruskal_p_values, marker='o')
plt.xlabel('Sample Size')
plt.ylabel('Kruskal-Wallis P-Value')
#plt.title('Kruskal-Wallis P-Values for Different Sample Sizes')
plt.yscale('log')  # Log scale to better visualize p-values
plt.show()

# Plotting the Dunn's test p-values for specific comparisons
plt.figure(figsize=(14, 10))
for pair, p_values in dunn_p_values.items():
    plt.plot(sample_sizes, p_values, label=f'{pair[0]} vs {pair[1]}')
plt.xlabel('Sample Size')
plt.ylabel('Dunn Test P-Value')
#plt.title('Dunn Test P-Values for Specific Comparisons Across Sample Sizes')
plt.legend()
#plt.yscale('log')  # Log scale for better visualization
plt.axhline(y=0.05, color='red', linestyle='dotted', label='p = 0.05')  # Horizontal line at p = 0.05
plt.show()
