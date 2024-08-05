import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scikit_posthocs as sp
import scipy.stats as stats
import seaborn as sns
from scipy.optimize import curve_fit


# Define the decaying exponential function
def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c


def normalize_well_format(well):
    match = re.match(r"([A-Za-z])([0-9]+)", well)
    return f"{match[1]}{int(match[2]):02d}" if match else well


def load_and_prepare_data(file_path, plate_number):
    df = (pd.read_csv(file_path).query('Plate == @plate_number').assign(
        **{
            'Control Type': lambda x: x['Control Type'].fillna('Treated').replace('', 'Treated'),
            'Well': lambda x: x['Well'].apply(normalize_well_format)
        }
    ))
    return df


def download_csv(file_path, url):
    if not os.path.exists(file_path):
        response = requests.get(url)
        if response.status_code == 200:
            with open(file_path, 'wb') as file:
                file.write(response.content)
            print(f"File downloaded successfully: {file_path}")
        else:
            print(f"Failed to download the file. Status code: {response.status_code}")
    else:
        print("File already exists.")


def prepare_data(nuc_data, cyto_data, image_data, treatments, treatments_to_compounds, compounds, selected_wells):
    # Rename columns
    nuc_data = nuc_data.rename(columns=lambda x: 'Nuclear_' + x if 'Intensity' in x else x)
    cyto_data = cyto_data.rename(columns=lambda x: 'Cyto_' + x if 'Intensity' in x else x)

    # Merge data on ImageNumber and ObjectNumber
    combined_data = nuc_data.merge(cyto_data, on=['ImageNumber', 'ObjectNumber'], how='left')
    combined_data = combined_data.merge(image_data, on='ImageNumber', how='left')

    # Calculate ratios for Fascin and NuclearActin
    for compartment in ['Fascin', 'NuclearActin']:
        combined_data[f'{compartment}_Ratio'] = (
                combined_data[f'Nuclear_Intensity_MeanIntensity_{compartment}'] /
                (combined_data[f'Cyto_Intensity_MeanIntensity_{compartment}'] +
                 combined_data[f'Nuclear_Intensity_MeanIntensity_{compartment}'])
        )

    # Extract well information and map treatments
    combined_data['Well'] = combined_data['FileName_DNA'].str.extract(r'_(.*?)_')[0]
    combined_data = map_wells_to_treatments(combined_data, treatments, treatments_to_compounds, compounds)

    # Filter by selected wells if specified
    if selected_wells:
        combined_data = combined_data[combined_data['Well'].isin(selected_wells)]

    return combined_data.sort_values(by=['Treatment', 'Well'])


def map_wells_to_treatments(data, treatments, treatments_to_compounds, compounds):
    # Map wells to treatment names and then to compound names
    data['Treatment'] = data['Well'].map(treatments).map(treatments_to_compounds)

    # Handle cases where the treatment is 'Treated' differently
    treated_mask = data['Treatment'] == 'Treated'

    # Update the 'Treatment' column for treated data
    data.loc[treated_mask, 'Treatment'] = data.loc[treated_mask, 'Well'].map(compounds)

    # Fill any missing values with 'Unknown' or another appropriate default
    data['Treatment'] = data['Treatment'].fillna('Unknown')

    return data


def generate_swarmplot(fig_width, fig_height, plot_rows, plot_cols, plot_order, n_samples, sample_size, data,
                       color_dict, treatment_col, variable_of_interest, y_label):
    """
    Generates and saves swarm plots for the variable of interest across different treatments.

    Parameters:
    - fig_width, fig_height: Dimensions of the figure.
    - plot_rows, plot_cols: Number of rows and columns in the subplot grid.
    - plot_order: Order in which treatments are displayed on the x-axis.
    - n_samples: Number of sample groups to plot.
    - sample_size: Number of samples to take per treatment (if > 0).
    - data: DataFrame containing the data.
    - color_dict: Dictionary mapping treatments to colors for the plot.
    - treatment_col: Column name indicating the treatment type in the data.
    - variable_of_interest: The dependent variable to be plotted.
    - y_label: Label for the y-axis.
    """

    plt.figure(figsize=(fig_width, fig_height))

    # Sample the data if sample_size > 0
    if sample_size > 0:
        sampled_data = pd.concat([
            data[data[treatment_col] == 'Untreated'].sample(n=sample_size, replace=False),
            data[data[treatment_col] == 'DMSO'].sample(n=sample_size, replace=False),
            data[data[treatment_col] == 'SN0212398523'].sample(n=sample_size, replace=False),
            data[data[treatment_col] == 'Leptomycin b'].sample(n=sample_size, replace=False)
        ])

        # Save sampled data if needed (optional)
        for sample_index in range(n_samples):
            sampled_data.to_csv(f'./outputs/data/sampled_data_{sample_size:03}_{sample_index:01}.csv', index=False)
    else:
        sampled_data = data

    # Plot the data
    for sample_index in range(n_samples):
        ax = plt.subplot(plot_rows, plot_cols, sample_index + 1)
        sns.boxplot(x=treatment_col, y=variable_of_interest, data=sampled_data, order=plot_order, color='white',
                    showfliers=False, ax=ax, linewidth=2)
        sns.swarmplot(x=treatment_col, y=variable_of_interest, data=sampled_data, order=plot_order, palette=color_dict,
                      hue=treatment_col, size=5, alpha=0.9, ax=ax)
        ax.set_ylabel(y_label)
        ax.set_xlabel('')
        ax.set_ylim(bottom=0, top=1.0)

    plt.tight_layout()
    plt.savefig(f'./outputs/plots/selected_wells_{sample_size:03}.png', format='png', bbox_inches='tight')
    plt.show()


def generate_table(data):
    # Calculate descriptive statistics for each well
    descriptive_stats = data.groupby('Well')['Fascin_Ratio'].describe(percentiles=[.25, .5, .75])

    # Calculate 95% confidence interval
    descriptive_stats['95% CI lower'] = descriptive_stats['mean'] - 1.96 * (
            descriptive_stats['std'] / np.sqrt(descriptive_stats['count']))
    descriptive_stats['95% CI upper'] = descriptive_stats['mean'] + 1.96 * (
            descriptive_stats['std'] / np.sqrt(descriptive_stats['count']))

    # Calculate inter-quartile range
    descriptive_stats['IQR'] = descriptive_stats['75%'] - descriptive_stats['25%']

    # Reset index to merge with treatment information
    descriptive_stats.reset_index(inplace=True)

    # Merge with treatment information
    descriptive_stats = pd.merge(descriptive_stats, data[['Well', 'Treatment']].drop_duplicates(), on='Well',
                                 how='left')

    plt.figure(num=1, figsize=(28, 8))
    ax = plt.subplot(1, 1, 1)

    ax.axis('tight')  # turns off the axis lines and labels
    ax.axis('off')  # changes x and y axis limits such that all data is shown

    # plotting data
    table = ax.table(cellText=descriptive_stats.values,
                     colLabels=descriptive_stats.columns,
                     loc="center")
    table.set_fontsize(50)
    # table.scale(1, 2)
    plt.show()


def plot_mean_v_sample_size(sample_sizes, num_iterations, data, treatment_col, variable_of_interest, y_label):
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

    # Calculate the mean, minimum, and maximum for the mean values
    mean_values_mean = {treatment: np.nanmean(mean_values[treatment], axis=1) for treatment in
                        data[treatment_col].unique()}
    mean_values_25th = {treatment: np.nanpercentile(mean_values[treatment], 25, axis=1) for treatment in
                        data[treatment_col].unique()}
    mean_values_75th = {treatment: np.nanpercentile(mean_values[treatment], 75, axis=1) for treatment in
                        data[treatment_col].unique()}

    # Plotting the mean Fascin_Ratio for each treatment with uncertainty ranges
    plt.figure(figsize=(14, 10))
    for treatment in data[treatment_col].unique():
        plt.plot(sample_sizes, mean_values_mean[treatment], label=treatment)
        plt.fill_between(sample_sizes, mean_values_25th[treatment], mean_values_75th[treatment], alpha=0.2)

    plt.xlabel('Sample Size')
    plt.ylabel(y_label)
    plt.legend(fontsize=20)
    plt.show()


def plot_effect_size_v_sample_size(sample_sizes, num_iterations, data, treatment_col, variable_of_interest, y_label,
                                   treatments):
    # Initialize dictionaries to store multiple mean values per sample size for each treatment
    mean_values = {treatment: [[] for _ in range(len(sample_sizes))] for treatment in data[treatment_col].unique()}
    for sample_size_index, sample_size in enumerate(sample_sizes):
        for _ in range(num_iterations):
            combined_data = pd.DataFrame()
            for treatment in treatments:
                subsample = data[data[treatment_col] == treatment].sample(n=sample_size, replace=False)
                control_subsample = data[data[treatment_col] == 'Untreated'].sample(n=sample_size, replace=False)
                mean = subsample[variable_of_interest].mean() - control_subsample[variable_of_interest].mean()
                mean_values[treatment][sample_size_index].append(mean)
                combined_data = pd.concat([combined_data, subsample])

    # Calculate the mean, minimum, and maximum for the mean values
    median_values_mean = {treatment: np.nanmedian(mean_values[treatment], axis=1) for treatment in
                          data[treatment_col].unique()}
    mean_values_25th = {treatment: np.nanpercentile(mean_values[treatment], 25, axis=1) for treatment in
                        data[treatment_col].unique()}
    mean_values_75th = {treatment: np.nanpercentile(mean_values[treatment], 75, axis=1) for treatment in
                        data[treatment_col].unique()}

    # Plotting the mean Fascin_Ratio for each treatment with uncertainty ranges

    for t in range(len(treatments)):
        plt.figure(figsize=(15, 10))
        for treatment in treatments[:t + 1]:
            plt.plot(sample_sizes, median_values_mean[treatment], label=treatment)
            plt.fill_between(sample_sizes, mean_values_25th[treatment], mean_values_75th[treatment], alpha=0.2)

        plt.xlabel('Sample Size')
        plt.ylabel(y_label)
        plt.legend(fontsize=20)
        plt.axhline(y=0.0, color='black', linestyle='dotted')
        plt.show()


def plot_iqr_v_sample_size(sample_sizes, num_iterations, data, treatment_col, variable_of_interest, y_label):
    # Initialize dictionaries to store multiple mean values per sample size for each treatment
    mean_values = {treatment: [[] for _ in range(len(sample_sizes))] for treatment in data[treatment_col].unique()}
    for sample_size_index, sample_size in enumerate(sample_sizes):
        for _ in range(num_iterations):
            combined_data = pd.DataFrame()
            for treatment in data[treatment_col].unique():
                subsample = data[data[treatment_col] == treatment].sample(n=sample_size, replace=False)
                q1, q3 = np.percentile(subsample[variable_of_interest], [25, 75])
                iqr = q3 - q1
                mean_values[treatment][sample_size_index].append(iqr)
                combined_data = pd.concat([combined_data, subsample])

    # Calculate the mean, minimum, and maximum for the mean values
    iqr_values_min = {treatment: np.nanmin(mean_values[treatment], axis=1) for treatment in
                      data[treatment_col].unique()}
    iqr_values_max = {treatment: np.nanmax(mean_values[treatment], axis=1) for treatment in
                      data[treatment_col].unique()}

    # Plotting
    plt.figure(figsize=(14, 10))
    for treatment in data[treatment_col].unique():
        diff_iqr = iqr_values_max[treatment] - iqr_values_min[treatment]
        plt.scatter(sample_sizes, diff_iqr, label=treatment, alpha=0.5)
        initial_guesses = [1, 0.01, np.median(diff_iqr)]
        # Fit the decaying exponential function
        params, _ = curve_fit(exp_decay, sample_sizes, diff_iqr, p0=initial_guesses, maxfev=5000)

        # Generate the fitted curve
        fitted_curve = exp_decay(np.array(sample_sizes), *params)
        plt.plot(sample_sizes, fitted_curve, label=f"{treatment} exp fit", linestyle='--')

    plt.xlabel('Sample Size')
    plt.ylabel(y_label)
    plt.legend(fontsize=20)
    plt.show()


def plot_cumulative_histogram_samples(data, variable_of_interest, treatment_col, treatment):
    total_samples = []
    max_samples = 500
    step = 10

    subsample = data[data[treatment_col] == treatment]

    median_values = []
    mean_values = []
    std_values = []
    iqr_values = []
    sample_sizes = []

    for sample_size in range(step, max_samples + 1, step):
        # Determine the number of new samples to add
        new_samples_count = sample_size - len(total_samples)

        # Ensure we don't sample more than what's available in the dataframe
        remaining_samples = subsample[~subsample.index.isin(total_samples)].shape[0]
        new_samples_count = min(new_samples_count, remaining_samples)

        # Sample additional data and add it to the total_samples list
        if new_samples_count > 0:
            new_samples = subsample[~subsample.index.isin(total_samples)].sample(n=new_samples_count,
                                                                                 replace=False).index.tolist()
            total_samples.extend(new_samples)

        # Extract the data for the current total samples
        sample_data = subsample.loc[total_samples, variable_of_interest]

        median = sample_data.median()
        mean = sample_data.mean()
        std = sample_data.std()
        q1 = sample_data.quantile(0.25)
        q3 = sample_data.quantile(0.75)
        iqr = q3 - q1

        median_values.append(median)
        mean_values.append(mean)
        std_values.append(std)
        iqr_values.append(iqr)
        sample_sizes.append(sample_size)

        if sample_size in (20, 50, 100, 200, 300, 500):
            # Plot histogram
            plt.figure(figsize=(14, 10))
            n, bins, patches = plt.hist(sample_data, bins=50, alpha=0.75, density=True)
            plt.axvline(x=median, color='r', linestyle='--', label='Median')
            plt.axvline(x=q1, color='g', linestyle='-', label='Q1')
            plt.axvline(x=q3, color='b', linestyle='-', label='Q3')
            # Calculate the density
            bin_maxes = np.maximum.reduceat(n, np.digitize([q1, q3], bins[:-1]) - 1)
            max_density = max(bin_maxes)

            # Shade the IQR region
            plt.fill_betweenx(np.arange(0, max_density, 0.01), q1, q3, color='grey', alpha=0.3, label='IQR')

            plt.title(f'{len(total_samples)} random samples from {treatment}')
            plt.xlabel(variable_of_interest)
            plt.ylabel('Frequency')
            plt.ylim(bottom=0, top=20)
            plt.xlim(left=0, right=1)
            plt.grid(True)
            plt.show()

        # print(
        #     f'Median: {np.median(sample_data)} IQR: {np.percentile(sample_data, 75) - np.percentile(sample_data, 25)}')

        # Break the loop if we have included all available samples
        if remaining_samples <= new_samples_count:
            break

    plt.figure(figsize=(14, 10))
    ax1 = plt.gca()
    ax1.scatter(sample_sizes, mean_values, label='_Mean', alpha=0.5, color='blue')
    ax1.scatter(sample_sizes, median_values, label='_Median', alpha=0.5, color='orange')
    ax1.plot(sample_sizes, mean_values, label='Mean', color='blue')
    ax1.plot(sample_sizes, median_values, label='Median', color='orange')
    ax1.set_ylabel('Mean, Median % Nuclear Fascin')
    ax1.set_xlabel('Sample Size')
    ax2 = ax1.twinx()
    ax2.scatter(sample_sizes, std_values, label='_Standard Deviation', alpha=0.5, color='gray')
    ax2.scatter(sample_sizes, iqr_values, label='_IQR', alpha=0.5, color='purple')
    ax2.plot(sample_sizes, std_values, label='Standard Deviation', color='gray')
    ax2.plot(sample_sizes, iqr_values, label='IQR', color='purple')
    ax2.set_ylabel('SD, IQR of % Nuclear Fascin')
    # Create a single legend
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)
    plt.show()


# Assuming you have a dataframe 'df' loaded with the column 'your_column_name',
# you would call the function like this:
# plot_cumulative_histogram_samples(df, 'your_column_name')


def plot_p_v_sample_size(sample_sizes, num_iterations, data, treatment_col, variable_of_interest, dunn_pairs):
    # Modify the dictionary initialization to store multiple p-values per sample size
    dunn_p_values = {pair: [[] for _ in range(len(sample_sizes))] for pair in dunn_pairs}
    for sample_size_index, sample_size in enumerate(sample_sizes):
        for _ in range(num_iterations):
            combined_data = pd.DataFrame()

            for treatment in data[treatment_col].unique():
                subsample = data[data[treatment_col] == treatment].sample(n=sample_size, replace=False)
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
    plt.legend(fontsize=20)
    plt.show()


def generate_swarmplot_of_well_means(fig_width, fig_height, plot_order, treatments, data, color_dict, treatment_col,
                                     variable_of_interest, y_label, dunn_pairs, sample_size):
    # treatments = data[treatment_col].unique()
    mean_data = pd.DataFrame()
    dunn_p_values = {pair: [] for pair in dunn_pairs}
    plt.figure(num=1, figsize=(fig_width, fig_height))

    for t in treatments:
        tdata = data[data[treatment_col] == t]
        wells = tdata['Well'].unique()
        # if len(wells) > 3:
        #     wells = np.random.choice(wells, 3)
        for w in wells:
            if sample_size > 0:
                sdata = tdata[tdata['Well'] == w].sample(n=sample_size, replace=False)
            else:
                sdata = tdata[tdata['Well'] == w]
            smean = np.mean(sdata[variable_of_interest])
            print(f'Treatment: {t} Well: {w} Mean: {smean}')
            mean_data = pd.concat([mean_data,
                                   pd.DataFrame(data=[[t, smean]],
                                                columns=['TREATMENT', 'MEAN_FASCIN_RATIO'])])

    # ax = plt.subplot(1, 1, 1)
    sns.boxplot(x='TREATMENT', y='MEAN_FASCIN_RATIO', data=mean_data, order=plot_order, color='white',
                showfliers=False, linecolor='black', linewidth=2)
    sns.swarmplot(x='TREATMENT', y='MEAN_FASCIN_RATIO', data=mean_data, order=plot_order, palette=color_dict,
                  hue='TREATMENT', size=18, alpha=0.8)
    plt.ylabel(y_label)
    plt.xlabel('')
    plt.ylim(bottom=0.42, top=0.60)
    plt.title(f'Each dot represents the mean of {sample_size} cells')

    # _, p_value = stats.kruskal(
    #     *(mean_data[mean_data['TREATMENT'] == t]['MEAN_FASCIN_RATIO'] for t in
    #       mean_data['TREATMENT'].unique()))
    #
    # dunn_result = sp.posthoc_dunn(mean_data, val_col='MEAN_FASCIN_RATIO', group_col='TREATMENT')
    # for pair in dunn_pairs:
    #     dunn_p_values[pair].append(dunn_result.loc[pair[0], pair[1]])
    #
    # y, h, col = mean_data['MEAN_FASCIN_RATIO'].max() + 0.005, 0.005, 'black'
    #
    # ymax = []
    # for t in range(len(mean_data['TREATMENT'].unique()) - 1):
    #     ymax.append(0)
    #
    # for pair in dunn_pairs:
    #     x1, x2 = pair
    #     x1 = [label.get_text() for label in ax.get_xticklabels()].index(x1)
    #     x2 = [label.get_text() for label in ax.get_xticklabels()].index(x2)
    #
    #     y = mean_data[mean_data['TREATMENT'].isin(pair)].loc[:, 'MEAN_FASCIN_RATIO'].max() + 0.02
    #
    #     for x in range(min(x1, x2), max(x1, x2)):
    #         if y <= ymax[x] + 0.01:
    #             y = ymax[x] + 0.01
    #         ymax[x] = y
    #
    #     if x1 < x2:
    #         x1 += 0.02
    #         x2 -= 0.02
    #     else:
    #         x1 -= 0.02
    #         x2 += 0.02
    #
    #     # Draw line
    #     plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
    #
    #     str_p_value = f'p = {dunn_p_values[pair][0]:.3f}'
    #
    #     # Annotate line with p-value
    #     plt.text((x1 + x2) * .5, y + h, str_p_value, ha='center', va='bottom', color=col, fontsize=20)

    plt.show()
