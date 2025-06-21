import itertools
import re
from typing import List, Any, Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sph
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
from tqdm.notebook import tqdm


def generate_pairs(input_list: List[Any]) -> List[Tuple[Any, Any]]:
    """Generate all possible unique pairs from the elements of the input list.

    Parameters:
    input_list (List[Any]): A list of elements from which pairs are to be generated.

    Returns:
    List[Tuple[Any, Any]]: A list of tuples, where each tuple is a unique pair of elements.
    """
    return list(itertools.combinations(input_list, 2))


def exp_decay(x: float, a: float, b: float, c: float) -> float:
    """
    Compute the value of an exponential decay function f(x) = a * e^(-b * x) + c.

    Parameters:
    x (float): The point at which to evaluate the function.
    a (float): Initial amplitude.
    b (float): Decay rate.
    c (float): Asymptotic value.

    Returns:
    float: The value of the exponential decay function at point x.
    """
    return a * np.exp(-b * x) + c


def normalize_well_format(well: str) -> str:
    """
    Normalize a well identifier to a standard format with a letter followed by a two-digit number.

    Parameters:
    well (str): The well identifier to normalize (e.g., "A1" or "B24").

    Returns:
    str: The normalized well identifier (e.g., "A01" or "B24"). Returns the original string if no match.
    """
    match = re.match(r"([A-Za-z])([0-9]+)", well)
    return f"{match[1]}{int(match[2]):02d}" if match else well


def load_and_prepare_data(file_path: str, plate_number, column: str, fill: str) -> pd.DataFrame:
    """
    Load a CSV file, filter by plate number, and prepare specified columns.

    This function reads a CSV file, filters rows by plate number, fills missing values
    in a specified column, and normalizes well identifiers.

    Parameters:
    file_path (str): Path to the CSV file.
    plate_number (int): Plate number to filter by.
    column (str): Column name where missing values will be filled.
    fill (str): Value to fill missing or empty entries in the specified column.

    Returns:
    pd.DataFrame: A DataFrame filtered by plate number with the specified transformations applied.
    """
    df = (pd.read_csv(file_path)
          .query('Plate == @plate_number')
          .assign(**{
        column: lambda x: x[column].fillna(fill).replace('', fill),
        'Well': lambda x: x['Well'].apply(normalize_well_format)
    }))
    return df


def prepare_data(nuc_data, cyto_data, image_data, image_indices, treatments, treatments_to_compounds,
                 compounds, selected_wells, proteins_of_interest):
    """
    Prepare and merge data from nuclear, cytoplasmic, and image datasets, calculating ratios and mapping treatments.

    Parameters:
    nuc_data (pd.DataFrame): Nuclear data.
    cyto_data (pd.DataFrame): Cytoplasmic data.
    image_data (pd.DataFrame): Image data.
    image_indices (pd.DataFrame): DataFrame containing image indices and well names.
    treatments (Dict): Dictionary mapping well identifiers to treatment names.
    treatments_to_compounds (Dict): Dictionary mapping treatment names to compound names.
    compounds (Dict): Dictionary mapping well identifiers to compound names for treated cases.
    selected_wells (List): List of wells to filter by.
    proteins_of_interest (List): List of protein names for which ratios are calculated.

    Returns:
    pd.DataFrame: Combined and processed DataFrame.
    """
    # Rename columns
    nuc_data = nuc_data.rename(columns=lambda x: 'Nuclear_' + x if 'Intensity' in x else x)
    cyto_data = cyto_data.rename(columns=lambda x: 'Cyto_' + x if 'Intensity' in x else x)

    # Merge data
    combined_data = nuc_data.merge(cyto_data, on=['ImageNumber', 'ObjectNumber'], how='left')
    combined_data = combined_data.merge(image_data, on='ImageNumber', how='left')

    # Calculate ratios
    for compartment in proteins_of_interest:
        nuclear_intensity = f'Nuclear_Intensity_MeanIntensity_{compartment}'
        cyto_intensity = f'Cyto_Intensity_MeanIntensity_{compartment}'
        combined_data[f'{compartment}_Ratio'] = (
                combined_data[nuclear_intensity] / (combined_data[cyto_intensity] + combined_data[nuclear_intensity])
        )

    # Map well names
    if image_indices is not None:
        filename_to_well = dict(zip(image_indices['sourcefilename'], image_indices['WellName']))
        combined_data['Well'] = combined_data['FileName_Hoechst'].map(filename_to_well)
    else:
        combined_data['Well'] = combined_data['FileName_DNA'].str.extract(r'_(.*?)_')[0]

    # Normalize well format and map treatments
    combined_data['Well'] = combined_data['Well'].apply(normalize_well_format)
    combined_data = map_wells_to_treatments(combined_data, treatments, treatments_to_compounds, compounds)

    # Filter by selected wells
    if selected_wells:
        combined_data = combined_data[combined_data['Well'].isin(selected_wells)]

    return combined_data.sort_values(by=['Treatment', 'Well'])


def map_wells_to_treatments(data: pd.DataFrame, treatments: Dict, treatments_to_compounds: Dict,
                            compounds: Dict) -> pd.DataFrame:
    """
    Map well identifiers to treatments and compounds, handling special cases and filling missing values.

    Parameters:
    data (pd.DataFrame): DataFrame containing well data.
    treatments (Dict): Dictionary mapping well identifiers to treatment names.
    treatments_to_compounds (Dict): Dictionary mapping treatment names to compound names.
    compounds (Dict): Dictionary mapping well identifiers to compound names for treated cases.

    Returns:
    pd.DataFrame: DataFrame with an added 'Treatment' column, mapping wells to treatments or compounds.
    """
    # Map wells to treatment names and then to compound names
    if treatments_to_compounds:
        data['Treatment'] = data['Well'].map(treatments).map(treatments_to_compounds)
    else:
        data['Treatment'] = data['Well'].map(treatments)

    # Handle cases where the treatment is 'Treated' differently
    treated_mask = data['Treatment'] == 'Treated'

    # Update the 'Treatment' column for treated data
    data.loc[treated_mask, 'Treatment'] = data.loc[treated_mask, 'Well'].map(compounds)

    # Fill any missing values with 'Unknown' or another appropriate default
    data['Treatment'] = data['Treatment'].fillna('Unknown')

    return data


def generate_swarmplot(plot_order, data, color_dict, treatment_col, variable_of_interest, y_label,
                       point_size=2, p_values=False, random_seed=42, fig_width=14, fig_height=10,
                       plot_rows=1, plot_cols=1, n_samples=1, sample_size=-1, filename=None, title=None):
    """
    Generates and saves swarm plots for the variable of interest across different treatments.

    Parameters:
    - plot_order: Order in which treatments are displayed on the x-axis.
    - data: DataFrame containing the data.
    - color_dict: Dictionary mapping treatments to colors for the plot.
    - treatment_col: Column name indicating the treatment type in the data.
    - variable_of_interest: The dependent variable to be plotted.
    - y_label: Label for the y-axis.
    - point_size: Size of points in the swarm plot.
    - p_values: Whether to calculate and display p-values.
    - random_seed: Random seed for reproducibility.
    - fig_width, fig_height: Dimensions of the figure.
    - plot_rows, plot_cols: Number of rows and columns in the subplot grid.
    - n_samples: Number of sample groups to plot.
    - sample_size: Number of samples to take per treatment (if > 0).
    - filename: Filename to save the plot.
    - title: Title of the plot.
    """

    # Initialize the plot with specified dimensions
    plt.figure(figsize=(fig_width, fig_height))
    dunn_pairs = generate_pairs(plot_order)
    dunn_p_values = {pair: [] for pair in dunn_pairs}

    # Sample the data if a sample size is specified
    if sample_size > 0:
        # Sample data for each treatment and concatenate into a single DataFrame
        sampled_data = pd.concat([
            data[data[treatment_col] == treatment].sample(n=sample_size, replace=False, random_state=random_seed)
            for treatment in plot_order
        ])
    else:
        # Use the entire dataset if no sampling is required
        sampled_data = data

    # Generate plots for each sample
    for sample_index in range(n_samples):
        ax = plt.subplot(plot_rows, plot_cols, sample_index + 1)

        # Create a swarm plot for the sampled data
        sns.swarmplot(
            x=treatment_col, y=variable_of_interest, data=sampled_data, order=plot_order,
            palette=color_dict, hue=treatment_col, size=point_size, alpha=0.9, ax=ax, zorder=1
        )

        # Overlay a box plot on the swarm plot
        sns.boxplot(
            x=treatment_col, y=variable_of_interest, data=sampled_data, order=plot_order,
            boxprops=dict(facecolor='none', zorder=2),
            whiskerprops=dict(color="black", linewidth=2, zorder=2),
            capprops=dict(color="black", linewidth=2, zorder=2),
            medianprops=dict(color="black", linewidth=2, zorder=2),
            showfliers=False, ax=ax
        )

        ax.set_ylabel(y_label)
        ax.set_xlabel('')

        if p_values:
            # Perform Kruskal-Wallis test to compare samples
            _, p_value = stats.kruskal(*(
                sampled_data[sampled_data[treatment_col] == t][variable_of_interest]
                for t in sampled_data[treatment_col].unique()
            ))

            # Perform Dunn's test for post-hoc analysis
            dunn_result = sph.posthoc_dunn(sampled_data, val_col=variable_of_interest, group_col=treatment_col)
            for pair in dunn_pairs:
                dunn_p_values[pair].append(dunn_result.loc[pair[0], pair[1]])

            # Prepare to annotate p-values on the plot
            ymax = []
            for t in range(len(sampled_data[treatment_col].unique()) - 1):
                ymax.append(0)

            for pair in dunn_pairs:
                x1, x2 = pair
                # Find the x positions of the pair on the plot
                x1 = [label.get_text() for label in ax.get_xticklabels()].index(x1)
                x2 = [label.get_text() for label in ax.get_xticklabels()].index(x2)

                # Determine the y position for the annotation
                y = sampled_data[sampled_data[treatment_col].isin(pair)].loc[:, variable_of_interest].max() + 0.02

                # Adjust y position to avoid overlap
                for x in range(min(x1, x2), max(x1, x2)):
                    if y <= ymax[x] + 0.05:
                        y = ymax[x] + 0.05
                    ymax[x] = y

                # Adjust x positions for better visualization
                if x1 < x2:
                    x1 += 0.02
                    x2 -= 0.02
                else:
                    x1 -= 0.02
                    x2 += 0.02

                # Draw a line connecting the pair
                plt.plot([x1, x1, x2, x2], [y, y + 0.005, y + 0.005, y], lw=1.5, c='black')

                # Format the p-value string based on its significance
                str_p_value = f'p = {dunn_p_values[pair][0]:.3f}'
                if dunn_p_values[pair][0] < 0.0001:
                    str_p_value = '****'
                elif dunn_p_values[pair][0] < 0.001:
                    str_p_value = '***'
                elif dunn_p_values[pair][0] < 0.01:
                    str_p_value = '**'
                elif dunn_p_values[pair][0] < 0.05:
                    str_p_value = '*'

                # Annotate the plot with the p-value
                plt.text((x1 + x2) * 0.5, y + 0.005, str_p_value, ha='center', va='bottom', color='black', fontsize=20)

    if title is not None:
        plt.title(title)

    plt.tight_layout()

    if filename is not None:
        plt.savefig(f'../plots/{filename}')

    plt.show()


def plot_effect_size_v_sample_size(sample_sizes, num_iterations, data, treatment_col, variable_of_interest, y_label,
                                   treatments, control_name, initial_random_seed=42):
    """
    Plot the effect size of a variable of interest across different sample sizes for each treatment compared to a control.

    Parameters:
    - sample_sizes: List of sample sizes to iterate over.
    - num_iterations: Number of iterations to perform for each sample size.
    - data: DataFrame containing the data.
    - treatment_col: Column name indicating the treatment type in the data.
    - variable_of_interest: The variable for which the effect size is calculated and plotted.
    - y_label: Label for the y-axis of the plot.
    - treatments: List of treatments to analyze.
    - control_name: Name of the control treatment.
    - initial_random_seed: Initial random seed for reproducibility.
    """

    # Initialize dictionaries to store multiple mean values per sample size for each treatment
    mean_values = {treatment: [[] for _ in range(len(sample_sizes))] for treatment in treatments}
    random_seed = initial_random_seed

    # Use tqdm.notebook.tqdm to show progress over all iterations
    total_iterations = len(sample_sizes) * num_iterations
    with tqdm(total=total_iterations, desc="Processing") as pbar:
        for sample_size_index, sample_size in enumerate(sample_sizes):
            for _ in range(num_iterations):
                combined_data = pd.DataFrame()
                for treatment in treatments:
                    # Subsample data for each treatment and control
                    subsample = data[data[treatment_col] == treatment].sample(n=sample_size, replace=False,
                                                                              random_state=random_seed)
                    control_subsample = data[data[treatment_col] == control_name].sample(n=sample_size, replace=False,
                                                                                         random_state=random_seed)

                    # Calculate the effect size as the difference in means divided by the control's standard deviation
                    mean = (subsample[variable_of_interest].mean() - control_subsample[variable_of_interest].mean()) / \
                           control_subsample[variable_of_interest].std()
                    mean_values[treatment][sample_size_index].append(mean)
                    combined_data = pd.concat([combined_data, subsample])
                    random_seed += 1

                # Update the progress bar
                pbar.update(1)

    # Calculate the median, 25th percentile, and 75th percentile for the effect sizes
    median_values_mean = {treatment: np.nanmedian(mean_values[treatment], axis=1) for treatment in treatments}
    mean_values_25th = {treatment: np.nanpercentile(mean_values[treatment], 25, axis=1) for treatment in treatments}
    mean_values_75th = {treatment: np.nanpercentile(mean_values[treatment], 75, axis=1) for treatment in treatments}

    # Plotting the effect sizes for each treatment with uncertainty ranges
    for t in range(len(treatments)):
        plt.figure(figsize=(15, 10))
        for treatment in treatments[:t + 1]:
            # Plot the median effect size for each treatment
            plt.plot(sample_sizes, median_values_mean[treatment], label=treatment)
            # Fill the area between the 25th and 75th percentiles to show uncertainty
            plt.fill_between(sample_sizes, mean_values_25th[treatment], mean_values_75th[treatment], alpha=0.2)

        plt.xlabel('Number of Cells')
        plt.ylabel(y_label)
        plt.legend(fontsize=20)
        plt.axhline(y=0.0, color='black', linestyle='dotted')
        plt.show()


def plot_iqr_v_sample_size(sample_sizes, num_iterations, data, treatment_col, variable_of_interest, y_label,
                           initial_random_seed=42, filename=None):
    """
    Plot the interquartile range (IQR) of a variable of interest across different sample sizes for each treatment.

    Parameters:
    - sample_sizes: List of sample sizes to iterate over.
    - num_iterations: Number of iterations to perform for each sample size.
    - data: DataFrame containing the data.
    - treatment_col: Column name indicating the treatment type in the data.
    - variable_of_interest: The variable for which the IQR is calculated and plotted.
    - y_label: Label for the y-axis of the plot.
    - initial_random_seed: Initial random seed for reproducibility.
    - filename: Filename to save the plot.
    """

    # Initialize dictionaries to store IQR values per sample size for each treatment
    iqr_values = {treatment: [[] for _ in range(len(sample_sizes))] for treatment in data[treatment_col].unique()}
    random_seed = initial_random_seed

    # Use tqdm to show progress over all iterations
    total_iterations = len(sample_sizes) * num_iterations
    with tqdm(total=total_iterations, desc="Processing") as pbar:
        for sample_size_index, sample_size in enumerate(sample_sizes):
            for _ in range(num_iterations):
                combined_data = pd.DataFrame()
                for treatment in data[treatment_col].unique():
                    # Subsample data for each treatment
                    subsample = data[data[treatment_col] == treatment].sample(n=sample_size, replace=False,
                                                                              random_state=random_seed)

                    # Calculate the interquartile range (IQR)
                    q1, q3 = np.percentile(subsample[variable_of_interest], [25, 75])
                    iqr = q3 - q1
                    iqr_values[treatment][sample_size_index].append(iqr)
                    combined_data = pd.concat([combined_data, subsample])

                random_seed += 1
                pbar.update(1)  # Update the progress bar

    # Calculate the minimum and maximum IQR values for each treatment
    iqr_values_min = {treatment: np.nanmin(iqr_values[treatment], axis=1) for treatment in data[treatment_col].unique()}
    iqr_values_max = {treatment: np.nanmax(iqr_values[treatment], axis=1) for treatment in data[treatment_col].unique()}

    # Plotting the IQR differences for each treatment
    plt.figure(figsize=(14, 10))
    for treatment in data[treatment_col].unique():
        # Calculate the difference between max and min IQR values
        diff_iqr = iqr_values_max[treatment] - iqr_values_min[treatment]
        plt.scatter(sample_sizes, diff_iqr, label=treatment, alpha=0.5)

        # Initial guesses for fitting the decaying exponential function
        initial_guesses = [1, 0.01, np.median(diff_iqr)]

        # Fit the decaying exponential function to the IQR differences
        params, _ = curve_fit(exp_decay, sample_sizes, diff_iqr, p0=initial_guesses, maxfev=5000)

        # Generate the fitted curve
        fitted_curve = exp_decay(np.array(sample_sizes), *params)
        plt.plot(sample_sizes, fitted_curve, label=f"{treatment} exp fit", linestyle='--')

    plt.xlabel('Number of Cells')
    plt.ylabel(y_label)
    plt.legend(fontsize=20)
    plt.tight_layout()

    # Save the plot to a file if a filename is provided
    if filename is not None:
        plt.savefig(f'../plots/{filename}')

    # Display the plot
    plt.show()


def plot_cumulative_histogram_samples(data, variable_of_interest, treatment_col, treatment, x_label,
                                      initial_random_seed=42, filenames=None):
    """
    Generate cumulative histograms and plot statistical measures for a given treatment across different sample sizes.

    Parameters:
    - data: DataFrame containing the data.
    - variable_of_interest: The variable for which the statistics are calculated and plotted.
    - treatment_col: Column name indicating the treatment type in the data.
    - treatment: Specific treatment to analyze.
    - x_label: Label for the x-axis of the plot.
    - initial_random_seed: Initial random seed for reproducibility.
    - filenames: List of filenames to save the plots.
    """

    total_samples = []
    max_samples = 500
    step = 10
    filecount = 0
    random_seed = initial_random_seed

    # Filter data for the specific treatment
    subsample = data[data[treatment_col] == treatment]

    # Initialize lists to store statistical measures
    median_values = []
    mean_values = []
    std_values = []
    iqr_values = []
    sample_sizes = []

    # Iterate over sample sizes and calculate statistics
    for sample_size in range(step, max_samples + 1, step):
        # Determine the number of new samples to add
        new_samples_count = sample_size - len(total_samples)

        # Ensure we don't sample more than what's available in the dataframe
        remaining_samples = subsample[~subsample.index.isin(total_samples)].shape[0]
        new_samples_count = min(new_samples_count, remaining_samples)

        # Sample additional data and add it to the total_samples list
        if new_samples_count > 0:
            new_samples = subsample[~subsample.index.isin(total_samples)].sample(
                n=new_samples_count, replace=False, random_state=random_seed).index.tolist()
            total_samples.extend(new_samples)
            random_seed += 1

        # Extract the data for the current total samples
        sample_data = subsample.loc[total_samples, variable_of_interest]

        # Calculate statistics
        median = sample_data.median()
        mean = sample_data.mean()
        std = sample_data.std()
        q1 = sample_data.quantile(0.25)
        q3 = sample_data.quantile(0.75)
        iqr = q3 - q1

        # Store statistics
        median_values.append(median)
        mean_values.append(mean)
        std_values.append(std)
        iqr_values.append(iqr)
        sample_sizes.append(sample_size)

        # Plot histogram for specific sample sizes
        if sample_size in (20, 50, 100, 200, 300, 500):
            plt.figure(figsize=(14, 10))
            n, bins, patches = plt.hist(sample_data, bins=50, alpha=0.75, density=True)
            plt.axvline(x=median, color='r', linestyle='--', label='Median')
            plt.axvline(x=q1, color='g', linestyle='-', label='Q1')
            plt.axvline(x=q3, color='b', linestyle='-', label='Q3')

            # Calculate the density
            bin_maxes = np.maximum.reduceat(n, np.digitize([q1, q3], bins[:-1]) - 1)
            max_density = max(bin_maxes)

            plt.title(f'{len(total_samples)} {treatment} Cells')
            plt.xlabel(x_label)
            plt.ylabel('Frequency (%)')
            plt.ylim(bottom=0, top=20)
            plt.xlim(left=0, right=1)
            plt.tight_layout()

            # Save the plot if filenames are provided
            if filenames is not None:
                plt.savefig(f'../plots/{filenames[filecount]}')
            plt.show()
            filecount += 1

        # Break the loop if we have included all available samples
        if remaining_samples <= new_samples_count:
            break

    # Calculate differences relative to the full dataset statistics
    mean_values = [x - subsample[variable_of_interest].mean() for x in mean_values]
    median_values = [x - subsample[variable_of_interest].median() for x in median_values]
    std_values = [x - subsample[variable_of_interest].std() for x in std_values]
    iqr_values = [x - (subsample[variable_of_interest].quantile(0.75) - subsample[variable_of_interest].quantile(0.25))
                  for x in iqr_values]

    # Plot the differences in statistical measures
    plt.figure(figsize=(14, 10))
    plt.scatter(sample_sizes, mean_values, label='_Mean', alpha=0.5, color='blue')
    plt.scatter(sample_sizes, median_values, label='_Median', alpha=0.5, color='orange')
    plt.plot(sample_sizes, mean_values, label='Mean', color='blue')
    plt.plot(sample_sizes, median_values, label='Median', color='orange')
    plt.ylabel(f'Difference relative to statistic evaluated\n for all cells in this well')
    plt.xlabel('Number of Cells')
    plt.axhline(y=0.0, color='black', linestyle='dotted', linewidth=2)
    plt.scatter(sample_sizes, std_values, label='_Standard Deviation', alpha=0.5, color='gray')
    plt.scatter(sample_sizes, iqr_values, label='_IQR', alpha=0.5, color='purple')
    plt.plot(sample_sizes, std_values, label='Standard Deviation', color='gray')
    plt.plot(sample_sizes, iqr_values, label='IQR', color='purple')

    # Create a single legend
    plt.legend()
    plt.tight_layout()

    # Save the plot if filenames are provided
    if filenames is not None:
        plt.savefig(f'../plots/{filenames[filecount]}')
    plt.show()


def generate_superplot(plot_order, treatments, data, treatment_col, replicate_col, variable_of_interest,
                       y_label, random_seed=42, fig_width=16, fig_height=8, sample_size=-1, point_size=3,
                       filename=None):
    """
    Generate a plot visualizing the distribution of a variable of interest across different treatments and replicates.

    Parameters:
    - plot_order: Order in which treatments are displayed on the x-axis.
    - treatments: List of treatments to include in the plot.
    - data: DataFrame containing the data.
    - treatment_col: Column name indicating the treatment type in the data.
    - variable_of_interest: The variable to be plotted.
    - y_label: Label for the y-axis of the plot.
    - random_seed: Random seed for reproducibility.
    - fig_width, fig_height: Dimensions of the figure.
    - sample_size: Number of samples to take per well (if > 0).
    - point_size: Size of points in the swarm plot.
    - filename: Filename to save the plot.
    """

    # Initialize an empty DataFrame to store mean data
    mean_data = pd.DataFrame()
    replicate_index = 'Replicate_Index'

    # Iterate over each treatment and well to sample data
    for t in treatments:
        # Filter data for the current treatment
        tdata = data[data[treatment_col] == t]
        replicates = tdata[replicate_col].unique()
        for w in range(len(replicates)):
            # Sample data for each well if a sample size is specified
            if sample_size > 0:
                sdata = tdata[tdata[replicate_col] == replicates[w]].sample(n=sample_size, replace=False, random_state=random_seed)
            else:
                sdata = tdata[tdata[replicate_col] == replicates[w]]

            # Add a replicate identifier
            sdata[replicate_index] = w
            mean_data = pd.concat([mean_data, sdata])

    # Calculate the average for each replicate within each treatment
    ReplicateAverages = mean_data.groupby([treatment_col, replicate_index], as_index=False).agg(
        {variable_of_interest: "mean"})

    # Generate the plot
    plt.figure(figsize=(fig_width, fig_height))
    ax = plt.subplot(1, 1, 1)

    # Create a box plot for the replicate averages
    sns.boxplot(x=treatment_col, y=variable_of_interest, data=ReplicateAverages, order=plot_order, color='white',
                showfliers=False, linecolor='black', linewidth=2, zorder=1, boxprops=dict(facecolor='none'))

    # Create a swarm plot for the individual data points
    sns.swarmplot(x=treatment_col, y=variable_of_interest, hue=replicate_index, data=mean_data, size=point_size,
                  order=plot_order, zorder=0, palette={0: 'cornflowerblue', 1: 'gray', 2: 'orange'})

    # Create a swarm plot for the replicate averages with larger points
    sns.swarmplot(x=treatment_col, y=variable_of_interest, hue=replicate_index, size=20, edgecolor="k", linewidth=2,
                  data=ReplicateAverages, order=plot_order, zorder=2,
                  palette={0: 'cornflowerblue', 1: 'gray', 2: 'orange'})

    # Set plot limits and labels
    plt.ylim(bottom=0.0, top=1.0)
    plt.xlabel('')
    plt.ylabel(y_label)
    plt.title(f'{sample_size} cells per population')

    # Remove the legend for clarity
    ax.legend_.remove()

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot to a file if a filename is provided
    if filename is not None:
        plt.savefig(f'../plots/{filename}')

    # Display the plot
    plt.show()
