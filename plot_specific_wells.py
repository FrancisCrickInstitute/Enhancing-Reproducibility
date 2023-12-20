import math
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import scipy.stats as stats
import seaborn as sns

plt.rcParams['font.size'] = 16


def normalize_well_format(well):
    match = re.match(r"([A-Za-z])([0-9]+)", well, re.I)
    if match:
        items = match.groups()
        return f"{items[0]}{int(items[1]):02d}"
    return well


def load_and_prepare_data(file_path, plate_number):
    df = pd.read_csv(file_path)
    df = df[df['Plate'] == plate_number]
    df['Control Type'] = df['Control Type'].fillna('Treated').replace('', 'Treated')
    df['Well'] = df['Well'].apply(normalize_well_format)
    return df


sample_size = 500
n_samples = 4
treatment_col = 'Treatment'
variable_of_interest = 'Fascin_Ratio'
dunn_pairs = [('Untreated', 'DMSO'), ('DMSO', 'SN0212398523'), ('SN0212398523', 'Leptomycin b'),
              ('Untreated', 'SN0212398523')]
dunn_p_values = {pair: [[] for _ in range(n_samples)] for pair in dunn_pairs}

annotations = load_and_prepare_data(
    'E:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/idr0139-screenA-annotation.csv',
    1093711385)
treatments = annotations.set_index('Well')['Control Type'].to_dict()

treatments_to_compounds = {'Treated': 'SN0212398523', 'Negative Control': 'Untreated', 'Neutral Control': 'DMSO',
                           'Stimulator Control': 'Leptomycin b'}

image_data = pd.read_csv('Z:/working/barryd/IDR/Outputs/Image.csv')
nuc_data = pd.read_csv('Z:/working/barryd/IDR/Outputs/Nuclei.csv')
cyto_data = pd.read_csv('Z:/working/barryd/IDR/Outputs/Cytoplasm.csv')

# Rename columns
nuc_data = nuc_data.rename(columns=lambda x: 'Nuclear_' + x if 'Intensity' in x else x)
cyto_data = cyto_data.rename(columns=lambda x: 'Cyto_' + x if 'Intensity' in x else x)

# Merge data
combined_data = pd.merge(pd.merge(nuc_data, cyto_data, on=['ImageNumber', 'ObjectNumber'], how='left'), image_data,
                         on='ImageNumber', how='left')

# Calculate ratios
for compartment in ['Fascin', 'NuclearActin']:
    combined_data[f'{compartment}_Ratio'] = combined_data[f'Nuclear_Intensity_MeanIntensity_{compartment}'] / (
            combined_data[f'Cyto_Intensity_MeanIntensity_{compartment}'] + combined_data[
        f'Nuclear_Intensity_MeanIntensity_{compartment}'])

# Extract well information and map treatments
combined_data['Well'] = combined_data['FileName_DNA'].str.extract(r'_(.*?)_')[0]
combined_data['Treatment'] = combined_data['Well'].map(treatments).map(treatments_to_compounds)

selected_wells = ['J05', 'O02', 'E22', 'L08']
selected_wells_data = combined_data[combined_data['Well'].isin(selected_wells)]
selected_wells_data = selected_wells_data.sort_values(by=['Treatment', 'Well'])

color_dict = {'SN0212398523': 'orange', 'Untreated': 'blue', 'DMSO': 'gray', 'Leptomycin b': 'purple'}

# Plotting
well_order = selected_wells_data['Well'].unique()
treatment_order = selected_wells_data['Treatment'].unique()
# fig, ax = plt.subplots(figsize=(14, 10))
#
# sns.swarmplot(x='Treatment', y='Fascin_Ratio', data=selected_wells_data,
#               order=['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'], palette=color_dict,
#               hue='Treatment', size=2.75, alpha=0.9, ax=ax)
# sns.boxplot(x='Treatment', y='Fascin_Ratio', data=selected_wells_data,
#             order=['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'], color='white',
#             showfliers=False, ax=ax)
#
# plt.ylabel('Log[Mean_Nuclear_Fascin_Intensity /\n (Mean_Nuclear_Fascin_Intensity + Mean_Cytoplasmic_Fascin_Intensity)]')
# plt.xlabel('')
# plt.show()

# plt.savefig("./plots/selected_wells_all_cells.pdf", format='pdf', bbox_inches='tight')

fig = plt.figure(num=1, figsize=(28, 20))

for sample_index, _ in enumerate(range(n_samples)):

    untreated_data = selected_wells_data[selected_wells_data['Treatment'] == 'Untreated'].sample(n=sample_size,
                                                                                                 replace=False)
    dmso_data = selected_wells_data[selected_wells_data['Treatment'] == 'DMSO'].sample(n=sample_size, replace=False)
    treated_data = selected_wells_data[selected_wells_data['Treatment'] == 'SN0212398523'].sample(n=sample_size,
                                                                                                  replace=False)
    stim_data = selected_wells_data[selected_wells_data['Treatment'] == 'Leptomycin b'].sample(n=sample_size,
                                                                                               replace=False)

    sampled_data = pd.concat([untreated_data, dmso_data, treated_data, stim_data])

    ax = plt.subplot(2, 2, sample_index + 1)

    sns.swarmplot(x='Treatment', y='Fascin_Ratio', data=sampled_data,
                  order=['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'], palette=color_dict,
                  hue='Treatment', size=2.75, alpha=0.9, ax=ax)
    sns.boxplot(x='Treatment', y='Fascin_Ratio', data=sampled_data,
                order=['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'], color='white',
                showfliers=False, ax=ax)

    plt.ylabel(
        'Log[Mean_Nuclear_Fascin_Intensity /\n (Mean_Nuclear_Fascin_Intensity + Mean_Cytoplasmic_Fascin_Intensity)]')
    plt.xlabel('')

    _, p_value = stats.kruskal(
        *(sampled_data[sampled_data[treatment_col] == t][variable_of_interest] for t in
          sampled_data[treatment_col].unique()))

    if p_value < 0.05:
        dunn_result = sp.posthoc_dunn(sampled_data, val_col=variable_of_interest, group_col=treatment_col)
        for pair in dunn_pairs:
            dunn_p_values[pair][sample_index].append(dunn_result.loc[pair[0], pair[1]])
    else:
        for pair in dunn_pairs:
            dunn_p_values[pair][sample_index].append(np.nan)

    y, h, col = sampled_data[variable_of_interest].max() + 0.02, 0.02, 'black'

    ymax = []
    for t in range(len(treatments_to_compounds) - 1):
        ymax.append(0)

    for pair in dunn_pairs:
        x1, x2 = pair
        x1 = [label.get_text() for label in ax.get_xticklabels()].index(x1)
        x2 = [label.get_text() for label in ax.get_xticklabels()].index(x2)

        y = sampled_data[sampled_data['Treatment'].isin(pair)].loc[:, variable_of_interest].max() + 0.02

        for x in range(min(x1, x2), max(x1, x2)):
            if y <= ymax[x] + 0.075:
                y = ymax[x] + 0.075
            ymax[x] = y

        if x1 < x2:
            x1 += 0.02
            x2 -= 0.02
        else:
            x1 -= 0.02
            x2 += 0.02

        # Draw line
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)

        raw_p_value = dunn_p_values[pair][sample_index][0]
        str_p_value = f'p = {raw_p_value:.4f}'
        if raw_p_value < 0.0001:
            str_p_value = 'p < 0.0001'
        elif raw_p_value < 0.001:
            str_p_value = 'p < 0.001'
        elif raw_p_value < 0.01:
            str_p_value = 'p < 0.01'
        elif raw_p_value < 0.05:
            str_p_value = 'p < 0.05'

        # Annotate line with p-value
        plt.text((x1 + x2) * .5, y + h, str_p_value, ha='center', va='bottom',
                 color=col)

plt.show()

# untreated_data.to_csv('./plots/untreated_raw_data.csv', index=False)
# dmso_data.to_csv('./plots/dmso_raw_data.csv', index=False)
# treated_data.to_csv('./plots/treated_raw_data.csv', index=False)
# stim_data.to_csv('./plots/stim_raw_data.csv', index=False)
#
# # Calculate descriptive statistics for each well
# descriptive_stats = selected_wells_data.groupby('Well')['Fascin_Ratio'].describe(percentiles=[.25, .5, .75])
#
# # Calculate 95% confidence interval
# descriptive_stats['95% CI lower'] = descriptive_stats['mean'] - 1.96 * (
#         descriptive_stats['std'] / np.sqrt(descriptive_stats['count']))
# descriptive_stats['95% CI upper'] = descriptive_stats['mean'] + 1.96 * (
#         descriptive_stats['std'] / np.sqrt(descriptive_stats['count']))
#
# # Calculate inter-quartile range
# descriptive_stats['IQR'] = descriptive_stats['75%'] - descriptive_stats['25%']
#
# # Reset index to merge with treatment information
# descriptive_stats.reset_index(inplace=True)
#
# # Merge with treatment information
# descriptive_stats = pd.merge(descriptive_stats, selected_wells_data[['Well', 'Treatment']].drop_duplicates(), on='Well',
#                              how='left')
#
# # Save the DataFrame with treatment information to a CSV file
# descriptive_stats.to_csv('./plots/selected_wells_descriptive_stats.csv', index=False)
# selected_wells_data.to_csv('./plots/selected_wells_raw_data.csv', index=False)
