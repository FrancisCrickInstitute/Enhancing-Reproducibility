import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def normalize_well_format(well):
    # Use regular expression to extract the letter and number parts
    match = re.match(r"([A-Za-z])([0-9]+)", well, re.I)
    if match:
        items = match.groups()
        letter = items[0]
        number = str(int(items[1])).zfill(2)  # Convert to integer to remove leading zeros, then back to string
        return letter + number
    return well  # Return the original string if no match is found


annotations = pd.read_csv(
    'D:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/idr0139-screenA-annotation.csv')
annotations = annotations[annotations['Plate'] == 1093711385]
annotations['Control Type'] = annotations['Control Type'].fillna('Treated').replace('', 'Treated')
annotations['Well'] = annotations['Well'].apply(normalize_well_format)

treatments = annotations.set_index('Well')['Control Type'].to_dict()
image_data = pd.read_csv(
    'Z:/working/barryd/IDR/Outputs/Image.csv')
nuc_data = pd.read_csv(
    'Z:/working/barryd/IDR/Outputs/Nuclei.csv')
cyto_data = pd.read_csv(
    'Z:/working/barryd/IDR/Outputs/Cytoplasm.csv')

nuc_data = nuc_data.rename(columns=lambda x: 'Nuclear_' + x if 'Intensity' in x else x)
cyto_data = cyto_data.rename(columns=lambda x: 'Cyto_' + x if 'Intensity' in x else x)

nuc_cyto_data = pd.merge(nuc_data, cyto_data, on=['ImageNumber', 'ObjectNumber'], how='left')
combined_data = pd.merge(nuc_cyto_data, image_data, on='ImageNumber', how='left')

combined_data['Fascin_Ratio'] = combined_data['Nuclear_Intensity_MeanIntensity_Fascin'] / combined_data[
    'Cyto_Intensity_MeanIntensity_Fascin']
combined_data['NuclearActin_Ratio'] = combined_data['Nuclear_Intensity_MeanIntensity_NuclearActin'] / combined_data[
    'Cyto_Intensity_MeanIntensity_NuclearActin']

combined_data['Well'] = combined_data['FileName_DNA'].str.extract(r'_(.*?)_')
combined_data['Treatment'] = combined_data['Well'].map(treatments)

# Sort DataFrame by Treatment
combined_data = combined_data.sort_values(by=['Treatment', 'Well'])

# Determine the order of x-axis categories
well_order = combined_data['Well'].unique()

# plt.figure(figsize=(20, 12))
# Boxplot to show distribution
# sns.boxplot(x='Well', y='Intensity_MeanIntensity_Fascin', data=combined_data, color='black', log_scale=10)

# Swarmplot to show individual data points
# sns.swarmplot(x='Well', y='Fascin_Ratio', data=combined_data, order=well_order,
#             palette="pastel",
#            hue='Treatment',
#           size=0.1,
#          log_scale=10)

# sns.boxplot(x='Well', y='Fascin_Ratio', data=combined_data, order=well_order,
#       color='white',
#      log_scale=10, showfliers=False, showmeans=True,
#     meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"})

# handles, labels = plt.gca().get_legend_handles_labels()
# large_markers = [mlines.Line2D([], [], color=handle.get_edgecolor(), marker='o', markersize=10, linestyle='None') for
#         handle in handles]
# legend = plt.legend(handles=large_markers, labels=labels, loc='center left', bbox_to_anchor=(1, 0.5), title='Treatment')

# plt.title('Stripplot of Intensity_MeanIntensity_Fascin grouped by Well (Log Scale)')
# plt.savefig("output_filename.pdf", format='pdf', bbox_inches='tight')

columns_to_keep = ['Well'] + combined_data.select_dtypes(include=np.number).columns.tolist()
df_numeric = combined_data[columns_to_keep]

# 2. Aggregate the intensity measurements at the well level (using mean as an example)
df_well_level = df_numeric.groupby('Well').median().reset_index()
df_well_level['Treatment'] = df_well_level['Well'].map(treatments)

plt.figure(figsize=(20, 12))
# Boxplot to show distribution
# sns.boxplot(x='Well', y='Intensity_MeanIntensity_Fascin', data=combined_data, color='black', log_scale=10)

# Swarmplot to show individual data points
sns.swarmplot(x='Treatment', y='NuclearActin_Ratio', data=df_well_level,
              size=10, palette="pastel", hue='Treatment')

sns.boxplot(x='Treatment', y='NuclearActin_Ratio', data=df_well_level, color='white',
            showfliers=False)

# plt.title('Stripplot of Intensity_MeanIntensity_Fascin grouped by Well (Log Scale)')
plt.savefig("agg_output_filename.pdf", format='pdf', bbox_inches='tight')

plt.figure(figsize=(20, 12))

sns.scatterplot(x='Fascin_Ratio', y='NuclearActin_Ratio', data=df_well_level, size=20,
                palette="pastel", hue='Treatment')

plt.savefig("scatter_output_filename.pdf", format='pdf', bbox_inches='tight')
