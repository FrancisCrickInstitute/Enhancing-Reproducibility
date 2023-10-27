import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

treatments = {
    # Compounds
    'B10': 'Treated',
    'B12': 'Treated',
    'B20': 'Treated',
    'B08': 'Treated',
    'D16': 'Treated',
    'D18': 'Treated',
    'E10': 'Treated',
    'E16': 'Treated',
    'E20': 'Treated',
    'E04': 'Treated',
    'E06': 'Treated',
    'H10': 'Treated',
    'H16': 'Treated',
    'H22': 'Treated',
    'I02': 'Treated',
    'I04': 'Treated',
    # Controls
    'B02': 'DMSO',
    'B22': 'DMSO',
    'C12': 'DMSO',
    'E17': 'DMSO',
    'E07': 'DMSO',
    'H12': 'DMSO',
    'H20': 'DMSO',
    'H04': 'DMSO',
    'L17': 'DMSO',
    'L07': 'DMSO',
    'N12': 'DMSO',
    'O02': 'DMSO',
    'O22': 'DMSO',
    # Stimulator
    'B23': 'Stimulator Control',
    'B03': 'Stimulator Control',
    'C13': 'Stimulator Control',
    'E18': 'Stimulator Control',
    'E08': 'Stimulator Control',
    'H13': 'Stimulator Control',
    'H21': 'Stimulator Control',
    'H05': 'Stimulator Control',
    'L18': 'Stimulator Control',
    'L08': 'Stimulator Control',
    'N13': 'Stimulator Control',
    'O23': 'Stimulator Control',
    'O03': 'Stimulator Control'
}

image_data = pd.read_csv(
    'E:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/Outputs/Image.csv')
nuc_data = pd.read_csv(
    'E:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/Outputs/Nuclei.csv')

combined_data = pd.merge(nuc_data, image_data, on='ImageNumber', how='left')

combined_data['Well'] = combined_data['FileName_DNA'].str.extract(r'_(.*?)_')
combined_data['Treatment'] = combined_data['Well'].map(treatments)

# Sort DataFrame by Treatment
combined_data = combined_data.sort_values(by=['Treatment', 'Well'])

# Determine the order of x-axis categories
well_order = combined_data['Well'].unique()

plt.figure(figsize=(20, 12))
# Boxplot to show distribution
# sns.boxplot(x='Well', y='Intensity_MeanIntensity_Fascin', data=combined_data, color='black', log_scale=10)

# Swarmplot to show individual data points
sns.swarmplot(x='Well', y='Intensity_IntegratedIntensity_NuclearActin', data=combined_data, order=well_order, palette="pastel",
              hue='Treatment',
              size=0.5,
              log_scale=10)

sns.boxplot(x='Well', y='Intensity_IntegratedIntensity_NuclearActin', data=combined_data, order=well_order, color='white',
            log_scale=10, showfliers=False, showmeans=True,
            meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"})

# plt.title('Stripplot of Intensity_MeanIntensity_Fascin grouped by Well (Log Scale)')
plt.savefig("output_filename.pdf", format='pdf', bbox_inches='tight')

columns_to_keep = ['Well'] + combined_data.select_dtypes(include=np.number).columns.tolist()
df_numeric = combined_data[columns_to_keep]

# 2. Aggregate the intensity measurements at the well level (using mean as an example)
df_well_level = df_numeric.groupby('Well').median().reset_index()
df_well_level['Treatment'] = df_well_level['Well'].map(treatments)

plt.figure(figsize=(20, 12))
# Boxplot to show distribution
# sns.boxplot(x='Well', y='Intensity_MeanIntensity_Fascin', data=combined_data, color='black', log_scale=10)

# Swarmplot to show individual data points
sns.swarmplot(x='Treatment', y='Intensity_IntegratedIntensity_NuclearActin', data=df_well_level,
              size=3)

sns.boxplot(x='Treatment', y='Intensity_IntegratedIntensity_NuclearActin', data=df_well_level, color='white',
            showfliers=False)

# plt.title('Stripplot of Intensity_MeanIntensity_Fascin grouped by Well (Log Scale)')
plt.savefig("agg_output_filename.pdf", format='pdf', bbox_inches='tight')
