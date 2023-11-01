import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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


annotations = load_and_prepare_data(
    'D:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/idr0139-screenA-annotation.csv',
    1093711385)
treatments = annotations.set_index('Well')['Control Type'].to_dict()

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
combined_data['Treatment'] = combined_data['Well'].map(treatments)

# Filter data
untreated_data = combined_data.query("Treatment == 'Negative Control'")
other_data = combined_data.query("Treatment != 'Negative Control'")

unique_wells = untreated_data['Well'].unique()
selected_wells = np.random.choice(unique_wells, min(10, len(unique_wells)), replace=False)
filtered_untreated_data = untreated_data[untreated_data['Well'].isin(selected_wells)]

filtered_data = pd.concat([filtered_untreated_data, other_data])
filtered_data = filtered_data.sort_values(by=['Treatment', 'Well'])

# Plotting
well_order = filtered_data['Well'].unique()
fig, ax = plt.subplots(figsize=(28, 12))

sns.swarmplot(x='Well', y='Fascin_Ratio', data=filtered_data, order=well_order, palette="pastel", hue='Treatment',
              size=0.4, log_scale=10, ax=ax)
sns.boxplot(x='Well', y='Fascin_Ratio', data=filtered_data, order=well_order, color='white', log_scale=10,
            showfliers=False, ax=ax)

plt.ylabel('Log[Mean_Nuclear_Fascin_Intensity / (Mean_Nuclear_Fascin_Intensity + Mean_Cytoplasmic_Fascin_Intensity)]')

handles, labels = ax.get_legend_handles_labels()
palette = sns.color_palette("pastel", n_colors=filtered_data['Treatment'].nunique())
legend_patches = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=palette[i], markersize=10) for i in
                  range(len(palette))]
plt.legend(handles=legend_patches, labels=np.unique(filtered_data['Treatment']).tolist())

plt.savefig("./plots/all_cells.pdf", format='pdf', bbox_inches='tight')

# ----------------------------------------------------------------------------------------------------------------------

columns_to_keep = ['Well'] + combined_data.select_dtypes(include=np.number).columns.tolist()
df_numeric = combined_data[columns_to_keep]
df_well_level = df_numeric.groupby('Well').median().reset_index()
df_well_level['Treatment'] = df_well_level['Well'].map(treatments)

plt.figure(figsize=(20, 12))
sns.swarmplot(x='Treatment', y='NuclearActin_Ratio', data=df_well_level,
              size=10, palette="pastel", hue='Treatment')
sns.boxplot(x='Treatment', y='NuclearActin_Ratio', data=df_well_level, color='white',
            showfliers=False)
plt.savefig("./plots/well_medians.pdf", format='pdf', bbox_inches='tight')

# ----------------------------------------------------------------------------------------------------------------------

plt.figure(figsize=(20, 12))
sns.scatterplot(x='Fascin_Ratio', y='NuclearActin_Ratio', data=df_well_level,
                palette="pastel", hue='Treatment')
plt.savefig("./plots/scatter_output_filename.pdf", format='pdf', bbox_inches='tight')
