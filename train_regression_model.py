import re
import numpy as np
import pandas as pd
import shap
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


def normalize_well_format(well):
    # Use regular expression to extract the letter and number parts
    match = re.match(r"([A-Za-z])([0-9]+)", well, re.I)
    if match:
        items = match.groups()
        letter = items[0]
        number = str(int(items[1])).zfill(2)  # Convert to integer to remove leading zeros, then back to string
        return letter + number
    return well  # Return the original string if no match is found


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
cyto_data = pd.read_csv(
    'E:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/Outputs/Cytoplasm.csv')
cell_data = pd.read_csv(
    'E:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/Outputs/Cells.csv')
scores = pd.read_csv(
    'E:/OneDrive - The Francis Crick Institute/Publications/2023_Dont_Trust_P_Values/idr0139-screenA-annotation.csv')

scores = scores[scores['Plate'] == 1093711385]

nuc_data = nuc_data.rename(columns=lambda x: 'Nuclear_' + x if 'Intensity' in x else x)
cyto_data = cyto_data.rename(columns=lambda x: 'Cyto_' + x if 'Intensity' in x else x)
cell_data = cell_data.rename(columns=lambda x: 'Cell_' + x if 'Intensity' in x else x)

#nuc_cell_data = pd.merge(nuc_data, cell_data, on=['ImageNumber', 'ObjectNumber'], how='left')

#nuc_cell_cyto_data = pd.merge(nuc_cell_data, cyto_data, on=['ImageNumber', 'ObjectNumber'], how='left')

#combined_data = pd.merge(nuc_cell_cyto_data, image_data, on='ImageNumber', how='left')
combined_data = pd.merge(nuc_data, image_data, on='ImageNumber', how='left')

combined_data['Well'] = combined_data['FileName_DNA'].str.extract(r'_(.*?)_')
combined_data['Treatment'] = combined_data['Well'].map(treatments)

# Sort DataFrame by Treatment
combined_data = combined_data.sort_values(by=['Treatment', 'Well'])

combined_data['Well'] = combined_data['Well'].apply(normalize_well_format)
scores['Well'] = scores['Well'].apply(normalize_well_format)

combined_data = pd.merge(combined_data, scores, on='Well', how='left')

# 1. Remove non-numeric columns (except for 'Well' and 'Score')
columns_to_keep = ['Well', 'Score'] + combined_data.select_dtypes(include=np.number).columns.tolist()
df_numeric = combined_data[columns_to_keep]

# 1. Separate the Score column and intensity measurements
score_data = combined_data[['Well', 'Score']].drop_duplicates()
intensity_data = df_numeric.drop(columns='Score')

# 2. Aggregate the intensity measurements at the well level (using mean as an example)
df_well_level = intensity_data.groupby('Well').median().reset_index()

# 3. Merge back the Score data
df_well_level = pd.merge(df_well_level, score_data, on='Well')

# Prepare the data
X = df_well_level.filter(like='Intensity')  # Select columns that contain 'Intensity'
y = df_well_level['Score']  # Target variable (scores)

# Split the data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

# Choose the model
model = LinearRegression()

# Train the model
model.fit(X_train, y_train)

# Initialize the SHAP explainer
explainer = shap.Explainer(model.predict, X)

# Calculate SHAP values for a specific well
# Here, I'm using the first well in the test set as an example
shap_values = explainer.shap_values(X)

# Show the feature importance for that specific well
shap.summary_plot(shap_values, X, feature_names=X.columns)
