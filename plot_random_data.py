from utility_functions import generate_swarmplot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.size'] = 28
plate_number = 1093711385
treatment_col = 'Treatment'
variable_of_interest = 'Fascin_Ratio'
treatments_to_compounds = {'Dataset 1': 'Negative Control', 'Dataset 2': 'Condition 1', 'Dataset 3': 'Condition 2'}
dunn_pairs = [('Negative Control', 'Condition 1'), ('Condition 1', 'Condition 2')]
color_dict = {'Negative Control': 'orange', 'Condition 1': 'blue', 'Condition 2': 'gray'}

# Generating three different datasets, each normally distributed around 1.0
data1 = np.random.normal(loc=0.995, scale=0.18, size=1500)
data2 = np.random.normal(loc=1.025, scale=0.22, size=1500)
data3 = np.random.normal(loc=1.01, scale=0.2, size=1500)

df = pd.DataFrame({
    'Value': np.concatenate([data1, data2, data3]),
    'Dataset': np.repeat(['Negative Control', 'Condition 1', 'Condition 2'], len(data1))
})

generate_swarmplot(14, 10, 1, 1, ['Negative Control', 'Condition 1', 'Condition 2'], 1, -1,
                   df, color_dict, 'Dataset',
                   'Value', dunn_pairs, treatments_to_compounds, 'Mean Nuclear Intensity /\n Mean Cytoplasmic Intensity')
