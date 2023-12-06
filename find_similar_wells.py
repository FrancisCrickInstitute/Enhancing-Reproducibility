import numpy as np
import pandas as pd
from scipy.stats import t


def calculate_t_test(df, well1, well2):
    # Extract data for the two wells
    data1 = df[df['Well'] == well1]
    data2 = df[df['Well'] == well2]

    # Means, standard deviations, and sample sizes
    mean1, mean2 = data1['mean'].values[0], data2['mean'].values[0]
    std1, std2 = data1['std'].values[0], data2['std'].values[0]
    n1, n2 = data1['count'].values[0], data2['count'].values[0]

    # Calculate the t-statistic
    se = np.sqrt((std1 ** 2 / n1) + (std2 ** 2 / n2))
    t_stat = (mean1 - mean2) / se

    # Degrees of freedom
    df = ((std1 ** 2 / n1 + std2 ** 2 / n2) ** 2) / (
                ((std1 ** 2 / n1) ** 2 / (n1 - 1)) + ((std2 ** 2 / n2) ** 2 / (n2 - 1)))

    # Calculate the two-tailed p-value
    p_value = t.sf(np.abs(t_stat), df) * 2

    return t_stat, p_value


# Load the dataframe (replace with your file path)
file_path = './plots/descriptive_stats.csv'
df = pd.read_csv(file_path)

# Step 1: Identify Negative, Neutral, and Treated Controls
negative_controls = df[df['Treatment'] == 'Negative Control']
neutral_controls = df[df['Treatment'] == 'Neutral Control']
treated_wells = df[df['Treatment'] == 'Treated']

# Step 2: Find tuples of wells that meet the criteria, quantify the differences, and calculate total difference
matching_tuples = []
for _, neg_row in negative_controls.iterrows():
    for _, neu_row in neutral_controls.iterrows():
        neg_neu_overlap = neu_row['95% CI lower'] - neg_row['95% CI upper']
        if (neg_neu_overlap <= 0) & (neu_row['mean'] > neg_row['mean']):  # Indicates overlap or just touching
            for _, tre_row in treated_wells.iterrows():
                neu_tre_overlap = tre_row['95% CI lower'] - neu_row['95% CI upper']
                if (neu_tre_overlap <= 0) & (tre_row['mean'] > neu_row['mean']):  # Indicates overlap or just touching
                    total_diff = abs(neg_neu_overlap) + abs(neu_tre_overlap)
                    matching_tuples.append({
                        'Negative Control': neg_row['Well'],
                        'Neutral Control': neu_row['Well'],
                        'Treated': tre_row['Well'],
                        'Neg-Neu Overlap': neg_neu_overlap,
                        'Neu-Tre Overlap': neu_tre_overlap,
                        'Total Difference': total_diff
                    })

# Adding t-test results to each tuple
for tuple in matching_tuples:
    neg_well, neu_well, tre_well = tuple['Negative Control'], tuple['Neutral Control'], tuple['Treated']

    # T-test between Negative Control and Neutral Control
    t_stat_neg_neu, p_value_neg_neu = calculate_t_test(df, neg_well, neu_well)

    # T-test between Neutral Control and Treated
    t_stat_neu_tre, p_value_neu_tre = calculate_t_test(df, neu_well, tre_well)

    # Add t-test results to the tuple
    tuple['T-stat Neg-Neu'] = t_stat_neg_neu
    tuple['P-value Neg-Neu'] = p_value_neg_neu
    tuple['T-stat Neu-Tre'] = t_stat_neu_tre
    tuple['P-value Neu-Tre'] = p_value_neu_tre

# Convert the results to a DataFrame and sort by total difference
results_df = pd.DataFrame(matching_tuples)
results_df = results_df.sort_values(by='Total Difference')

# Print the sorted results
print(results_df)

results_df.to_csv('./plots/significance.csv', index=False)

