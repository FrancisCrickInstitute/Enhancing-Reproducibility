import pandas as pd

df_1A = pd.read_csv('./LM2_GEFGAP_ONTARGETPlus_1A_instances.csv')
df_2A = pd.read_csv('./LM2_GEFGAP_ONTARGETPlus_2A_instances.csv')
df_2B = pd.read_csv('./LM2_GEFGAP_ONTARGETPlus_2B_instances.csv')

combined_df = pd.concat([df_1A, df_2A, df_2B])
combined_summary_df = combined_df.groupby('Treatment').agg({
    'YAPTAZ_Ratio': 'mean',  # Calculate the average 'YAPTAZ_Ratio'
    'Well': 'count'  # Count the number of instances
})

combined_summary_df.to_csv('./combined_summary.csv')
