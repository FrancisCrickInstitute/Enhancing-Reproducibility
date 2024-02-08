from utility_functions import *

plt.rcParams['font.size'] = 16
plate_number = 1093711385
treatment_col = 'Treatment'
variable_of_interest = 'Fascin_Ratio'
treatments_to_compounds = {'Treated': 'SN0212398523', 'Negative Control': 'Untreated', 'Neutral Control': 'DMSO',
                           'Stimulator Control': 'Leptomycin b'}
dunn_pairs = [('Untreated', 'DMSO'), ('DMSO', 'SN0212398523'), ('SN0212398523', 'Leptomycin b'),
              ('Untreated', 'SN0212398523')]
color_dict = {'SN0212398523': 'orange', 'Untreated': 'blue', 'DMSO': 'gray', 'Leptomycin b': 'purple'}

directories = ('./inputs/idr', './outputs/plots', './outputs/data')

for d in directories:
    if not os.path.exists(d):
        os.makedirs(d)

idr_annotations_file_path = './inputs/idr/idr0139-screenA-annotation.csv'
idr_annotations_url = 'https://raw.githubusercontent.com/IDR/idr0139-lawson-fascin/main/screenA/idr0139-screenA-annotation.csv'

download_csv(idr_annotations_file_path, idr_annotations_url)
annotations = load_and_prepare_data(idr_annotations_file_path, plate_number)
image_data = pd.read_csv('./inputs/cell_profiler_outputs/Image.csv')
nuc_data = pd.read_csv('./inputs/cell_profiler_outputs/Nuclei.csv')
cyto_data = pd.read_csv('./inputs/cell_profiler_outputs/Cytoplasm.csv')
treatments = annotations.set_index('Well')['Control Type'].to_dict()

data_subset = prepare_data(nuc_data, cyto_data, image_data, treatments, treatments_to_compounds,
                           ['J05', 'O02', 'E22', 'L08'])

# data_subset.to_csv('./all_data.csv')

# # FIGURE 1
# generate_swarmplot(14, 10, 1, 1, ['Untreated', 'SN0212398523'], 1, -1,
#                    data_subset[data_subset['Well'].isin(['J05', 'E22'])], color_dict, treatment_col,
#                    variable_of_interest, [dunn_pairs[3]], treatments_to_compounds)
#
# # FIGURE 2
# generate_swarmplot(14, 10, 1, 1, ['Untreated', 'DMSO', 'SN0212398523'], 1, -1,
#                    data_subset[data_subset['Well'].isin(['J05', 'O02', 'E22'])], color_dict, treatment_col,
#                    variable_of_interest, [dunn_pairs[0], dunn_pairs[1], dunn_pairs[3]], treatments_to_compounds)
#
# # FIGURE 3
# generate_swarmplot(14, 10, 1, 1, ['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'],
#                    1, -1, data_subset, color_dict, treatment_col, variable_of_interest, dunn_pairs,
#                    treatments_to_compounds)
#
# # TABLE 1
# generate_table(data_subset)
#
# # FIGURE 4
# plot_mean_v_sample_size([*range(10, 500, 10)], 100, data_subset, treatment_col, variable_of_interest)
#
# # FIGURE 5
# plot_p_v_sample_size([*range(10, 500, 10)], 100, data_subset, treatment_col, variable_of_interest, dunn_pairs)
#
# # FIGURE 6
# generate_swarmplot(28, 20, 2, 2, ['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'],
#                    4, 50, data_subset, color_dict, treatment_col, variable_of_interest, dunn_pairs,
#                    treatments_to_compounds)
#
# # FIGURE 7
# generate_swarmplot(28, 20, 2, 2, ['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'],
#                    4, 100, data_subset, color_dict, treatment_col, variable_of_interest, dunn_pairs,
#                    treatments_to_compounds)
#
# # FIGURE 8
# generate_swarmplot(28, 20, 2, 2, ['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'],
#                    4, 250, data_subset, color_dict, treatment_col, variable_of_interest, dunn_pairs,
#                    treatments_to_compounds)
#
# # FIGURE 9
# for i in range(4):
#     generate_swarmplot(14, 10, 1, 1, ['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'],
#                        1, 500, data_subset, color_dict, treatment_col, variable_of_interest, dunn_pairs,
#                        treatments_to_compounds, 'Relative Nuclear Fascin Localisation')

plot_iqr_v_sample_size([*range(10, 500, 10)], 100, data_subset, treatment_col, variable_of_interest, 'Error in Inter-Quartile Range')

plot_cumulative_histogram_samples(data_subset, variable_of_interest, treatment_col, 'Untreated')

print('All Done!')
