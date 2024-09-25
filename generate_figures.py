from utility_functions import *

plt.rcParams['font.size'] = 26
plt.rcParams['axes.linewidth'] = 2
plate_number = 1093711385
treatment_col = 'Treatment'
variable_of_interest = 'Fascin_Ratio'
dunn_pairs = [('Untreated', 'DMSO'), ('DMSO', 'SN0212398523'), ('SN0212398523', 'Leptomycin b'),
              ('Untreated', 'SN0212398523'), ('Untreated', 'Leptomycin b'), ('DMSO', 'Leptomycin b')]
color_dict = {'SN1066932540': 'red', 'SN1054616339': 'yellow', 'SN0212398523': 'orange', 'Untreated': 'blue',
              'DMSO': 'gray',
              'Leptomycin b': 'purple'}

directories = ('./inputs/idr', './outputs/plots', './outputs/data')

for d in directories:
    if not os.path.exists(d):
        os.makedirs(d)

idr_annotations_file_path = './inputs/idr/idr0139-screenA-annotation.csv'
idr_annotations_url = 'https://raw.githubusercontent.com/IDR/idr0139-lawson-fascin/main/screenA/idr0139-screenA-annotation.csv'

download_csv(idr_annotations_file_path, idr_annotations_url)
annotations = load_and_prepare_data(idr_annotations_file_path, plate_number)
compounds = annotations[annotations['Control Type'] == 'Treated'].set_index('Well')['Proprietary Compound'].to_dict()
treatments_to_compounds = {'Treated': 'Treated', 'Negative Control': 'Untreated', 'Neutral Control': 'DMSO',
                           'Stimulator Control': 'Leptomycin b'}
image_data = pd.read_csv('./inputs/cell_profiler_outputs/Image.csv')
nuc_data = pd.read_csv('./inputs/cell_profiler_outputs/Nuclei.csv')
cyto_data = pd.read_csv('./inputs/cell_profiler_outputs/Cytoplasm.csv')
treatments = annotations.set_index('Well')['Control Type'].to_dict()

data_subset = prepare_data(nuc_data, cyto_data, image_data, treatments, treatments_to_compounds, compounds,
                           ['J05', 'O02', 'E22', 'L08'])

# # FIGURE 2 A - F
# for s in [50, 200]:
#     for i in range(3):
generate_swarmplot(14, 10, 1, 1, ['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'],
                   1, -1, data_subset, color_dict, treatment_col, variable_of_interest,
                   '$ \\log \\left[ \\frac {I_{F_N}}{(I_{F_N} + I_{C_N})} \\right]$', dunn_pairs)
#
# # FIGURE 2 G - I
# plot_effect_size_v_sample_size([*range(10, 500, 10)], 100, data_subset, treatment_col, variable_of_interest,
#                                'Effect Size', ['SN0212398523', 'DMSO', 'Leptomycin b'])
#
# # FIGURE 3 A
# plot_iqr_v_sample_size([*range(10, 500, 10)], 100, data_subset, treatment_col, variable_of_interest,
#                        'Error in Inter-Quartile Range')
#
# # FIGURE 3 B - H
# plot_cumulative_histogram_samples(data_subset, variable_of_interest, treatment_col, 'Untreated')

# FIGURE 4 A
# for s in [50, 200, 500, -1]:
#     generate_swarmplot_of_well_means(24, 10,
#                                      ['Untreated', 'DMSO', 'SN0212398523', 'SN1054616339', 'SN1066932540',
#                                       'Leptomycin b'],
#                                      ['Untreated', 'DMSO', 'SN0212398523', 'SN1054616339', 'SN1066932540',
#                                       'Leptomycin b'],
#                                      prepare_data(nuc_data, cyto_data, image_data, treatments, treatments_to_compounds,
#                                                   compounds,
#                                                   ['J05', 'I19', 'G15', 'O02', 'B02', 'N12', 'L08', 'L18', 'H13', 'E22',
#                                                    'H10', 'B06']), color_dict, treatment_col, variable_of_interest,
#                                      '$ \\log \\left[ \\frac {I_{F_N}}{(I_{F_N} + I_{C_N})} \\right]$', dunn_pairs, s)

print('All Done!')

#                                      ['Untreated', 'DMSO', 'SN0212398523', 'SN1054616339', 'SN0022071220',
#                                       'Leptomycin b'],
#                                      ['Untreated', 'DMSO', 'SN0212398523', 'SN1054616339', 'SN0022071220',
#                                       'Leptomycin b'],
