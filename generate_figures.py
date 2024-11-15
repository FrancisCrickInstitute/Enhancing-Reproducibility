from utility_functions import *

plt.rcParams['font.size'] = 26
plt.rcParams['axes.linewidth'] = 2
plate_number = 1093711385
treatment_col = 'Treatment'
variable_of_interest = 'Fascin_Ratio'
color_dict = {'SN1066932540': 'red', 'SN1054616339': 'yellow', 'SN0212398523': 'orange', 'Untreated': 'blue',
              'DMSO': 'gray', 'Leptomycin b': 'purple'}

directories = ('./inputs/idr', './outputs/plots', './outputs/data')
y_label = 'Relative Nuclear Fascin Localisation'

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
output_dir = './outputs/plots'

data_subset = prepare_data(nuc_data, cyto_data, image_data, treatments, treatments_to_compounds, compounds,
                           ['J05', 'O02', 'E22', 'L08'])

# FIGURE 2 A - F
point_size = 8
filenames = ['Fig2A.png', 'Fig2B.png', 'Fig2C.png', 'Fig2D.png', 'Fig2E.png', 'Fig2F.png']
filecount = 0
random_seed = 42
for s in [50, 200]:
    if s > 50:
        point_size = 4
    for i in range(3):
        generate_swarmplot(['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'], data_subset, color_dict,
                           treatment_col, variable_of_interest, y_label, os.path.join(output_dir, filenames[filecount]),
                           point_size=point_size, random_seed=random_seed, sample_size=s)
        filecount = filecount + 1
        random_seed = random_seed + 1

# FIGURE 2 G - I
filenames = ['Fig2G.png', 'Fig2H.png', 'Fig2I.png']
plot_effect_size_v_sample_size([*range(10, 500, 10)], 100, data_subset, treatment_col, variable_of_interest,
                               'Median Effect Size Relative to Untreated', ['SN0212398523', 'DMSO', 'Leptomycin b'],
                               output_dir, filenames)

# FIGURE 3 A
plot_iqr_v_sample_size([*range(10, 500, 10)], 100, data_subset, treatment_col, variable_of_interest,
                       'Error in Inter-Quartile Range', os.path.join(output_dir, 'Fig3A.png'))

# FIGURE 3 B - H
filenames = ['Fig3B.png', 'Fig3C.png', 'Fig3D.png', 'Fig3E.png', 'Fig3F.png', 'Fig3G.png', 'Fig3H.png']
plot_cumulative_histogram_samples(data_subset, variable_of_interest, treatment_col, 'Untreated', output_dir, filenames,
                                  y_label)

# FIGURE 4 A - C
filenames = ['Fig4A.png', 'Fig4B.png', 'Fig4C.png']
filecount = 0
selected_treatments = ['Untreated', 'DMSO', 'SN0212398523', 'SN1054616339', 'SN1066932540', 'Leptomycin b']
point_size = {50: 8, 200: 4, 500: 2.5}
for s in [50, 200, 500]:
    generate_superplot(selected_treatments, selected_treatments,
                       prepare_data(nuc_data, cyto_data, image_data, treatments, treatments_to_compounds,
                                                  compounds,
                                                  ['J05', 'I19', 'G15', 'O02', 'B02', 'N12', 'L08', 'L18', 'H13', 'E22',
                                                   'H10', 'B06']), color_dict, treatment_col, variable_of_interest,
                       y_label, os.path.join(output_dir, filenames[filecount]), sample_size=s,
                       point_size=point_size[s])
    filecount = filecount + 1

# SUPP FIGURE 1
generate_swarmplot(['Untreated', 'DMSO', 'SN0212398523', 'Leptomycin b'], data_subset, color_dict, treatment_col,
                   variable_of_interest, y_label, os.path.join(output_dir, 'SupFig1.png'))

print('All Done!')
