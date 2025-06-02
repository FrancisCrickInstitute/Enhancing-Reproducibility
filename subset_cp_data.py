from notebooks.utility_functions import *

plate_IDs = ['1A', '2A', '2B']
wells_per_plate = {'1A': ['I23', 'J23', 'M23', 'N23', 'B06'], '2A': ['C13', 'J23', 'N23', 'K16'],
                   '2B': ['C13', 'C17', 'P13']}
all_data = pd.DataFrame()
for plate_ID in plate_IDs:
    plate_name = f'LM2_GEFGAP_ONTARGETPlus_{plate_ID}'
    plate_index = f'LM2_ONTARGETPlus_{plate_ID}'
    idr_annotations_file_path = 'inputs/idr/idr0028-screenB-annotation.csv'
    image_data = pd.read_csv(f'inputs/cell_profiler_outputs/idr0028/screenB/{plate_name}/Image.csv')
    nuc_data = pd.read_csv(f'inputs/cell_profiler_outputs/idr0028/screenB/{plate_name}/Nuclei.csv')
    cyto_data = pd.read_csv(f'inputs/cell_profiler_outputs/idr0028/screenB/{plate_name}/Cytoplasm.csv')
    image_indices = pd.read_csv(f'inputs/idr/{plate_name}_ImageIndex.ColumbusIDX.csv', delimiter='\t')
    well_to_filename = dict(zip(image_indices['WellName'], image_indices['sourcefilename']))
    filename_to_image_number = dict(zip(image_data['FileName_Hoechst'], image_data['ImageNumber']))
    well_to_image_number = []
    for well in wells_per_plate.get(plate_ID):
        # Step 1: Map well to filename
        filename = well_to_filename.get(well)
        # Step 2: Map filename to image number
        if filename is not None:
            image_number = filename_to_image_number.get(filename)
            well_to_image_number.append(image_number)

    output_dir = f'data_subsets/cell_profiler_outputs/idr0028/screenB/{plate_name}'
    os.makedirs(output_dir, exist_ok=True)
    image_subset = image_data[image_data['ImageNumber'].isin(well_to_image_number)]
    nuc_subset = nuc_data[nuc_data['ImageNumber'].isin(well_to_image_number)]
    cyto_subset = cyto_data[cyto_data['ImageNumber'].isin(well_to_image_number)]
    image_subset.to_csv(f'{output_dir}/Image.csv')
    nuc_subset.to_csv(f'{output_dir}/Nuclei.csv')
    cyto_subset.to_csv(f'{output_dir}/Cytoplasm.csv')
