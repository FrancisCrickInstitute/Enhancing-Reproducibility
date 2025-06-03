import pandas as pd
import os
from collections import defaultdict

plate_IDs = ['1A', '2A', '2B']
wells_per_plate = {'1A': ['I23', 'J23', 'M23', 'N23', 'B6'], '2A': ['C13', 'J23', 'N23', 'K16'],
                   '2B': ['C13', 'C17', 'P13']}
idr_output_dir = f'data_subsets/idr'
os.makedirs(idr_output_dir, exist_ok=True)
idr_annotations = pd.read_csv('inputs/idr/idr0028-screenB-annotation.csv')
idr_annotations_subset_list = []
for plate_ID in plate_IDs:
    plate_name = f'LM2_GEFGAP_ONTARGETPlus_{plate_ID}'
    plate_index = f'LM2_ONTARGETPlus_{plate_ID}'
    image_data = pd.read_csv(f'inputs/cell_profiler_outputs/idr0028/screenB/{plate_name}/Image.csv')
    nuc_data = pd.read_csv(f'inputs/cell_profiler_outputs/idr0028/screenB/{plate_name}/Nuclei.csv')
    cyto_data = pd.read_csv(f'inputs/cell_profiler_outputs/idr0028/screenB/{plate_name}/Cytoplasm.csv')
    image_indices = pd.read_csv(f'inputs/idr/{plate_name}_ImageIndex.ColumbusIDX.csv', delimiter='\t')
    # Assuming image_indices is your DataFrame
    well_to_filenames = defaultdict(list)
    # Iterate over the DataFrame and populate the dictionary
    for well_name, source_filename in zip(image_indices['WellName'], image_indices['sourcefilename']):
        well_to_filenames[well_name].append(source_filename)
    # Convert the defaultdict to a regular dict if needed
    well_to_filenames = dict(well_to_filenames)
    filename_to_image_number = dict(zip(image_data['FileName_Hoechst'], image_data['ImageNumber']))
    image_numbers = []
    for well in wells_per_plate.get(plate_ID):
        # Step 1: Map well to filename
        filenames = well_to_filenames.get(well)
        # Step 2: Map filename to image number
        for filename in filenames:
            image_number = filename_to_image_number.get(filename)
            image_numbers.append(image_number)

    cp_output_dir = f'data_subsets/cell_profiler_outputs/idr0028/screenB/{plate_name}'
    os.makedirs(cp_output_dir, exist_ok=True)
    image_subset = image_data[image_data['ImageNumber'].isin(image_numbers)]
    nuc_subset = nuc_data[nuc_data['ImageNumber'].isin(image_numbers)]
    cyto_subset = cyto_data[cyto_data['ImageNumber'].isin(image_numbers)]
    image_subset.to_csv(f'{cp_output_dir}/Image.csv')
    nuc_subset.to_csv(f'{cp_output_dir}/Nuclei.csv')
    cyto_subset.to_csv(f'{cp_output_dir}/Cytoplasm.csv')

    image_indices_subset = image_indices[image_indices['WellName'].isin(wells_per_plate.get(plate_ID))]
    image_indices_subset.to_csv(f'{idr_output_dir}/{plate_name}_ImageIndex.ColumbusIDX.csv', sep='\t')

    idr_annotations_subset = idr_annotations[idr_annotations['Plate'] == plate_index]
    idr_annotations_subset = idr_annotations_subset[idr_annotations_subset['Well'].isin(wells_per_plate.get(plate_ID))]
    idr_annotations_subset_list.append(idr_annotations_subset)

pd.concat(idr_annotations_subset_list).to_csv(f'{idr_output_dir}/idr0028-screenB-annotation.csv')
