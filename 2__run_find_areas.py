import nibabel as nib
import os 
import csv
import numpy as np
from skimage import measure
from scipy.ndimage import label
from tools.MRI_to_voxel import *
import tqdm


def get_coordinates(parcellation_data, label):
    """
    This function extracts the coordinates in voxels of the area labelled with 'label' 

    Returns:
        numpy.ndarray: all the coordinates
    """
    return np.argwhere(parcellation_data == label)

def get_center_of_mass(coords):
    """
    Determines the center of mass af a given area determined by a list of points (coords).
    
    Args:
        coords (list): a list of 3D points
    """
    x_coords, y_coords, z_coords = zip(*coords)
    # Calculating the means
    x_mean = sum(x_coords) / len(x_coords)
    y_mean = sum(y_coords) / len(y_coords)
    z_mean = sum(z_coords) / len(z_coords)

    # Center of mass
    center_of_mass = [[x_mean, y_mean, z_mean]]
    return center_of_mass

def read_existing_data(file_path):
    """
    Reads the existing data in the all_coordinates.csv file so we don't write the same information twice in the result csv file 
    """
    existing_data = set()
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                # Convert each row to a tuple and add to the set
                existing_data.add(tuple(row))
    except FileNotFoundError:
        # File not found, so we assume no existing data
        pass
    return existing_data

def write_results_to_file(file_path, identifier, coordinates, coord_type):
    """
    Writes the coordinates to a csv file if the coordinates are not already present in the file
    """
    # Read existing rows
    existing_rows = read_existing_data(file_path)
    
    new_rows = []
    
    for coord in coordinates:
        # Prepare the new row as a list
        row = [identifier, coord_type, *coord]
        # Convert the row to a tuple of strings for comparison
        row_str = tuple(map(str, row))
        if row_str not in existing_rows:
            new_rows.append(row)
    
    if new_rows:
        with open(file_path, 'a', newline='') as csvfile:  # Open in append mode
            writer = csv.writer(csvfile)
            writer.writerows(new_rows)

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

if __name__=='__main__':
    for file in os.listdir('data/MOTs'):
        filename = os.fsdecode(file)
        if filename.endswith('.nii') or filename.endswith('.nii.gz'): #ensure only processes files with .nii extension
            #Load the entire T1 
            t1 = nib.load(os.path.join('data/MOTs',filename))
            t1_data= t1.get_fdata()
            # Extract the identifier (MOTXX) from the filename (MOTXX.nii)
            identifier = filename.split('.')[0]  # Assumes the identifier is before the first underscore

            # Load the parcellation NIfTI file to find the areas 
            parcellation_img = nib.load(os.path.join( 'output' , f'{identifier}', '1__output',f'{identifier}_hcpmmp1_rein.nii.gz'))
            parcellation_data = parcellation_img.get_fdata()
        
            # Load the anatomical T1-weighted image to transform the coordinates 
            t1wBrain = nib.load(os.path.join('output', f'{identifier}','1__output', f'{identifier}T1wBrain.nii.gz'))
            t1wBrain_data = t1wBrain.get_fdata()

            #M1 and thalamus both hemispheres 
            # add the label id and name to the dictionary to find other areas --> go see in the ordered.txt file 

            label_values = {
                'L_hippocampus': 366,
                'R_hippocampus': 375,
                'L_caudate': 363,
                'R_caudate': 372,
                'L_ventralDC': 369,
                'R_ventralDC': 378,
                'brainstem': 379,
                'L_m1': 8, #label L_4
                'L_thalamus': 362,
                'R_m1': 188,
                'R_thalamus': 371, 
                'L_S1':9, #label 3b
                'R_S1':189,
                'L_SMAa': 44,#label 6ma
                'R_SMAa':224,
                'L_SMAp': 55, #label 6mp
                'R_SMAp':235,
                'L_SCEF':43, #label SCEF
                'R_SCEF': 223,
                'L_dorsal6': 54, #label 6dorsal
                'R_dorsal6': 234,
                'L_ventral6': 56, #label 6ventral
                'R_ventral6': 236,
                'L_rostral6': 78, #label 6rostral
                'R_rostral6': 258,
                'L_anterior6': 96, #label 6anterior
                'L_anterior6': 276
            }
            # Define file paths for MRI and voxel coordinates
            mri_results_file = f'output/{identifier}/2__output/coordinates_mri/all_coordinates.csv'
            ensure_directory_exists(os.path.dirname(mri_results_file))
            voxel_results_file = f'output/{identifier}/2__output/coordinates_voxel/all_coordinates.csv'
            ensure_directory_exists(os.path.dirname(voxel_results_file))
            # Define file paths for MRI and voxel center of mass coordinates
            cm_mri_results_file = f'output/{identifier}/2__output/coordinates_mri/cm_coordinates.csv'
            ensure_directory_exists(os.path.dirname(cm_mri_results_file))
            cm_voxel_results_file = f'output/{identifier}/2__output/coordinates_voxel/cm_coordinates.csv'
            ensure_directory_exists(os.path.dirname(cm_voxel_results_file))

            for label in label_values:
                coords= get_coordinates(parcellation_data, label_values[label])
                #Turn these coordinates back to MRI coordinates 
                coords_mri_final= [voxel_to_mri(t1wBrain, coord) for coord in coords]
                coords_mri= [voxel_to_mri(t1wBrain, coord)[:3] for coord in coords]
                write_results_to_file(mri_results_file, identifier, coords_mri, label)
                
                #Save the center of mass of these areas 
                center_mass_mri=get_center_of_mass(coords_mri)
                write_results_to_file(cm_mri_results_file, identifier, center_mass_mri, label)

                
                #turn into voxel coordinates back in T1 space
                final_coords = [mri_to_voxel(t1, coord_mri) for coord_mri in coords_mri_final]
                write_results_to_file(voxel_results_file, identifier, final_coords, label)
                
                #Save the center of mass of these areas 
                center_mass_voxel=get_center_of_mass(final_coords)
                write_results_to_file(cm_voxel_results_file, identifier, center_mass_voxel, label)
            
            