import numpy as np


def mri_to_voxel(file, target_coordinates):
    '''
    Get target voxel coordinates based on MRI coordinates
    '''

    #Get the affine transformation
    affine=file.affine
    inverse_affine = np.linalg.inv(affine)

    #get voxel coordinates 
    voxel_coords = inverse_affine @ target_coordinates

    return voxel_coords[:3]

def voxel_to_mri(file, voxel_coordinates):
    '''
    Get target MRI coordinates based on voxel coordinates
    '''

    # Get the affine transformation
    affine = file.affine

    # Convert voxel coordinates to homogeneous coordinates
    voxel_homogeneous = np.append(voxel_coordinates, 1)

    # Apply the affine transformation
    mri_coords = affine @ voxel_homogeneous

    return mri_coords

def physical_size_to_voxels(physical_size_cm, voxel_sizes_mm):
    '''
    Convert physical size in cm to voxel size.
    Assumes isotropic voxel sizes (i.e., all dimensions have the same size).
    '''
    voxel_size_mm = voxel_sizes_mm[0]  # Assuming all same size in every direction
    voxel_size = int(physical_size_cm * 10 / voxel_size_mm)
    return voxel_size
    
