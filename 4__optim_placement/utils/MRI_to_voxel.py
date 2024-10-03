import numpy as np
import math
from utils.MRI_to_voxel import *


def mri_to_voxel(file, target_coordinates):
    """
    Convert MRI coordinates to voxel coordinates using the inverse affine transformation matrix.

    Args:
        file (nib.Nifti1Image): The NIfTI image file containing the affine transformation matrix.
        target_coordinates (np.ndarray): The target coordinates in MRI space (RAS).

    Returns:
        np.ndarray: The corresponding voxel coordinates, rounded to the nearest integers.
    """
    #Get the affine transformation
    affine=file.affine
    inverse_affine = np.linalg.inv(affine)

    #get voxel coordinates 
    voxel_coords = inverse_affine @ target_coordinates
    voxel_coords=np.round(voxel_coords).astype(int)

    return voxel_coords[:3]

def voxel_to_mri(file, voxel_coordinates):
    """
    Convert voxel coordinates to MRI coordinates using the affine transformation matrix.

    Args:
        file (nib.Nifti1Image): The NIfTI image file containing the affine transformation matrix.
        voxel_coordinates (np.ndarray): The voxel coordinates

    Returns:
        np.ndarray: The corresponding MRI coordinates.
    """

    # Get the affine transformation
    affine = file.affine

    # Convert voxel coordinates to homogeneous coordinates
    voxel_homogeneous = np.append(voxel_coordinates, 1)

    # Apply the affine transformation
    mri_coords = affine @ voxel_homogeneous

    return mri_coords[:3]

def physical_size_to_voxels(physical_size_mm, voxel_sizes_mm):
    """
    Convert physical size in millimeters to voxel size, assuming isotropic voxel sizes.

    Args:
        physical_size_mm (float): The physical size in millimeters to be converted.
        voxel_sizes_mm (np.ndarray): The voxel sizes in millimeters for each dimension. Assumes isotropic voxels.

    Returns:
        int: The corresponding voxel size.
    """
    voxel_size_mm = voxel_sizes_mm[0]  # Assuming all same size in every direction
    voxel_size = int(physical_size_mm / voxel_size_mm)
    return voxel_size
    
def voxels_to_physical_size(voxel_size, voxel_sizes_mm):
    """
    Convert voxel size to physical size in millimeters, assuming isotropic voxel sizes.

    Args:
        voxel_size (int): The size in voxels to be converted.
        voxel_sizes_mm (np.ndarray): The voxel sizes in millimeters for each dimension. Assumes isotropic voxels.

    Returns:
        float: The corresponding physical size in millimeters.
    """
    voxel_sizes_mm=voxel_sizes_mm[0]
    physical_size_mm= voxel_size*voxel_sizes_mm
    return physical_size_mm