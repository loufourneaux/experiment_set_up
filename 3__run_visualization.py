import nibabel as nib
import numpy as np
import os 
import csv
import argparse
import plotly.graph_objects as go
from skimage import measure
from scipy.ndimage import label, binary_fill_holes
import scipy.ndimage
from tools.MRI_to_voxel import *
from tools.retrieve_coordinates import *

def get_biggest_component(binary_data):
    """
    Extracts the largest connected component from a binary image. 
    Used in our skull image so we only keep the skull and not potential artifacts.

    Args:
        binary_data (numpy.ndarray): A binary image where connected components
                                     are to be identified.

    Returns:
        numpy.ndarray: A binary image containing only the largest connected component.
    """
    # Label connected components in the binary image
    labeled_array, _ = label(binary_data)

    regions = measure.regionprops(labeled_array)

    # Find the region with the maximum area (biggest component)
    max_region = max(regions, key=lambda region: region.area)

    # Create a new binary image with only the largest component
    largest_component_binary = np.zeros_like(binary_data)
    largest_component_binary[tuple(max_region.coords.T)] = 1

    return largest_component_binary

def build_parser ():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="3D Plotting of Brain Areas")
    parser.add_argument("--subj", type=str, required=True, help="Subject identifier (e.g., 05)")
    parser.add_argument("--areas", type=str, required=True, help="Areas to visualize (m1, thalamus, both)")
    parser.add_argument("--surroundings", type=str, help="Surrounding areas to include, separated by commas")

    return parser


def brain_to_t1_space(t1_img, t1wBrain_img):
    """
    Aligns the brain-extracted T1-weighted image (T1wBrain aka brain wihtout skull) 
    to the original T1 space.

    This function computes the transformation matrix from the brain-extracted
    T1-weighted image (T1wBrain) to the original T1 space, applies the affine
    transformation, and returns the aligned brain data.

    Args:
        t1_img (nibabel.Nifti1Image): The original T1-weighted NIfTI image.
        t1wBrain_img (nibabel.Nifti1Image): The brain-extracted T1-weighted NIfTI image.

    Returns:
        numpy.ndarray: The brain data aligned to the original T1 space.
    """
    t1_data = t1_img.get_fdata()
    t1_affine = t1_img.affine

    t1wBrain_data = t1wBrain_img.get_fdata()
    t1wBrain_affine = t1wBrain_img.affine

    # Compute the transformation matrix from T1wBrain to T1
    transform_affine = np.linalg.inv(t1wBrain_affine) @ t1_affine

    # Apply the affine transformation
    aligned_brain_data = scipy.ndimage.affine_transform(
        t1wBrain_data,
        matrix=transform_affine[:3, :3],
        offset=transform_affine[:3, 3],
        output_shape=t1_data.shape
    )

    # Create a new NIfTI image with the aligned data
    #uncomment if you want to save the aligned brain image
    #aligned_brain_img = nib.Nifti1Image(aligned_brain_data, t1_affine)

    return aligned_brain_data


def plot_3D_figure(*args, directory):
    """
    Plot a 3D figure of all the args unpacked and save it in the directory path. 

    Args: 
        args (list): list of all the go objects that are going to be plot
        directory (path): place to store the figure

    Return:
        None, saves the image in the directory argument
    """
    fig = go.Figure(data=list(args))

    # Update layout for better visualization
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        ),
        scene_aspectmode='data'
    )

    #fig.show()
    fig.write_html(directory) #saves the figures in html
     

if __name__=='__main__':
    parser=build_parser()
    args = parser.parse_args()
    

    identifier = f'MOT{args.subj}'
    areas = args.areas.lower().split(',')
    surroundings = args.surroundings.split(',') if args.surroundings else []


    #Load the entire pCT for the plotting 
    ct= nib.load(f'data/PCTs/{identifier}_pCT.nii')
    ct_data= ct.get_fdata()
    
    #Load the t1 for space registration
    t1= nib.load(f'data/MOTs/{identifier}.nii')
    t1_data= t1.get_fdata()

    # Load the parcellation NIfTI file with the parcellated areas
    parcellation_img = nib.load(os.path.join( 'output', f'{identifier}','1__output', f'{identifier}_hcpmmp1_rein.nii.gz'))
    parcellation_data = parcellation_img.get_fdata()
    
    # Load the anatomical T1-weighted image for plotting 
    t1wBrain = nib.load(os.path.join('output', f'{identifier}','1__output', f'{identifier}T1wBrain.nii.gz'))
    t1wBrain_data = t1wBrain.get_fdata()

    # PROCESS THE PCT FOR HEAD PLOTTING 
    # Create binary image from 0 and other values image
    binary_data = (ct_data > -700).astype(int)
    binary_data = binary_fill_holes(binary_data)
    # Get the biggest component which is the skull using connected components (get rid of artifacts)
    binary_data=get_biggest_component(binary_data)
    # Extract the surface mesh using the Marching Cubes algorithm
    verts, faces, _, _ = measure.marching_cubes(binary_data, level=0.5) #returns a set of vertices and triangular faces

    # Register the brain to t1 space to plot and determine the areas' placement
    aligned_brain_data= brain_to_t1_space(t1, t1wBrain)
    verts_s, faces_s, _, _ =measure.marching_cubes(aligned_brain_data, level=0.5)

    #CREATE 3D OBJECTS FOR PLOTTING
    mesh3d = go.Mesh3d(
        x=verts[:, 0],
        y=verts[:, 1],
        z=verts[:, 2],
        i=faces[:, 0],
        j=faces[:, 1],
        k=faces[:, 2],
        opacity=0.3,
        color='gray',
    )
    mesh3d1 = go.Mesh3d(
        x=verts_s[:, 0],
        y=verts_s[:, 1],
        z=verts_s[:, 2],
        i=faces_s[:, 0],
        j=faces_s[:, 1],
        k=faces_s[:, 2],
        opacity=0.3,
        color='gray',
    )
    #all the points to be scattered 
    scatters = [mesh3d, mesh3d1]

    #if add areas then need to add the areas to the color dict as well. 
    color_dict = {       
        'hippocampus': 'blue',                
        'caudate': 'green',        
        'ventralDC': 'chocolate',                
        'brainstem': 'cyan',  
        'SMAa': 'mistyrose',
        'SMAp': 'olive',
        'SCEF': 'sienna', 
        'S1': 'gold', 
        'anterior6': 'salmon',
        'dorsal6': 'orchid',
        'rostral6':'lightblue',
        'ventral6':'peru'      
    }
    
    # Plot surrounding areas if specified
    for sur in surroundings:
        if sur=='brainstem':
            coords = read_coordinates(os.path.join('results', 'coordinates_voxel', 'all_coordinates.csv'), identifier,sur)
            scatter = go.Scatter3d(
                x=[point[0] for point in coords],
                y=[point[1] for point in coords],
                z=[point[2] for point in coords],
                mode='markers',
                marker=dict(size=3, color=color_dict[sur]),
                name=sur
            )
            scatters.append(scatter)
        else:
            right_sur=f'R_{sur}'
            left_sur=f'L_{sur}'
            coords = read_coordinates(os.path.join('output',f'{identifier}/2__output', 'coordinates_voxel', 'all_coordinates.csv'), identifier,left_sur)
            coords2= read_coordinates(os.path.join('output',f'{identifier}/2__output', 'coordinates_voxel', 'all_coordinates.csv'), identifier,right_sur)
            scatter = go.Scatter3d(
                x=[point[0] for point in coords],
                y=[point[1] for point in coords],
                z=[point[2] for point in coords],
                mode='markers',
                marker=dict(size=3, color=color_dict[sur]), 
                showlegend=False
            )
            scatter2 = go.Scatter3d(
                x=[point[0] for point in coords2],
                y=[point[1] for point in coords2],
                z=[point[2] for point in coords2],
                mode='markers',
                marker=dict(size=3, color=color_dict[sur]),
                name=sur
            )
            scatters.append(scatter)
            scatters.append(scatter2)

    if 'm1' in areas:
        coords_l=read_coordinates(os.path.join('output',f'{identifier}/2__output', 'coordinates_voxel', 'all_coordinates.csv'), identifier, 'L_m1')
        coords_r=read_coordinates(os.path.join('output',f'{identifier}/2__output', 'coordinates_voxel', 'all_coordinates.csv'), identifier, 'R_m1')
        left_m1=go.Scatter3d(
            x=[point[0] for point in coords_r],
            y=[point[1] for point in coords_r],
            z=[point[2] for point in coords_r],
            mode='markers',
            marker=dict(size=3, color='red'),
            showlegend=False
        )

        right_m1=go.Scatter3d(
            x=[point[0] for point in coords_l],
            y=[point[1] for point in coords_l],
            z=[point[2] for point in coords_l],
            mode='markers',
            marker=dict(size=3, color='red'),
            name='M1'
        )
        scatters.append(left_m1)
        scatters.append(right_m1)
    
    if 'thalamus' in areas :
        coords_l=read_coordinates(os.path.join('output',f'{identifier}/2__output', 'coordinates_voxel', 'all_coordinates.csv'), identifier, 'L_thalamus')
        coords_r=read_coordinates(os.path.join('output',f'{identifier}/2__output', 'coordinates_voxel', 'all_coordinates.csv'), identifier, 'R_thalamus')
        left=go.Scatter3d(
            x=[point[0] for point in coords_r],
            y=[point[1] for point in coords_r],
            z=[point[2] for point in coords_r],
            mode='markers',
            marker=dict(size=3, color='orange'),
            showlegend=False
        )

        right=go.Scatter3d(
            x=[point[0] for point in coords_l],
            y=[point[1] for point in coords_l],
            z=[point[2] for point in coords_l],
            mode='markers',
            marker=dict(size=3, color='orange'),
            name='thalamus'
        )
        scatters.append(left)
        scatters.append(right)

    # Plot the figure
    plot_3D_figure(*scatters, directory=f'output/{identifier}/3__output/{areas}_{surroundings}_visualization.html')