import os
import sys 
import pandas as pd
import numpy as np
from tqdm import tqdm
from skimage import measure
import plotly.graph_objects as go
from scipy.spatial import KDTree
from utils.MRI_to_voxel import *
from scipy.ndimage import label

def create_sphere(file, surface_points, voxel_coord, interactive, radius):
    """
    This function determines the appropriate radius for the search sphere based on user input.
    It calculates the Euclidean distance from surface points to the
    target voxel coordinates and creates a mask to identify points within the sphere. 

    Args:
        file (nibabel.Nifti1Image): The NIfTI file containing the voxel dimensions information.
        surface_points (numpy.ndarray): The coordinates of all the surface points.
        voxel_coord (numpy.ndarray): The target voxel coordinates.
        interactive (bool): Whether to prompt the user for input interactively.
        radius (float): The radius of the sphere in millimeters.

    Returns:
        tuple: A tuple containing:
            - search_area_points (numpy.ndarray): The points within the sphere.
            - sphere_scatter (plotly.graph_objects.Surface): The Plotly surface object for the sphere.
            - message (str): A message describing the radius used for the search.
    """
    # Define the radius of the sphere
    min_radius= smallest_distance(file, surface_points, voxel_coord)
    search_area_points=[]
    while len(search_area_points)==0:
        if interactive==True:
            rad = float(input(f"Enter the radius of the sphere in which you want to search for positions (in mm), the minimum distance is {math.ceil(min_radius)}: ")) 
            message = f'The radius used is the one given as an input: {rad} mm '
        elif interactive==False:
            if radius < math.ceil(min_radius):
                rad= math.ceil(min_radius)
                message = f'There was no surface points in the given radius so it was set the the minimum distance of {math.ceil(min_radius)}'
            else :
                rad = radius
                message = f'The radius used is the one given as an input: {radius} mm '
        radius_voxel = physical_size_to_voxels(rad, file.header.get_zooms())
        # Calculate the Euclidean distance from each surface point to the target coordinates
        distances=np.linalg.norm(surface_points-voxel_coord, axis=1)
        # Create a mask for points inside the sphere
        sphere_mask = distances <= radius_voxel
        # Extract the points inside the sphere
        search_area_points = surface_points[sphere_mask]

        # Define the mesh grid for the sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = voxel_coord[0] + radius_voxel * np.cos(u) * np.sin(v)
        y = voxel_coord[1] + radius_voxel * np.sin(u) * np.sin(v)
        z = voxel_coord[2] + radius_voxel * np.cos(v)
        
        # Create the plotly figure
        sphere_scatter = go.Surface(x=x, y=y, z=z, opacity=0.3, colorscale='Viridis', showscale=False)
    
    return search_area_points, sphere_scatter, message

def smallest_distance(file, surface_points, coord):
    """
    Determines the smallest distance from the coord to a surface point for when creating the sphere of search area. Can't enter
    a radius smaller than the smallest distance otherwise the search area will be empty. 
    """
    #determine all the distances 
    distances=np.linalg.norm(surface_points - coord, axis=1)
    min_dist = np.min(distances) #in voxels 
    return voxels_to_physical_size(min_dist, file.header.get_zooms())



def plot_3D_figure(*args, directory):
    """
    Plot a 3D figure of all the args unpacked and save it in the directory path. 

    Args: 
        args (list): list of all the go objects that are going to be plotted.
        directory (path): place to store the figure.

    Return:
        None, saves the image in the directory argument.
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
    fig.write_html(directory)

def create_transducer_plot(normal, intersect_point, file):
    """
    Create a tangent disk plot at the optimal placement coordinate.

    Args:
        normal (numpy.ndarray): the normal of the tangent plane.
        intersect_point(numpy.ndarray): point around which create the transducer point.
        file (nib.Nifti1Image): contains infromation about voxel dimensions.
    
    Returns:
        go.Surface: A Plotly Surface object representing the tangent plane.
    """
    #transducer shape
    diameter=physical_size_to_voxels(56.32, file.header.get_zooms())
    # Compute two orthogonal vectors on the plane
    if np.allclose(normal, [1, 0, 0]):
        arbitrary_vector = np.array([0, 1, 0])
    else:
        arbitrary_vector = np.array([1, 0, 0])

    u = np.cross(normal, arbitrary_vector)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    v /= np.linalg.norm(v)

    # Define the radius
    radius = diameter / 2

    # Create points within the disk
    theta = np.linspace(0, 2 * np.pi, 100)
    r = np.linspace(0, radius, 50)
    r, theta = np.meshgrid(r, theta)
    x_circle = r * np.cos(theta)
    y_circle = r * np.sin(theta)

    # Transform these points onto the plane
    xx = intersect_point[0]  + x_circle * u[0] + y_circle * v[0]
    yy = intersect_point[1]  + x_circle * u[1] + y_circle * v[1]
    zz = intersect_point[2]  + x_circle * u[2] + y_circle * v[2]

    # Create a surface plot for the tangent plane
    tangent_plane = go.Surface(x=xx, y=yy, z=zz, colorscale='Blues', opacity=0.9, showscale=False)
    return tangent_plane

def get_biggest_component(binary_data):
    """
    Extracts the largest connected component from a binary image. We only want to plot the scalp's surface 
    and not potential artifacts extracted from the binary image. 

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



def cosine_similarity(vec1, vec2):
    """
    Compute the cosine similarity between two vectors.

    Args:
        vec1 (numpy.ndarray): The first vector.
        vec2 (numpy.ndarray): The second vector.

    Returns:
        float: The cosine similarity between vec1 and vec2, which indicates how close the vectors are to being parallel.
               A value close to -1 or 1 indicates that the vectors are close to parallel.
    """
    dot_product = np.dot(vec1, vec2)
    norm_vec1 = np.linalg.norm(vec1)
    norm_vec2 = np.linalg.norm(vec2)
    return dot_product / (norm_vec1 * norm_vec2)


def tangent_to_surface(verts, intersect_point, radius=30):
    '''
    Determine the tangent to the local surface 

    Inputs:
        verts: represent the surface
        intersect_point (numpy.ndarray): the point around which determine the surface and tangent

    Return:
        normal (numpy.ndarray): normal to the tangent plane 
    '''
    # Create a KDTree for fast nearest neighbor search
    kdtree = KDTree(verts)

    # Find all vertices within the specified radius
    indices = kdtree.query_ball_point(intersect_point, r=radius)
    local_verts = verts[indices]

    # Fit a plane to the local vertices using Singular Value Decomposition (SVD)
    # Center the local vertices
    centroid = local_verts.mean(axis=0)
    centered_verts = local_verts - centroid

    # Perform SVD
    U, S, Vt = np.linalg.svd(centered_verts)
    normal = Vt[2, :]  # The normal to the plane is the last row of Vt

    return normal

def point_translation(intersect_point, voxel_coord, normal, length, file):
    """
    Move the transducer placement along the direction of the beam to get the right coordinates for neural navigation and simulation.

    Args:
        intersect_point (numpy.ndarray): The initial intersection point coordinates.
        voxel_coord (numpy.ndarray): The target voxel coordinates.
        normal (numpy.ndarray): The normal vector at the intersection point.
        length (float): The distance to move along the normal direction.
        file (nib.Nifti1Image): The NIfTI file object containing the voxel dimensions information.

    Returns:
        numpy.ndarray: The new intersection point after translation.
    """
    
    length=physical_size_to_voxels(length, file.header.get_zooms())
    
    # Calculate the initial distance before moving
    dist_bf_mv = np.linalg.norm(intersect_point - voxel_coord)
    
    # Translate the intersect point along the normal direction by 'length'
    new_intersect = intersect_point - normal * length

    # Calculate the distance after moving
    dist_after_mv = np.linalg.norm(new_intersect - voxel_coord)

    # If the distance decreased, move in the opposite direction bc that means we're inside the skull
    if dist_after_mv < dist_bf_mv:
        new_intersect = intersect_point + normal * length

    return new_intersect


def determine_thickness(binary_image, voxel_coord, end_point, direction):
    """
    Traverse a line from voxel_coord to end_point in the binary image and calculate the thickness.

    Parameters:
        binary_image (numpy.ndarray): 3D numpy array representing the binary image (1 for skull, 0 for background)
        voxel_coord (numpy.ndarray): Starting voxel coordinate 
        end_point (numpy.ndarray): End point voxel coordinate
        direction (numpy.ndarray): Direction vector from voxel_coord to end_point.

    Returns:
        thickness (float): The thickness of the skull along the line
    """
    # Convert coordinates to numpy arrays for easy calculations
    voxel_coord = np.array(voxel_coord, dtype=np.float64)
    end_point = np.array(end_point, dtype=np.float64)

    length = np.linalg.norm(direction)
    direction = direction / length  # Normalize the direction vector

    # Initialize variables to track the first and last skull points
    first_skull_point = None
    last_skull_point = None
    current_point = voxel_coord.copy()

    # Traverse the line
    while True:
        # Round the current point to the nearest voxel index
        voxel_index = np.round(current_point).astype(int)

        # Check if the voxel index is within the image boundaries
        if (0 <= voxel_index[0] < binary_image.shape[0] and
            0 <= voxel_index[1] < binary_image.shape[1] and
            0 <= voxel_index[2] < binary_image.shape[2]):
            
            # Check if the voxel is part of the skull
            if binary_image[voxel_index[0], voxel_index[1], voxel_index[2]] == 1:
                if first_skull_point is None:
                    first_skull_point = current_point.copy()  # Set the first skull point
                last_skull_point = current_point.copy()  # Update the last skull point
                
        else:
            # Exit if the voxel index is out of bounds
            break

        # Move to the next voxel along the line
        current_point += direction

        # Break if we have reached the end point
        if np.linalg.norm(current_point - voxel_coord) >= length:
            break

    if first_skull_point is not None and last_skull_point is not None:
        # Calculate the thickness as the Euclidean distance between the first and last skull points
        thickness = np.linalg.norm(last_skull_point - first_skull_point)
    else:
        # If no skull points are found, thickness is zero
        thickness = 0

    return thickness

def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def write_results_to_file(file, file_path, target_coord, voxel_coord, search_area_points,search_area, best_placements, time, mode, transducer_coord,transducer_model,  message, gel_pad):
    """
    Writes the results of the placement analysis to a CSV file and appends additional information.

    This function creates a DataFrame from the placement results and saves it to a CSV file. 
    It also appends various details about the input parameters, transducer settings, and script execution time.

    Args:
        file (nibabel.Nifti1Image): The NIfTI file containing the image data.
        file_path (str): Path to the CSV file where results will be saved.
        target_coord (list): Target coordinates in RAS space.
        voxel_coord (list): Target coordinates in voxel space.
        search_area_points (numpy.ndarray): Points within the search area sphere.
        search_area (str): Description of the search area.
        best_placements (list): List of best placements with details such as coordinates and scores.
        time (float): Time taken by the script to run, in seconds.
        mode (int): Mode of operation (e.g., 1 or 2).
        transducer_coord (list): Transducer coordinates in RAS space.
        transducer_model (str): Model of the transducer used.
        message (str): Message related to the radius or search area.
        gel_pad (float): Thickness of the gel pad in mm.

    Returns:
        None
    """
    # Initialize lists to store the coordinates and scores
    exit_plane_voxels = []
    exit_plane_RAS = []
    center_of_mass_voxels = []
    center_of_mass_RAS = []
    similarity_scores = []
    distances=[]
    thicknesses=[]
    
    for placement in best_placements:
        exit_plane_voxel = placement[1]
        exit_plane_ras = voxel_to_mri(file, placement[1])
        center_of_mass_voxel = placement[4]
        center_of_mass_ras = voxel_to_mri(file, placement[4])
        score = np.round(np.abs(placement[3]),4)
        distance =np.round(placement[5],2)
        thickness =np.round(placement[6],3)
        
        exit_plane_voxels.append(exit_plane_voxel)
        exit_plane_RAS.append(exit_plane_ras)
        center_of_mass_voxels.append(center_of_mass_voxel)
        center_of_mass_RAS.append(center_of_mass_ras)
        similarity_scores.append(score)
        distances.append(distance)
        thicknesses.append(thickness)
    
    # Create a DataFrame to hold the data
    df = pd.DataFrame({
        'Exit Plane Voxel Coordinates': exit_plane_voxels,
        'Exit Plane RAS Coordinates': exit_plane_RAS,
        'Center of Mass Voxel Coordinates': center_of_mass_voxels,
        'Center of Mass RAS Coordinates': center_of_mass_RAS,
        'Similarity Scores': similarity_scores, 
        'Distances': distances,
        'Thicknesses': thicknesses
    })
    
    # Save the DataFrame to a CSV file
    df.index = df.index + 1
    df.to_csv(file_path, index=True)

    with open(file_path, 'a') as f:
        f.write("\n")
        f.write(f"The given inputs were:\n")
        f.write(f"Area of choice: {search_area}\n")
        f.write(f"Area coordinates: {target_coord[:3]} in RAS and in voxels: {voxel_coord} \n")
        f.write(f"Transducer model: {transducer_model}\n")
        if mode ==2:
            f.write(f"Transducer coordinates: {transducer_coord[:3]} in RAS and in voxels: {mri_to_voxel(file, transducer_coord)}\n")
        f.write(f"Gel pad thickness: {gel_pad}\n")
        if mode ==1 or mode ==2:
            f.write(f'{message}\n')
        f.write("\n")
        f.write(f"Searching through {len(search_area_points)} possible placement positions\n")
        f.write(f"The script ran for {time} seconds\n")



def find_optimal_placement( binary_image, verts, voxel_coord, closest_vertices, file, transducer_mode, gel_pad, nb_points, radius):
    """
    Test how similar the tangent planes are for each closest point, then filter by skull thickness and distance.

    Parameters:
        binary_image (numpy.ndarray): Binary image data.
        verts (numpy.ndarray): Array of vertices from the surface mesh.
        voxel_coord (numpy.ndarray): Voxel coordinates of the target point.
        closest_vertices (numpy.ndarray): Array of the closest vertices to the target point.
        file (Nifti1Image): NIfTI image file object.
        transducer_mode (str): Transducer model identifier (e.g., 'CTX250', 'CTX500').
        gel_pad (float): Thickness of the gel pad in mm.
        nb_points (int): Number of points to consider for optimal placement.
        radius (float): Radius for tangent plane calculation.

    Returns:
        list: List of tuples containing the tangent plane, optimal point, direction, similarity score,
              simulated intersection point, distance, and thickness for each optimal placement.
    """
    similarities = np.array([])
    normals = []
    directions = []

    for vertice in tqdm(closest_vertices, desc="Processing vertices"):
        # Find the line between target point and closest point
        direction = vertice - voxel_coord
        normal = tangent_to_surface(verts, vertice, radius)
        
        # Store them to get the optimal ones in the end
        directions.append(direction)
        normals.append(normal)

        # Determine the similarity for each point in search area
        similarities = np.append(similarities, cosine_similarity(direction, normal))
    
    # Get the indices of the top nb_points similarities
    top_indices = np.argsort(np.abs(similarities))[-nb_points:][::-1]  # descending indices
    
    # Prepare results
    best_placements = []
    for idx in top_indices:
        intersect_point = closest_vertices[idx]
        # Calculate the distance
        intersect_point = point_translation(intersect_point, voxel_coord, normals[idx], gel_pad, file)
        distance = np.linalg.norm(intersect_point - voxel_coord)
        
        # Calculate the skull thickness at this point
        thickness = determine_thickness(binary_image, voxel_coord, intersect_point,  direction)

        if transducer_mode == 'CTX250':
            simu_intersect_pt = point_translation(intersect_point, voxel_coord, normals[idx], 3.5, file)
        elif transducer_mode == 'CTX500':
            simu_intersect_pt = point_translation(intersect_point, voxel_coord, normals[idx], 4.215, file)
        
        normal = normals[idx]
        direction = directions[idx]
        similarity = similarities[idx]
        
        tangent_plane = create_transducer_plot(normal, intersect_point, file)
        
        best_placements.append((tangent_plane, intersect_point, direction, similarity, simu_intersect_pt, distance, thickness))
    
    # Sort placements by thickness and then by distance
    best_placements = sorted(best_placements, key=lambda x: (x[6], x[5]))  # Sort by thickness, then by distance

    return best_placements
