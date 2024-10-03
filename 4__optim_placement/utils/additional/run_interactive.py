
#import subprocess
import os 
import nibabel as nib
import time
from scipy.ndimage import binary_fill_holes
from utils.utils import *


if __name__ =='__main__':
    subject=input("For which subject do you want to run the optimization? Enter number (eg. 05 for MOT05): ")
    subject=f'MOT{subject}'
    filename=f'../../../data/PCTs/{subject}_pCT.nii'
    filename2 = f'../../input/pCTs_thresholded/{subject}_pCT_threshold_300_2000.nii' 

    while os.path.isfile(filename)==False:
        subject =input("You need to enter a number (eg. 05) that exists in the folder pCTs. Enter number: ")
        filename=f'../../../PCTs/MOT{subject}_pCT.nii'
        filename2 = f'../../input/pCTs_thresholded/{subject}_pCT_threshold_300_2000.nii' 
        

    #LOAD DATA 
    start_time=time.time()
    # Load the thresholded MRI aka skull
    file=nib.load(filename)
    data= file.get_fdata()
    file_skull=nib.load(filename2)
    data_skull=file_skull.get_fdata()

    # Create binary image from 0 and other values image
    binary_data = (data > -700).astype(int)
    binary_data = binary_fill_holes(binary_data)
    binary_data_skull= (data_skull !=0).astype(int)
    binary_data_skull=binary_fill_holes(binary_data_skull)

    # Get the biggest component which is the skull using connected components (get rid of artifacts)
    binary_data=get_biggest_component(binary_data)
    binary_data_skull=get_biggest_component(binary_data_skull)

    # EXTRACT SURFACE MESH using the Marching Cubes algorithm
    verts, faces, _, _ = measure.marching_cubes(binary_data, level=0.5) #returns a set of vertices and triangular faces
    #FIND SURFACE POINTS
    surface_points = np.round(verts).astype(int)

    scatters=[]
    mesh3d = go.Mesh3d(
        x=verts[:, 0],
        y=verts[:, 1],
        z=verts[:, 2],
        i=faces[:, 0],
        j=faces[:, 1],
        k=faces[:, 2],
        opacity=0.5,
        color='gray',
    )
    scatters.append(mesh3d)


    # Prompt user for area choice
    area_choice = input("Enter area to target (eg. M1 / Thalamus): ")
    # Prompt user for voxel coordinates of area 
    target_coord = np.array([float(x.strip()) for x in input("Enter target coordinates in RAS (comma-separated x,y,z):").split(',')]+[1])
    # Transform the coordinates into voxel coordinates for this script
    voxel_coord = mri_to_voxel(file, target_coord) 

    #Mode of transducer for translation
    transducer_mode= input("Which transducer are you using (CTX500/CTX250)? Enter the model number: ").upper()
    while transducer_mode!='CTX250' and transducer_mode!= 'CTX500':
        transducer_mode=input("Please enter either 250 or 500 depending on your CTX model: ")
    #Mode of transducer for translation
    gel_pad= float(input("What is the thickness of the gel pad (in mm)? "))
    transducer_coord=[]
    
    #Target
    point_scatter = go.Scatter3d(
        x=[voxel_coord[0]],
        y=[voxel_coord[1]],
        z=[voxel_coord[2]],
        mode='markers',
        marker=dict(size=10, color='red'),
        name='Target'
    )
    scatters.append(point_scatter)

    #attention droite est gauche maybe dans header 
    #CHOSE MODE
    print("There are three modes of optimization of placement:\n [1] around the target \n [2] around an approximate transducer placement \n [3] search for every possibility of placement in the same side of the brain as the target")
    mode =5 #default to enter the loop
    while mode!= 1 and mode !=2 and mode != 3: 
        mode=int(input("Which placement optimization mode do you want to use (1/2/3)? "))
        #RADIUS BASED
        if mode ==1:
            # Define the radius of the sphere around target point
            #give min distance of sphere 
            search_area_points, sphere_scatter, message = create_sphere(file, surface_points, voxel_coord, interactive=True, radius=0)
            best_placements= find_optimal_placement(binary_data_skull, verts, voxel_coord, search_area_points, file, transducer_mode,gel_pad, nb_points=10, radius=5)
            scatters.append(sphere_scatter)
            #WHERE TO STORE RESULTS 
            dir= f"results/{subject}/target_based"
            ensure_directory_exists(dir)
            result_file = os.path.join(dir, f"{area_choice}_{target_coord[:3]}_TargetBasedPlacement.csv")
            visu_dir=os.path.join(dir, f"{area_choice}_{target_coord[:3]}_visualization.html")
            
        elif mode==2:
            # Prompt user for voxel coordinates
            print("Enter approximate transducer coordinates (comma-separated x,y,z):")
            transducer_coord = np.array([float(x.strip()) for x in input("Transducer coordinates: ").split(',')]+[1])
            # Transform the coordinates into voxel coordinates for this script
            transducer_voxel_coord = mri_to_voxel(file, transducer_coord)
            transducer_scatter= go.Scatter3d(
                x=[transducer_voxel_coord[0]],
                y=[transducer_voxel_coord[1]],
                z=[transducer_voxel_coord[2]],
                mode='markers',
                marker=dict(size=2, color='black'),
                name='Approximate transducer coordinates given'
            )
            scatters.append(transducer_scatter)

            #define the radius of the sphere around transducer coord 
            search_area_points, sphere_scatter, message= create_sphere(file, surface_points, transducer_voxel_coord, interactive=True, radius=0)
            scatters.append(sphere_scatter)
            best_placements= find_optimal_placement(binary_data_skull, verts, voxel_coord, search_area_points, file, transducer_mode,gel_pad, nb_points=10, radius=5)
            
            #WHERE TO STORE RESULTS
            dir= f"results/{subject}/transducer_based"
            ensure_directory_exists(dir)
            result_file = os.path.join(dir, f"{area_choice}_{target_coord[:3]}_TransducerBasedPlacement.csv")
            visu_dir=os.path.join(dir, f"{area_choice}_{target_coord[:3]}_visualization.html")
        
        elif mode==3:
            # The MRI is of shape (208,256,256) 
            # Remove wrong/useless surface points like the neck and face
            surface_points= surface_points[surface_points[:, 2] > 200]#approx. 110
            surface_points= surface_points[surface_points[:, 1] < 35]#approx. 160
            # Define the search area (hemisphere) split into left an right hemisphere along the sagittal plane
            if voxel_coord[0] < 104:
                search_area_points= surface_points[surface_points[:, 0] < 104]
            else:
                search_area_points= surface_points[surface_points[:, 0] >= 104]
            
            surface_scatter= go.Scatter3d(
                x=[point[0] for point in search_area_points],
                y=[point[1] for point in search_area_points],
                z=[point[2] for point in search_area_points],
                mode='markers',
                marker=dict(size=2, color='salmon', symbol='circle', opacity=0.7),
                name= 'Search area'
            )

            best_placements = find_optimal_placement(binary_data_skull, verts, voxel_coord, search_area_points, file, transducer_mode,gel_pad,  nb_points=100, radius=5)
            
            # Create a scatter plot for the intersect points
            # Split best_placements into the first 10 and the rest
            first_10 = best_placements[:10]
            remaining = best_placements[10:]

            annotations = [str(i+1) for i in range(len(first_10))]

            # Create a Scatter3d plot for the first 10 points
            scatter_first_10 = go.Scatter3d(
                x=[intersect_point[0] for (_, intersect_point, _, _, _,_,_) in first_10],
                y=[intersect_point[1] for (_, intersect_point, _, _, _,_,_) in first_10],
                z=[intersect_point[2] for (_, intersect_point, _, _, _,_,_) in first_10],
                mode='markers+text',
                marker=dict(size=4, color='blue', symbol='circle'), 
                text=annotations,
                textposition='top center',
                name='First 10 placements'
            )

            # Create a Scatter3d plot for the remaining points
            scatter_remaining = go.Scatter3d(
                x=[intersect_point[0] for (_, intersect_point, _, _, _,_,_) in remaining],
                y=[intersect_point[1] for (_, intersect_point, _, _, _,_,_) in remaining],
                z=[intersect_point[2] for (_, intersect_point, _, _, _,_,_) in remaining],
                mode='markers',
                marker=dict(size=4, color='orange', symbol='circle-open'),  
                name='Rest of the placements '
            )
            scatters.append(scatter_first_10)
            scatters.append(scatter_remaining)

            #WHERE TO SAVE RESULTS 
            dir=f"results/{subject}/all_placements"
            ensure_directory_exists(dir)
            result_file = os.path.join(dir, f"{area_choice}_{target_coord[:3]}_AllOptimizedPlacements.csv")
            visu_dir=os.path.join(dir, f"{area_choice}_{target_coord[:3]}_visualization.html")
            message=''
        else:
            print('Please enter 1, 2 or 3')

    surface_scatter= go.Scatter3d(
        x=[point[0] for point in search_area_points],
        y=[point[1] for point in search_area_points],
        z=[point[2] for point in search_area_points],
        mode='markers',
        marker=dict(size=2, color='salmon', symbol='circle', opacity=0.7),
        name= 'Search area'
    )
    scatters.append(surface_scatter)
    
    end_time=time.time()
    elapsed_time=end_time-start_time
    
    if mode==1 or mode==2:
        tangent_plane=best_placements[0][0]
        intersect_point=best_placements[0][1]
        direction=best_placements[0][2]
        similarity=best_placements[0][3]
        simu_intersect_point=best_placements[0][4]
        #Line scatter
        line = go.Scatter3d(
            x=[voxel_coord[0], voxel_coord[0] + direction[0] * 2],
            y=[voxel_coord[1], voxel_coord[1] + direction[1] * 2],
            z=[voxel_coord[2], voxel_coord[2] + direction[2] * 2],
            mode='lines',
            line=dict(color='black', width=1),
            name='US propagation direction'
        )

        #Intersection points
        intersect_scatter = go.Scatter3d(
            x=[intersect_point[0] ],
            y=[intersect_point[1] ],
            z=[intersect_point[2] ],
            mode='markers',
            marker=dict(size=5, color='blue'),
            name='Optimal placement of the transducer'
        )
        scatters.append(tangent_plane)
        scatters.append(line)
        scatters.append(intersect_scatter)

   
    write_results_to_file(file, result_file, target_coord, voxel_coord,search_area_points, area_choice, best_placements, elapsed_time, mode, transducer_coord, transducer_mode, message, gel_pad)
    
    #PLOT FINAL 
    plot_3D_figure(*scatters, directory=visu_dir)
