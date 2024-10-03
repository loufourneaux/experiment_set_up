
import os 
import nibabel as nib
import time
from scipy.ndimage import binary_fill_holes
from utils.utils import *
from utils.retrieve_coordinates import *


if __name__ =='__main__':
    
    # Read the Excel file
    df = pd.read_excel('input/subject_data.xlsx', engine='openpyxl')

    for index, row in df.iterrows():
        subject = row['Subject']
        filename = f'../data/PCTs/{subject}_pCT.nii'
        filename2 = f'input/pCTs_thresholded/{subject}_pCT_threshold_300_2000.nii' 

        if not os.path.isfile(filename):
            print(f"The file {filename} does not exist.")
            continue
        
        if not os.path.isfile(filename2):
            print(f"The file {filename2} does not exist.")
            continue
        
        start_time=time.time()

        # LOAD DATA 
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
        #level: specifies the threshold value at which to find the surface. 
        #Isovalue in the 3D volume data where the surface mesh will be extracted.
        #level=0.5 to find the boundary between 0 (background) and 1 (foreground)
        verts, faces, _, _ = measure.marching_cubes(binary_data, level=0.5) #returns a set of vertices and triangular faces
        # FIND SURFACE POINTS
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

        #get all of the information from the excel sheet
        area_choice = row['Area Choice']
        glasser = row['Use Glasser found target?']
        if glasser=='yes':
            target_coord= read_coordinates(f'../output/{subject}/2__output/coordinates_mri/cm_coordinates.csv', subject, area_choice)[0]
            target_coord.append(1) #for homogeneous 
        elif glasser =='no':
            target_coordinates_str = row['Target Coordinates (RAS)']
            target_coord = np.array([float(x.strip()) for x in target_coordinates_str.strip('()').split(',')] + [1])
        voxel_coord = mri_to_voxel(file, target_coord)

        transducer_mode = row['Transducer Model']
        transducer_coordinates_str = row['Transducer Coordinates (RAS)']
        transducer_coord = np.array([float(x.strip()) for x in transducer_coordinates_str.strip('()').split(',')] + [1])
        gel_pad = float(row['Gel Pad Thickness (mm)'])
        optimization_mode = int(row['Optimization Mode'])
        radius= float(row['Radius (mm)'])
        
        #Target plot
        point_scatter = go.Scatter3d(
            x=[voxel_coord[0]],
            y=[voxel_coord[1]],
            z=[voxel_coord[2]],
            mode='markers',
            marker=dict(size=5, color='black'),
            name='Target center of mass'
        )
        scatters.append(point_scatter)

        area_coords=read_coordinates(f'../output/{subject}/2__output/coordinates_voxel/all_coordinates.csv', subject, area_choice)
        area_scatter= go.Scatter3d(
            x=[voxel_coord[0] for voxel_coord in area_coords],
            y=[voxel_coord[1] for voxel_coord in area_coords],
            z=[voxel_coord[2] for voxel_coord in area_coords],
            mode='markers',
            marker=dict(size=2, color='red', opacity=0.1, symbol='circle-open'),
            name='Target volume'
        )
        scatters.append(area_scatter)


        # run differs depending the mode chosen
        #RADIUS BASED
        if optimization_mode ==1:
            # Define the radius of the sphere around target point
            #give min distance of sphere 
            search_area_points, sphere_scatter, message = create_sphere(file, surface_points, voxel_coord, False, radius)
            best_placements= find_optimal_placement(binary_data_skull, verts, voxel_coord, search_area_points, file, transducer_mode,gel_pad, nb_points=10, radius=5)
            scatters.append(sphere_scatter)
            #WHERE TO STORE RESULTS 
            dir= f"results/{subject}/target_based"
            ensure_directory_exists(dir)
            result_file = os.path.join(dir, f"{area_choice}_{target_coord[:3]}_TargetBasedPlacement.csv")
            visu_dir=os.path.join(dir, f"{area_choice}_{target_coord[:3]}_visualization.html")
            
        elif optimization_mode==2:
            # Transform the coordinates into voxel coordinates for this script
            transducer_voxel_coord = mri_to_voxel(file, transducer_coord)
            #plot the estimated transducer placement
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
            search_area_points, sphere_scatter, message= create_sphere(file, surface_points, transducer_voxel_coord,False, radius)
            scatters.append(sphere_scatter)
            best_placements= find_optimal_placement(binary_data_skull, verts, voxel_coord, search_area_points, file, transducer_mode,gel_pad, nb_points=10, radius=5)
            
            #WHERE TO STORE RESULTS
            dir= f"results/{subject}/transducer_based"
            ensure_directory_exists(dir)
            result_file = os.path.join(dir, f"{area_choice}_{target_coord[:3]}_TransducerBasedPlacement.csv")
            visu_dir=os.path.join(dir, f"{area_choice}_{target_coord[:3]}_visualization.html")

        elif optimization_mode==3:
            # The MRI is of shape (208,256,256) 
            # Remove wrong/useless surface points like the neck and face
            surface_points= surface_points[surface_points[:, 2] > 110]#approx. 110, for test 200
            surface_points= surface_points[surface_points[:, 1] < 160]#approx. 160, for test 25
            # Define the search area (hemisphere) split into left an right hemisphere along the sagittal plane
            if voxel_coord[0] < 104:
                search_area_points= surface_points[surface_points[:, 0] < 104]
            else:
                search_area_points= surface_points[surface_points[:, 0] >= 104]
            
            best_placements = find_optimal_placement(binary_data_skull, verts, voxel_coord, search_area_points, file, transducer_mode,gel_pad,  nb_points=1000, radius=5)
            
            # Create a scatter plot for the intersect points
            # Split best_placements into the first 10 and the rest to see if there are some clusters
            first_10 = best_placements[:100]
            remaining = best_placements[100:]

            annotations = [str(i+1) for i in range(len(first_10))]

            # Create a Scatter3d plot for the first 10 points
            scatter_first_10 = go.Scatter3d(
                x=[intersect_point[0] for (_, intersect_point, _, _, _,_, _) in first_10],
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
                x=[intersect_point[0] for (_, intersect_point, _, _, _, _,_) in remaining],
                y=[intersect_point[1] for (_, intersect_point, _, _, _,_,_) in remaining],
                z=[intersect_point[2] for (_, intersect_point, _, _, _,_,_) in remaining],
                mode='markers',
                marker=dict(size=4, color='orange', symbol='circle-open'),  
                name='Rest of the placements'
            )
            scatters.append(scatter_first_10)
            scatters.append(scatter_remaining)

            #WHERE TO SAVE RESULTS OF MODE 3
            dir=f"results/{subject}/all_placements"
            ensure_directory_exists(dir)
            result_file = os.path.join(dir, f"{area_choice}_{target_coord[:3]}_AllOptimizedPlacements.csv")
            visu_dir=os.path.join(dir, f"{area_choice}_{target_coord[:3]}_visualization.html")
            message=''

        #plot the search area surface 
        surface_scatter= go.Scatter3d(
            x=[point[0] for point in search_area_points],
            y=[point[1] for point in search_area_points],
            z=[point[2] for point in search_area_points],
            mode='markers',
            marker=dict(size=2, color='salmon', symbol='circle', opacity=0.7),
            name= 'Search area'
        )
        scatters.append(surface_scatter)
        
        #timed the script 
        end_time=time.time()
        elapsed_time=end_time-start_time
        
        #if the mode is 1 or 2, the plot will be a bit different 
        if optimization_mode==1 or optimization_mode==2:
            #retrieve the first optimal point, its tangent, the US propagation direction...
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

        #save all of the information to a csv file
        write_results_to_file(file, result_file, target_coord, voxel_coord,search_area_points, area_choice, best_placements, elapsed_time, optimization_mode, transducer_coord, transducer_mode, message, gel_pad)
        
        #PLOT FINAL 
        plot_3D_figure(*scatters, directory=visu_dir)
