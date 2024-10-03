
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.ndimage import binary_fill_holes
from utils.utils import *
import seaborn as sns


#insert the file path for the points you want to plot and corresponding pCT
file_path= '/Users/loufourneaux/Desktop/EPFL/projet_II/optim_placement/results/MOT01/all_placements/R_thalamus_[7.372 5.527 1.476]_AllOptimizedPlacements.csv'
filename=f'pCTs/MOT01_pCT.nii'

file=nib.load(filename)
data= file.get_fdata()

binary_data = (data > -700).astype(int)
binary_data = binary_fill_holes(binary_data)
binary_data=get_biggest_component(binary_data)
# EXTRACT SURFACE MESH using the Marching Cubes algorithm
verts, faces, _, _ = measure.marching_cubes(binary_data, level=0.5)

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

df = pd.read_csv(file_path, index_col=0)


# Convert voxel coordinates from strings to lists of rounded integers directly
def parse_and_round_voxel_coordinates(voxel_str):
    try:
        # Remove square brackets and split by spaces
        voxel_coords = voxel_str.strip('[]').split()
        # Convert to floats and then round to integers
        return [int(round(float(coord))) for coord in voxel_coords]
    except:
        return np.nan

df['Exit Plane Voxel Coordinates'] = df['Exit Plane Voxel Coordinates'].apply(parse_and_round_voxel_coordinates)
df = df.dropna(subset=['Exit Plane Voxel Coordinates'])


# Sort the DataFrame by 'Distances' or 'Thicknesses'
sort_by = 'Thicknesses'  # Change to 'Thicknesses' if needed
sorted_df = df.sort_values(by=sort_by)

# Extract sorted voxel coordinates
sorted_voxels = sorted_df['Exit Plane Voxel Coordinates'].values

# Define color scale based on sorted 'Distances' or 'Thicknesses'
distances = sorted_df[f'{sort_by}'].values  # or 'Thicknesses'
colorscale = plt.cm.inferno(distances)

# Calculate the average thickness and its standard deviation
average_thickness = np.mean(distances)
std_thickness = np.std(distances)

# Create a bar plot of the distances
plt.figure(figsize=(10, 6))
sns.barplot(x=np.arange(len(distances)), y=distances, palette='inferno')
plt.xlabel('Index')
plt.ylabel('Thickness (mm)')
plt.title('Bar Plot of Thickness')
# Hide x-axis ticks
plt.xticks([])

# Overlay the average thickness and standard deviation
plt.axhline(y=average_thickness, color='blue', linestyle='--', label=f'Avg Thickness: {average_thickness:.2f}')
plt.fill_between(np.arange(len(distances)), average_thickness - std_thickness, average_thickness + std_thickness, color='blue', alpha=0.2, label=f'Std Dev: {std_thickness:.2f}')

# Add a legend
plt.legend()

# Show the plot
plt.show()

# Create Scatter3d plot
point_scatter = go.Scatter3d(
    x=[voxel[0] for voxel in sorted_voxels],
    y=[voxel[1] for voxel in sorted_voxels],
    z=[voxel[2] for voxel in sorted_voxels],
    mode='markers',
        marker=dict(
        size=2,
        color=distances,
        colorscale='Inferno',
        colorbar=dict(title=sort_by),
        opacity=0.8
    ),
    name='Points'
)

#plot_3D_figure( mesh3d, point_scatter, directory=f'points_sorted_by_{sort_by}.html')