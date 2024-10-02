# Parcellation

Parcellation of T1w with Glasser 360 Atlas and extract the coordinates from the parcellated brain files. 

### Dependencies

# Softwares already on server!

To install the packages: 
```
pip install -r requirements.py
```

# Information on the results data 

* The hcpmmp1_ordered.txt file contains the label id, label name, and color+opacity of each area found in the parcellated MRI
* Parcellated data contains a value between 1 and 379 refering to an area in each of the voxels in the 3D space (256, 256, 256). This is the dimensions in which it is stored by freesurfer after skull-stripping.
* The T1w is in the shape (208, 256, 256) because it was acquired by the MRI machine and these are the usual dimensions.

### Executing program

The program needs to be executed in the following way to retrieve all of the information. 

# Execute the run_parcellation.py file using the following command 

```
python 1__run_parcellation.py --subj MOT01 --data_path data --output_path output -f  -v
```

* --subj : list of the MRIs want to parcellate
* --data_path : the folder where the MRIs are 
* --output_path : the folder where the recon-all result will be stored and the fsaverage is stored 
* -f: only to overwrite the files already found 
* -v: verbose is true 

# Excecute the run_find_areas.py file to find the list of the coordinates af the areas. 


```
python 2__run_find_areas.py 
```

Results are stored in a dictionary.


# Execute the run_plot_one.py to get the coordinates of left and right M1 and thalamus of 1 subject of choice and a plot of these areas in 3D

Before running, copy all the corresponding pseudoCTs file needed for the plotting the names are of the form MOTXX_pCT.nii. 

Takes the results from the 2__run_find_areas so need to run it before. 

```
python 3__run_visualization.py --subj 05 --areas m1 --surroundings S1,SMAa,SMAp,SCEF,anterior6,ventral6,dorsal6,rostral6
```
* --subj: put the subject you want to visualize
* -- areas: m1 or thalamus or both 
* -- surroundings: all of the surrounding areas of interest you want to see (they are listed in the color-dict in the script)


## Architecture of folder example


├── output
│   ├── MOT05
│   │   ├── recon-all
│   │   │   ├── mri
│   │   │   └── ...
│   │   ├── 1__output
│   │   ├── 2__output
│   │   └── 3__output
│   └── MOT...
│       ├── ...
│       └── ...
├── data
│    ├──MOTs
│    │    ├──MOT05
│    │    └──...
│    └── PCTs
│        └── MOT05_pCTs
├── tools  
│    ├── utils.py
│    ├──MRI_to_voxels.py
│    ├── retrieve_coordinates.py
│    ├── fsaverage
│    │   ├── label
│    │   ├── mri
│    │   ├──mri.2mm
│    │   ├──scripts
│    │   ├──surf
│    │   └──xhemi
│    └──registration_ants.py
├── 1__run_parcellation.py
├── 2__run_find_areas.py
├── 3__run_visualization.py
├── 4__optim_placement.py
│    ├──.... see ReadMe of folder
├── utils.py
└── README.md


## Authors

Contributors names and contact info: 

original script by Fabienne Windel, UPHummel lab 
modified by Lou Fourneaux, lou.fourneaux@epfl.ch