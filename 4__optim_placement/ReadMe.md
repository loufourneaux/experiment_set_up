# Optimizer of transducer placement

Optimization of placement of the transducer for FUS based on the beam direction (1), the skull thickness (2), the distance of the transducer to the target area (3). 

## Description

The transducer placement needs to be optimized in order to reach the target. 
There are three different ways you can otpimize the placement:
    [1]: chose the radius of a sphere around the target area in which you want to search for the optimal points--> the points will also be the closest possible
    [2]: chose the radius of a sphere around your chosen transducer placement in whch you want to search for the optimal placement --> optimize an already given placement
    [3]: compute all the best placements situated in the same hemisphere of the scalp as the target point.

## Dependencies
Run the following command: 

```
pip install -r requirements.txt
```

### Executing program
# Run only one optimization

If you want to run the optimization for only one subject:
```
python run_interactive.py
```

# Run multiple optimizations

Fill in the excel with all the options you want for each subject and it will run the optimization adapted to every subject's parameters. 


## Architecture of folder example


├── input
│   ├── subject_data
│   └── pCTs_thresholded
│        ├──MOT05_thresholded_300.2000....
├── results
│    ├──MOTXX   
│       ├──all_placements
│       ├──transducer_based
│       └──target_based
├── utils
│    ├── utils.py
│    ├── retrieve_coordinates.py
│    ├── additional 
│    │    ├──run_interactive.py
│    │    └──all_viusalization.py
│    └── MRI_to_voxels.py
├── 
├── run_all.py
└── README.md


## Authors

Semester Project
Lou Fourneaux, lou.fourneaux@epfl.ch

## example points in RAS

R_thalamus: (7.372,5.527,1.476)
1rst placement: (39.231,14.931,61.684)
2nd placement: (57.554,12.169,47.391)
3rd placement: (5.370,-44.128,72.493)
4th placement: (5.184,-71.087,53.836)