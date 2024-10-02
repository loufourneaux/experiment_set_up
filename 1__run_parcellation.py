
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is using the transplanted lesion T1 and applys freesurfer's recon all function = all parts of the cortical reconstruction process
# Reconstructing a two dimensional cortical surface from a three dimensional volume (tranplanted T1 in this case)
# General freesurfer (fs): any 3D volumes are stored in "mri" folder
# First step: skull strip to create brainmask.mgz (specific fs extension)
# Second step: estimation of interface between white and grey matter for the two hemispheres > saved as lh.orig and rh.orig
# Refined and saved under lh.white and rh.white
# Edge of grey matter detected and saved under lh.pial and rh.pial (all if the above can be visualised with freeview)
# The pial files can be expanded (to better determine locations along the banks and ridges of the gyri) and are saved under lh.inflated and rh.inflated
# Inflated surfaces can be inflated again into a sphere, normalised to a template image (fsaverage, 40 subjects)
# Once individual surface map is normalised to the template, two atlases can be used for parcellation: Desikan-Killiany and Destrieux (more detailed)
# The -all flag instructs fs to run all processing steps 
# -openmp 12 specifies the number of OpenMP (multi-platform shared-memory multiprocessing programming) threads to be used -> here 12 to speed up the processing time
# -brainstem-structures tells fs to include segmentation of braistem structures in the processing pipeline

from __future__ import division

import argparse
import logging
import os
import json
import subprocess
import nibabel as nib
import numpy as np
import os.path
import shutil
import time
from copy import copy
from genericpath import isfile
from tools.utils import *
from tools.registration_ants import *



def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subj', nargs='+', dest='subj', help="Subject index (e.g., MOT05).")
    p.add_argument('--data_path', default='data', dest='data_path',
        help="Input folder path. ['%(default)s']")
    p.add_argument('--output_path', default='output', dest='output_path',
        help="Output folder path. ['%(default)s']")

    p.add_argument('-f', action='store_true', dest='isForce',
        help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p

def freesurfer_func(data_path: str, output_path: str, subj: str, isForce: bool):  
    """
    Run FreeSurfer processing for a given subject.

    This function performs several steps for FreeSurfer processing:
    1. Sets up the environment and paths.
    2. Copies T1-weighted MRI image to the appropriate FreeSurfer directory.
    3. Converts the T1 image to FreeSurfer's .mgz format.
    4. Runs the FreeSurfer `recon-all` command to process the MRI data.
    5. Segments the brainstem structures if not already done.

    Args:
        data_path (str): Path to the input data directory.
        output_path (str): Path to the FreeSurfer output directory.
        subj (str): Subject identifier.
        isForce (bool): If True, forces reprocessing even if outputs already exist.

    Returns:
        None
    """ 
    logging.info(f'Running freesurfer for subject: {subj}.')
    
    anat_folder = os.path.join(data_path, subj, "anat") 
    freesurfer_folder = os.path.join(output_path)
    
    os.environ["SUBJECTS_DIR"] = freesurfer_folder #where the recon-all stores result, otherwise will sotre in the freesurfer file

    # Check if folder exists, if not - create it, copy T1 image
    T1_image = os.path.join(anat_folder, f'{subj}_T1w.nii.gz') 
    T1_for_fs = os.path.join(data_path, 'MOTs', f"{subj}.nii") 

    #if os.path.isfile(T1_image):
    if not os.path.exists(os.path.join(freesurfer_folder, subj, "recon-all", "mri", "orig")):
        os.makedirs(os.path.join(freesurfer_folder, subj,"recon-all", "mri", "orig")) 

    if not os.path.isfile(T1_for_fs) or isForce:
        logging.info(f'Copying {os.path.basename(T1_image)} into fs folder.')
        shutil.copyfile(T1_image, T1_for_fs)

    # mri convert nii in mgz which is the type used in freesurfer 
    t1_mgz = os.path.join(freesurfer_folder, subj, "recon-all","mri", "orig", "001.mgz") 
    json_file = os.path.join(freesurfer_folder, subj, "recon-all","mri", "orig", "001.json")
    if os.path.isfile(t1_mgz) and not isForce:
        logging.info('T1 for freesurfer already converted')
    else:
        #convert from T1 to mgz
        mr_convert_cmd = f"mri_convert {T1_for_fs} {t1_mgz}"
        logging.info(f'mr convert command: "{mr_convert_cmd}".')
        subprocess.call(mr_convert_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': mr_convert_cmd,
                'Description': 'Conversion of T1 for freesurfer',
                'Anat_filename': t1_mgz,
                'Time': time.asctime()
            }
            json.dump(j, outfile)

    # recon all: the freesurfer command which amongst other will extract the brain from the skull
    # takes a long time to run (approx. 5 hours full server)
    fs_output = os.path.join(freesurfer_folder, subj, "recon-all", "mri", "brain.mgz")
    if os.path.isfile(fs_output) and not isForce:
        logging.info('Freesurfer already run')
    else:
        reconall_cmd = f"recon-all -all -subjid {subj} -openmp 12 -brainstem-structures"
        logging.info(f'recon all command: "{reconall_cmd}".')
        subprocess.call(reconall_cmd, shell=True)

    # segment BS, is not included in the recon all and parcellation otherwise
    fs_bs_output = os.path.join(freesurfer_folder, subj,"recon-all", "mri", "brainstemSsLabels.v13.mgz")
    if os.path.isfile(fs_bs_output) and not isForce:
        logging.info('Brainstem segmentation already done')
    else:
        segment_bs_cmd = f"segmentBS.sh {subj} {freesurfer_folder}"
        logging.info(f'segment BS command: "{segment_bs_cmd}".')
        subprocess.call(segment_bs_cmd, shell=True)


def glasserbuild_func(data_path:str, output_path:str,  subj:str, isForce:bool):
    """
    Run Glasser processing for a given subject.

    This function performs several steps for processing Glasser parcellation data:
    1. Sets up the environment and paths.
    2. Runs the `mri_surf2surf` command to map the annotation files for both hemispheres.
    3. Runs the `mri_aparc2aseg` command to map the annotations onto a volumetric image.
    4. Converts the brain.mgz file to NIfTI format.
    5. Transforms the volumetric image to T1 space.
    6. Converts the labels using `labelconvert`.

    Args:
        data_path (str): Path to the input data directory.
        output_path (str): Path to the FreeSurfer output directory.
        subj (str): Subject identifier.
        isForce (bool): If True, forces reprocessing even if outputs already exist.

    Returns:
        None
    """
    session_folder = os.path.join(output_path, subj, "1__output")
    logging.info('Processing Glasser subject: {0}.'.format(subj))


    freesurfer_folder = os.path.join(output_path, subj, "recon-all")
    os.environ["SUBJECTS_DIR"] = 'output'  #output of the freesurfer commands
    freesurfer_path='tools/fsaverage' 

    # Run the surf command for both hemispheres, map the annotation files so get a new annot file for next step
    for side in ("lh", "rh"):
        # target annotation
        annot_label_file = os.path.join(freesurfer_folder, "label", side + ".hcpmmp1.annot")
        if os.path.isfile(annot_label_file) and not isForce:
            logging.info('mri_surf2surf already done')
        else:
            # source annotation 
            annot_avg_file = os.path.join(freesurfer_path, "label", side + ".hcpmmp1.annot")
            mri_surf2surf_cmd = "mri_surf2surf --srcsubject fsaverage --trgsubject " + subj + " --hemi " + side + \
                " --sval-annot " + annot_avg_file + " --tval " + annot_label_file
            logging.info('mri_surf2surf command: "{0}".'.format(mri_surf2surf_cmd))
            subprocess.call(mri_surf2surf_cmd, shell=True)


    hcpmmp_file = os.path.join(session_folder, subj + "_hcpmmp1.mgz")
    if os.path.isfile(hcpmmp_file) and not isForce:
        logging.info('aparc2aseg already done')
    else:
        # map annotations onto volumetric image 
        mri_aparc2aseg_cmd = "mri_aparc2aseg --old-ribbon --s " + subj + " --annot hcpmmp1 --o " + hcpmmp_file
        logging.info('aparc2aseg command: "{0}".'.format(mri_aparc2aseg_cmd))
        subprocess.call(mri_aparc2aseg_cmd, shell=True)

    
    # Convert brain.mgz to NIfTI for the vol2vol 
    brain_mgz = os.path.join(freesurfer_folder, "mri", "brain.mgz")
    t1_brain_nii = os.path.join(session_folder, f"{subj}T1wBrain.nii.gz")
    if os.path.isfile(t1_brain_nii) and not isForce:
        logging.info(f'Brain-extracted T1-weighted image already exists: {t1_brain_nii}')
    else:
        brain_convert_cmd = f"mri_convert {brain_mgz} {t1_brain_nii}"
        logging.info(f'Brain extraction command: "{brain_convert_cmd}".')
        subprocess.call(brain_convert_cmd, shell=True)
    
    # Convert to NIfTI and transform to T1 space
    t1_brain_filename = os.path.join(session_folder, subj + "T1wBrain.nii.gz")
    hcpmmp_anat_file = os.path.join(session_folder, subj + "_hcpmmp1.nii.gz")
    if os.path.isfile(hcpmmp_anat_file) and not isForce:
        logging.info('aparc+aseg already in t1 space: "{0}".'.format(hcpmmp_anat_file))
        
    else:
        vol2vol_cmd = "mri_vol2vol --targ " + t1_brain_filename + " --mov " + hcpmmp_file + " --o " + hcpmmp_anat_file + " --regheader --interp nearest"
        logging.info('VOL2VOL command: "{0}".'.format(vol2vol_cmd))
        subprocess.call(vol2vol_cmd, shell=True)
    
    #this is already present in the server 
    mrtrix_labelconvert_folder = "/opt/mrtrix3/share/mrtrix3/labelconvert/"
    lookuptable_in = os.path.join(mrtrix_labelconvert_folder, "hcpmmp1_original.txt")
    lookuptable_out = os.path.join(mrtrix_labelconvert_folder, "hcpmmp1_ordered.txt")
    hcpmmp_anat_file_rein = os.path.join(session_folder, subj + "_hcpmmp1_rein.nii.gz")
    if os.path.isfile(hcpmmp_anat_file_rein) and not isForce:
        logging.info('labels already converted')
    else:
        label_convert_cmd = "labelconvert -force " + hcpmmp_anat_file + " " + lookuptable_in + " " + lookuptable_out + " " + hcpmmp_anat_file_rein
        logging.info('label convert command: "{0}".'.format(label_convert_cmd))
        subprocess.call(label_convert_cmd, shell=True)
    

if __name__ == "__main__":
    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce
    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)

    subj_list = [subj for subj in args.subj]

    data_path = args.data_path
    output_path=args.output_path

    if "all" in subj_list:
        subjects = [s for s in os.listdir(os.path.join(data_path, 'MOTs')) if os.path.isdir(os.path.join(data_path, s))]
    else:
        subjects = [subj for subj in subj_list]

    for subj in subjects:
        freesurfer_func(data_path, output_path, subj, isForce)
        glasserbuild_func(data_path, output_path, subj, isForce)
