

#---------------------------------
# New invocation of recon-all Wed Jul 24 12:24:05 CEST 2024 
#--------------------------------------------
#@# MotionCor Wed Jul 24 12:24:06 CEST 2024

 cp /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/orig/001.mgz /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/rawavg.mgz 


 mri_convert /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/rawavg.mgz /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/transforms/talairach.xfm /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/orig.mgz /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/orig.mgz 

#--------------------------------------------
#@# Talairach Wed Jul 24 12:24:19 CEST 2024

 mri_nu_correct.mni --no-rescale --i orig.mgz --o orig_nu.mgz --n 1 --proto-iters 1000 --distance 50 


 talairach_avi --i orig_nu.mgz --xfm transforms/talairach.auto.xfm 

talairach_avi log file is transforms/talairach_avi.log...

 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Wed Jul 24 12:26:20 CEST 2024

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /opt/freesurfer/bin/extract_talairach_avi_QA.awk /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/transforms/talairach_avi.log 


 tal_QC_AZS /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Nu Intensity Correction Wed Jul 24 12:26:20 CEST 2024

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 


 mri_add_xform_to_header -c /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/transforms/talairach.xfm nu.mgz nu.mgz 

#--------------------------------------------
#@# Intensity Normalization Wed Jul 24 12:28:51 CEST 2024

 mri_normalize -g 1 -mprage nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Wed Jul 24 12:30:58 CEST 2024

 mri_em_register -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mri_em_register.skull.dat -skull nu.mgz /opt/freesurfer/average/RB_all_withskull_2016-05-10.vc700.gca transforms/talairach_with_skull.lta 


 mri_watershed -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mri_watershed.dat -T1 -brain_atlas /opt/freesurfer/average/RB_all_withskull_2016-05-10.vc700.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Wed Jul 24 12:47:01 CEST 2024

 mri_em_register -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mri_em_register.dat -uns 3 -mask brainmask.mgz nu.mgz /opt/freesurfer/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Wed Jul 24 12:54:15 CEST 2024

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /opt/freesurfer/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Wed Jul 24 12:56:10 CEST 2024

 mri_ca_register -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mri_ca_register.dat -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /opt/freesurfer/average/RB_all_2016-05-10.vc700.gca transforms/talairach.m3z 

#--------------------------------------
#@# SubCort Seg Wed Jul 24 15:40:25 CEST 2024

 mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz transforms/talairach.m3z /opt/freesurfer/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/mri/transforms/cc_up.lta MOT01 

#--------------------------------------
#@# Merge ASeg Wed Jul 24 16:20:22 CEST 2024

 cp aseg.auto.mgz aseg.presurf.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Wed Jul 24 16:20:22 CEST 2024

 mri_normalize -mprage -aseg aseg.presurf.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Wed Jul 24 16:22:30 CEST 2024

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Wed Jul 24 16:22:31 CEST 2024

 mri_segment -mprage brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.presurf.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Wed Jul 24 16:24:08 CEST 2024

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Wed Jul 24 16:24:36 CEST 2024

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Tessellate rh Wed Jul 24 16:24:40 CEST 2024

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Wed Jul 24 16:24:44 CEST 2024

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Smooth1 rh Wed Jul 24 16:24:49 CEST 2024

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Wed Jul 24 16:24:54 CEST 2024

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# Inflation1 rh Wed Jul 24 16:25:09 CEST 2024

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Wed Jul 24 16:25:25 CEST 2024

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# QSphere rh Wed Jul 24 16:27:09 CEST 2024

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology Copy lh Wed Jul 24 16:28:48 CEST 2024

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 

#--------------------------------------------
#@# Fix Topology Copy rh Wed Jul 24 16:28:48 CEST 2024

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 

#@# Fix Topology lh Wed Jul 24 16:28:48 CEST 2024

 mris_fix_topology -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_fix_topology.lh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 MOT01 lh 

#@# Fix Topology rh Wed Jul 24 16:39:15 CEST 2024

 mris_fix_topology -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_fix_topology.rh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 MOT01 rh 


 mris_euler_number ../surf/lh.orig 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make White Surf lh Wed Jul 24 16:54:39 CEST 2024

 mris_make_surfaces -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs MOT01 lh 

#--------------------------------------------
#@# Make White Surf rh Wed Jul 24 16:58:26 CEST 2024

 mris_make_surfaces -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs MOT01 rh 

#--------------------------------------------
#@# Smooth2 lh Wed Jul 24 17:02:12 CEST 2024

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white.preaparc ../surf/lh.smoothwm 

#--------------------------------------------
#@# Smooth2 rh Wed Jul 24 17:02:16 CEST 2024

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white.preaparc ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Wed Jul 24 17:02:21 CEST 2024

 mris_inflate -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_inflate.lh.dat ../surf/lh.smoothwm ../surf/lh.inflated 

#--------------------------------------------
#@# Inflation2 rh Wed Jul 24 17:02:35 CEST 2024

 mris_inflate -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_inflate.rh.dat ../surf/rh.smoothwm ../surf/rh.inflated 

#--------------------------------------------
#@# Curv .H and .K lh Wed Jul 24 17:02:51 CEST 2024

 mris_curvature -w lh.white.preaparc 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 lh.inflated 

#--------------------------------------------
#@# Curv .H and .K rh Wed Jul 24 17:03:48 CEST 2024

 mris_curvature -w rh.white.preaparc 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 rh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Wed Jul 24 17:04:45 CEST 2024

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm MOT01 lh curv sulc 


#-----------------------------------------
#@# Curvature Stats rh Wed Jul 24 17:04:49 CEST 2024

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm MOT01 rh curv sulc 

#--------------------------------------------
#@# Sphere lh Wed Jul 24 17:04:52 CEST 2024

 mris_sphere -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_sphere.lh.dat -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Sphere rh Wed Jul 24 17:14:52 CEST 2024

 mris_sphere -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_sphere.rh.dat -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg lh Wed Jul 24 17:23:27 CEST 2024

 mris_register -curv -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_register.lh.dat ../surf/lh.sphere /opt/freesurfer/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Surf Reg rh Wed Jul 24 17:40:41 CEST 2024

 mris_register -curv -rusage /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/touch/rusage.mris_register.rh.dat ../surf/rh.sphere /opt/freesurfer/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Wed Jul 24 18:01:48 CEST 2024

 mris_jacobian ../surf/lh.white.preaparc ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# Jacobian white rh Wed Jul 24 18:01:49 CEST 2024

 mris_jacobian ../surf/rh.white.preaparc ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Wed Jul 24 18:01:51 CEST 2024

 mrisp_paint -a 5 /opt/freesurfer/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#--------------------------------------------
#@# AvgCurv rh Wed Jul 24 18:01:52 CEST 2024

 mrisp_paint -a 5 /opt/freesurfer/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Wed Jul 24 18:01:53 CEST 2024

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MOT01 lh ../surf/lh.sphere.reg /opt/freesurfer/average/lh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.annot 

#-----------------------------------------
#@# Cortical Parc rh Wed Jul 24 18:02:04 CEST 2024

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MOT01 rh ../surf/rh.sphere.reg /opt/freesurfer/average/rh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf lh Wed Jul 24 18:02:14 CEST 2024

 mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs MOT01 lh 

#--------------------------------------------
#@# Make Pial Surf rh Wed Jul 24 18:12:52 CEST 2024

 mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs MOT01 rh 

#--------------------------------------------
#@# Surf Volume lh Wed Jul 24 18:24:11 CEST 2024
#--------------------------------------------
#@# Surf Volume rh Wed Jul 24 18:24:14 CEST 2024
#--------------------------------------------
#@# Cortical ribbon mask Wed Jul 24 18:24:17 CEST 2024

 mris_volmask --aseg_name aseg.presurf --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon MOT01 

#-----------------------------------------
#@# Parcellation Stats lh Wed Jul 24 18:30:43 CEST 2024

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab MOT01 lh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.pial.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab MOT01 lh pial 

#-----------------------------------------
#@# Parcellation Stats rh Wed Jul 24 18:31:30 CEST 2024

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab MOT01 rh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.pial.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab MOT01 rh pial 

#-----------------------------------------
#@# Cortical Parc 2 lh Wed Jul 24 18:32:19 CEST 2024

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MOT01 lh ../surf/lh.sphere.reg /opt/freesurfer/average/lh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Cortical Parc 2 rh Wed Jul 24 18:32:32 CEST 2024

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MOT01 rh ../surf/rh.sphere.reg /opt/freesurfer/average/rh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Wed Jul 24 18:32:46 CEST 2024

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab MOT01 lh white 

#-----------------------------------------
#@# Parcellation Stats 2 rh Wed Jul 24 18:33:10 CEST 2024

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab MOT01 rh white 

#-----------------------------------------
#@# Cortical Parc 3 lh Wed Jul 24 18:33:34 CEST 2024

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MOT01 lh ../surf/lh.sphere.reg /opt/freesurfer/average/lh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Cortical Parc 3 rh Wed Jul 24 18:33:45 CEST 2024

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MOT01 rh ../surf/rh.sphere.reg /opt/freesurfer/average/rh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Parcellation Stats 3 lh Wed Jul 24 18:33:56 CEST 2024

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.DKTatlas.stats -b -a ../label/lh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab MOT01 lh white 

#-----------------------------------------
#@# Parcellation Stats 3 rh Wed Jul 24 18:34:20 CEST 2024

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.DKTatlas.stats -b -a ../label/rh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab MOT01 rh white 

#-----------------------------------------
#@# WM/GM Contrast lh Wed Jul 24 18:34:44 CEST 2024

 pctsurfcon --s MOT01 --lh-only 

#-----------------------------------------
#@# WM/GM Contrast rh Wed Jul 24 18:34:49 CEST 2024

 pctsurfcon --s MOT01 --rh-only 

#-----------------------------------------
#@# Relabel Hypointensities Wed Jul 24 18:34:54 CEST 2024

 mri_relabel_hypointensities aseg.presurf.mgz ../surf aseg.presurf.hypos.mgz 

#-----------------------------------------
#@# AParc-to-ASeg aparc Wed Jul 24 18:35:10 CEST 2024

 mri_aparc2aseg --s MOT01 --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /opt/freesurfer/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt 

#-----------------------------------------
#@# AParc-to-ASeg a2009s Wed Jul 24 18:38:39 CEST 2024

 mri_aparc2aseg --s MOT01 --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /opt/freesurfer/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt --a2009s 

#-----------------------------------------
#@# AParc-to-ASeg DKTatlas Wed Jul 24 18:42:04 CEST 2024

 mri_aparc2aseg --s MOT01 --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /opt/freesurfer/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt --annot aparc.DKTatlas --o mri/aparc.DKTatlas+aseg.mgz 

#-----------------------------------------
#@# APas-to-ASeg Wed Jul 24 18:45:30 CEST 2024

 apas2aseg --i aparc+aseg.mgz --o aseg.mgz 

#--------------------------------------------
#@# ASeg Stats Wed Jul 24 18:45:35 CEST 2024

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --brainmask mri/brainmask.mgz --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --euler --ctab /opt/freesurfer/ASegStatsLUT.txt --subject MOT01 

#-----------------------------------------
#@# WMParc Wed Jul 24 18:46:02 CEST 2024

 mri_aparc2aseg --s MOT01 --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject MOT01 --surf-wm-vol --ctab /opt/freesurfer/WMParcStatsLUT.txt --etiv 

#--------------------------------------------
#@# BA_exvivo Labels lh Wed Jul 24 18:51:08 CEST 2024

 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA1_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA1_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA2_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA2_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA3a_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA3a_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA3b_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA3b_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA4a_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA4a_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA4p_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA4p_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA6_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA6_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA44_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA44_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA45_exvivo.label --trgsubject MOT01 --trglabel ./lh.BA45_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.V1_exvivo.label --trgsubject MOT01 --trglabel ./lh.V1_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.V2_exvivo.label --trgsubject MOT01 --trglabel ./lh.V2_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.MT_exvivo.label --trgsubject MOT01 --trglabel ./lh.MT_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.entorhinal_exvivo.label --trgsubject MOT01 --trglabel ./lh.entorhinal_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.perirhinal_exvivo.label --trgsubject MOT01 --trglabel ./lh.perirhinal_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA1_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA1_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA2_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA2_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA3a_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA3a_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA3b_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA3b_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA4a_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA4a_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA4p_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA4p_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA6_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA6_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA44_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA44_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.BA45_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.BA45_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.V1_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.V1_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.V2_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.V2_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.MT_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.MT_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.entorhinal_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.entorhinal_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/lh.perirhinal_exvivo.thresh.label --trgsubject MOT01 --trglabel ./lh.perirhinal_exvivo.thresh.label --hemi lh --regmethod surface 


 mris_label2annot --s MOT01 --hemi lh --ctab /opt/freesurfer/average/colortable_BA.txt --l lh.BA1_exvivo.label --l lh.BA2_exvivo.label --l lh.BA3a_exvivo.label --l lh.BA3b_exvivo.label --l lh.BA4a_exvivo.label --l lh.BA4p_exvivo.label --l lh.BA6_exvivo.label --l lh.BA44_exvivo.label --l lh.BA45_exvivo.label --l lh.V1_exvivo.label --l lh.V2_exvivo.label --l lh.MT_exvivo.label --l lh.entorhinal_exvivo.label --l lh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose 


 mris_label2annot --s MOT01 --hemi lh --ctab /opt/freesurfer/average/colortable_BA.txt --l lh.BA1_exvivo.thresh.label --l lh.BA2_exvivo.thresh.label --l lh.BA3a_exvivo.thresh.label --l lh.BA3b_exvivo.thresh.label --l lh.BA4a_exvivo.thresh.label --l lh.BA4p_exvivo.thresh.label --l lh.BA6_exvivo.thresh.label --l lh.BA44_exvivo.thresh.label --l lh.BA45_exvivo.thresh.label --l lh.V1_exvivo.thresh.label --l lh.V2_exvivo.thresh.label --l lh.MT_exvivo.thresh.label --l lh.entorhinal_exvivo.thresh.label --l lh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose 


 mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.stats -b -a ./lh.BA_exvivo.annot -c ./BA_exvivo.ctab MOT01 lh white 


 mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.thresh.stats -b -a ./lh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab MOT01 lh white 

#--------------------------------------------
#@# BA_exvivo Labels rh Wed Jul 24 18:54:57 CEST 2024

 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA1_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA1_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA2_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA2_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA3a_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA3a_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA3b_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA3b_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA4a_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA4a_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA4p_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA4p_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA6_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA6_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA44_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA44_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA45_exvivo.label --trgsubject MOT01 --trglabel ./rh.BA45_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.V1_exvivo.label --trgsubject MOT01 --trglabel ./rh.V1_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.V2_exvivo.label --trgsubject MOT01 --trglabel ./rh.V2_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.MT_exvivo.label --trgsubject MOT01 --trglabel ./rh.MT_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.entorhinal_exvivo.label --trgsubject MOT01 --trglabel ./rh.entorhinal_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.perirhinal_exvivo.label --trgsubject MOT01 --trglabel ./rh.perirhinal_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA1_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA1_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA2_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA2_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA3a_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA3a_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA3b_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA3b_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA4a_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA4a_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA4p_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA4p_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA6_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA6_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA44_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA44_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.BA45_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.BA45_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.V1_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.V1_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.V2_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.V2_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.MT_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.MT_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.entorhinal_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.entorhinal_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /data/MOTUS/Folder_lou/parcellation_V2/output/fsaverage/label/rh.perirhinal_exvivo.thresh.label --trgsubject MOT01 --trglabel ./rh.perirhinal_exvivo.thresh.label --hemi rh --regmethod surface 


 mris_label2annot --s MOT01 --hemi rh --ctab /opt/freesurfer/average/colortable_BA.txt --l rh.BA1_exvivo.label --l rh.BA2_exvivo.label --l rh.BA3a_exvivo.label --l rh.BA3b_exvivo.label --l rh.BA4a_exvivo.label --l rh.BA4p_exvivo.label --l rh.BA6_exvivo.label --l rh.BA44_exvivo.label --l rh.BA45_exvivo.label --l rh.V1_exvivo.label --l rh.V2_exvivo.label --l rh.MT_exvivo.label --l rh.entorhinal_exvivo.label --l rh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose 


 mris_label2annot --s MOT01 --hemi rh --ctab /opt/freesurfer/average/colortable_BA.txt --l rh.BA1_exvivo.thresh.label --l rh.BA2_exvivo.thresh.label --l rh.BA3a_exvivo.thresh.label --l rh.BA3b_exvivo.thresh.label --l rh.BA4a_exvivo.thresh.label --l rh.BA4p_exvivo.thresh.label --l rh.BA6_exvivo.thresh.label --l rh.BA44_exvivo.thresh.label --l rh.BA45_exvivo.thresh.label --l rh.V1_exvivo.thresh.label --l rh.V2_exvivo.thresh.label --l rh.MT_exvivo.thresh.label --l rh.entorhinal_exvivo.thresh.label --l rh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose 


 mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.stats -b -a ./rh.BA_exvivo.annot -c ./BA_exvivo.ctab MOT01 rh white 


 mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.thresh.stats -b -a ./rh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab MOT01 rh white 

#--------------------------------------------
#@# Segmentation of brainstem substructures  Wed Jul 24 18:58:48 CEST 2024

 /opt/freesurfer/bin/segmentBS.sh /opt/freesurfer/MCRv80 /opt/freesurfer MOT01 /data/MOTUS/Folder_lou/parcellation_V2/output 

See log file: /data/MOTUS/Folder_lou/parcellation_V2/output/MOT01/scripts/brainstem-structures.log
