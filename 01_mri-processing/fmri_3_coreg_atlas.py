########################################################
########################################################
# Author: Lionel Butry (@lbutry)
# Purpose: 
#  (1) Performs atlas individualization with the following steps:
#       - Coregister atlas to T1w space
#       - Perform subcortical segmentation (fsl.FLIRT)
#       - Combine subcortical & GM mask
#       - Mask atlas with combined mask
# Requirements:
#   - fmriprep output
#   - Installation of 'FSL', 'ANTs'
########################################################
########################################################

# Import modules
import os
import shutil
from os.path import join
from nipype import Node, Workflow
from nipype.interfaces import fsl, ants
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.utility import IdentityInterface
import multiprocessing as mp

###################### USER INPUT ######################
#id_list = id_list # define subjects to be processed: ['PREDICT12345', ...]; DO NOT define, if you use the run_preproc_wrapper # type: ignore
# Paths
folder_fmriprep = "f01_fmriprep"
wf_name = "f03_coreg_atlas"
path_derivatives = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives"
path_atlas = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/utils/BN_Atlas_246_1mm.nii"
# Computational parameters
n_procs = os.cpu_count() - 1 # Number of CPUs
memory_gb = 15 #mac
n_gpu = 0 #mac
ngpuproc = 0 #mac
########################################################

# Get IDs which are already processed
already_processed = os.listdir(f"{path_derivatives}/{wf_name}")
already_processed = [folder for folder in already_processed if "sub-PREDICT" in folder] # select only subject output folders

#Get all IDs from fmriprep output
fmriprep_id_list = os.listdir(f"{path_derivatives}/{folder_fmriprep}") # list of all folders
fmriprep_id_list = [folder for folder in fmriprep_id_list if "PREDICT" in folder and not folder.endswith('.html')] # select only subject output folders

# IDs to be processed
id_list = set(fmriprep_id_list) - set(already_processed)
id_list = [f.strip('sub-') for f in id_list] # strip the str to obtain ID

# Define input & iterate over all subs
idsource = Node(IdentityInterface(fields=['id']), name="infosource")
idsource.iterables = [('id', id_list)]

path_t1w        = join(folder_fmriprep, "sub-{id}", "anat", "sub-{id}_acq-anat_desc-preproc_T1w.nii.gz")
path_transforms = join(folder_fmriprep, "sub-{id}", "anat", "sub-{id}_acq-anat_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5")
path_dseg       = join(folder_fmriprep, "sub-{id}", "anat", "sub-{id}_acq-anat_dseg.nii.gz")

path_dict = {"atlas": path_atlas,
             "t1w": path_t1w,
             "transform": path_transforms,
             "dseg": path_dseg}

selectfiles = Node(SelectFiles(templates=path_dict, base_directory=path_derivatives), name="selectfiles")

# Create Workflow
wf = Workflow(name=wf_name, base_dir=path_derivatives)
wf.connect([(idsource, selectfiles, [("id", "id")])])

## (1) Coregister Atlas -> T1w
applytransform = Node(ants.ApplyTransforms(), name = "applytransform")
applytransform.inputs.interpolation = "NearestNeighbor"
applytransform.inputs.output_image = "_atlas_coreg.nii.gz"

wf.connect([
    (selectfiles, applytransform, [("t1w", "reference_image"),
                                   ("atlas", "input_image"),
                                   ("transform", "transforms")])])

## (2) Subcortical segmentation via FIRST (fsl)
first = Node(fsl.FIRST(), name = "first")

wf.connect([
    (selectfiles, first, [("t1w", "in_file")])])

## (3) Create GM mask; Combine FIRST & dseg (thresholded at 0.75-1.25)
fsl_math_multi = Node(fsl.MultiImageMaths(), name="gm_mask")
fsl_math_multi.inputs.op_string = "-thr 0.75 -uthr 1.25 -bin -add %s -bin"
fsl_math_multi.inputs.out_file = "_gm_mask.nii.gz"

wf.connect([
    (selectfiles, fsl_math_multi, [("dseg", "in_file")]),
    (first, fsl_math_multi, [("segmentation_file", "operand_files")])])

## (4) Apply combined mask to atlas
applymask = Node(fsl.ApplyMask(), name = "applymask")
applymask.inputs.out_file = "_atlas_coreg_masked.nii.gz"

wf.connect([
    (applytransform, applymask, [("output_image", "in_file")]),
    (fsl_math_multi, applymask, [("out_file", "mask_file")])])

# Define output/naming structure
datasink = Node(DataSink(), name="datasink")
datasink.inputs.base_directory = "."
datasink.inputs.container = "."

substitutions = [("_id_", "sub-"),
                ("datasink", "")]
subjNaming = [("sub-%s/" % sub, "sub-%s" % sub) for sub in id_list]
substitutions.extend(subjNaming)
datasink.inputs.substitutions = substitutions

wf.connect([(applytransform, datasink, [("output_image", "atlas.@coreg")]),
            (fsl_math_multi, datasink, [("out_file", "anat.@combseg")]),
            (applymask, datasink, [("out_file", "atlas.@masked")])])

# Visualize the Workflow
wf.write_graph(join(path_derivatives, "graphs/flow_fmri_coreg_atlas.dot"))
wf.write_graph(join(path_derivatives, "graphs/flow_fmri_coreg_atlas.dot"), graph2use="flat")

# Set parallel computing
mp.set_start_method("fork") 
plugin_args = {'n_procs': n_procs,
               'memory_gb': memory_gb,
               'n_gpu': n_gpu,
               'ngpuproc': ngpuproc}

# Run the workflow
wf.run(plugin='MultiProc', plugin_args=plugin_args)

# Delete working directory
dir = os.listdir(f"{path_derivatives}/{wf_name}")
folders_to_delete = [item for item in dir if item.startswith("_")]
for folder in folders_to_delete:
    folder_path = os.path.join(f"{path_derivatives}/{wf_name}", folder)
    shutil.rmtree(folder_path)

print("The script has finished.")
