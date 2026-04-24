########################################################
########################################################
# Author: Lionel Butry (@lbutry)
# Purpose: 
#  (1) Creates structural connectivity matrices based iFOD2 & SIFT2
# Requirements:
#   - Output of 'dmripipe'
#   - Installation of 'mrtrix3'
########################################################
########################################################

# Import modules
import os
import subprocess

###################### USER INPUT ######################
path_derivatives = "/Volumes/DataDrive/01_Scanner_Comparison/02_bids_root/derivatives"
path_dmripipe = os.path.join(path_derivatives, "d01_dmripipe") # dmri pipeline folder
output_folder = os.path.join(path_derivatives, "d02_connectome_prob")
n_proc = os.cpu_count() - 1 # Number of CPUs
########################################################

id_list = [] # list of subject IDs (e.g., ["PREDICTxxxxx", "PREDICTxxxxx", ...])

# Create output_dir
if not os.path.exists(output_folder):
  os.mkdir(output_folder)

def create_connectome(id):

  path_fod = f"{path_dmripipe}/{id}/recon/sub-{id}_wm_fod_norm.mif"
  path_tracks = f"{path_dmripipe}/{id}/recon/sub-{id}_tracks_ifod2.tck"
  path_sift2 = f"{output_folder}/sub-{id}_sift2_weights.txt"
  path_atlas  = f"{path_dmripipe}/{id}/atlas/sub-{id}_space_dwi_atlas.nii.gz"
  path_connectome = f"{output_folder}/sub-{id}_connectome_sift2.csv"
  path_assignments = f"{output_folder}/sub-{id}_assignments_sift2.csv"

  # Run SIFT2
  cmd = f"tcksift2 {path_tracks} {path_fod} {path_sift2}"
  subprocess.run(cmd, shell=True, capture_output=True)

  # Create connectome based on SIFT2
  cmd = f"tck2connectome {path_tracks} {path_atlas} {path_connectome} -tck_weights_in {path_sift2} -out_assignments {path_assignments} -symmetric -zero_diagonal"
  subprocess.run(cmd, shell=True, capture_output=True)

  print("Finished connectome construction for: ", id)

# Run 'create_connectome' for each subject
for id in id_list:
  create_connectome(id)
