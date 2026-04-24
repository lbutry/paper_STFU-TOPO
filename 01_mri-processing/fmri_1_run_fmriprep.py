########################################################
########################################################
# Author: Lionel Butry (@lbutry)
# Purpose:
#   (1) Runs fMRIPrep using Docker
# Requirements:
#   - fMRI data in BIDS format
#   - Installation of 'fmriprep' (docker) & 'fmriprep-docker' (docker wrapper)
#   - Freesurfer license in .txt
########################################################
# Note: Recommended to set working directory to local directory and output to external drive.
#       Work dir is deleted after each successful Docker run.
########################################################

# Import modules
import os 
import subprocess
import sys
import datetime
import time
import shutil

########################################################
###################### USER INPUT ######################
# Subjects to preprocess
manual_id_list = None # define subjects to be processed: ['PREDICT12345', ...]; use None to process all subjects in the BIDS directory in batches
batch_size = 10 # Number of subjects to process per Docker run
exclude_id_list = [] # List of subject IDs to exclude from processing, e.g. ['PREDICT12345', 'PREDICT12346']; use None to include all subjects
# Options
rm_work_dir = True # Remove work directory after each Docker run (recommended)
# Paths
bids_root = "/media/team-tesla-linux/DataDrive2/00_PREDICT-LBP/bids_root" # Path to BIDS directory
output_folder_name = "f01_fmriprep" # Name of the output folder
output_dir = f"{bids_root}/derivatives/{output_folder_name}" # Path to output directory
work_dir = "/home/team-tesla-linux/Documents/Lionel/fmriprep_work" # Path to work directory
freesurfer_txt = "/home/team-tesla-linux/Documents/Lionel/freesurfer.txt" # freesurfer license
# Computational parameteres
nthreads = os.cpu_count() - 2 # maximum number of CPUs
memory_mb = 65000 # maximum MB of memory

########################################################
########################################################

########################################################
# Functions
########################################################

def get_unprocessed_subjects(bids_root, deri, num_subs, exclude_id=None):
    """
    Returns a list of subject IDs from the BIDS directory that do not have corresponding derivative outputs.
        
    Args:
    - bids_root (str): Path to the BIDS root directory.
    - deri (str): Path to the derivative directory.
    - num_subs (int): Number of subjects IDs
    """

    # all subject folders in BIDS root
    bids_id_list = [
        d.removeprefix('sub-') for d in os.listdir(bids_root) # remove 'sub-' prefix
        if d.startswith('sub-') and os.path.isdir(os.path.join(bids_root, d))
    ]

    # get IDs from HTML-reports (final output of fMRIPrep)
    deri_html_id = [
        f.removeprefix('sub-').removesuffix('.html') for f in os.listdir(deri)
        if f.startswith('sub-') and f.endswith('.html') and os.path.isfile(os.path.join(deri, f))
    ]

    # exclude a priori specified IDs
    bids_id_set = set(bids_id_list)
    if exclude_id is not None:
        bids_id_set -= set(exclude_id)

    # subjects in BIDS but not in derivatives
    unprocessed_ids = bids_id_set - set(deri_html_id)

    return sorted(unprocessed_ids)[:num_subs]

def get_runtime():
    end_time = datetime.datetime.now()
    run_time = end_time - start_time
    run_time = str(run_time).split('.')[0] + ' (HH:MM:SS)'
    return run_time

########################################################
# Main script
########################################################

# Initialize start time
start_time = datetime.datetime.now()

# Set id_list based on user input
if manual_id_list is not None: # Set id_list to manual_id_list
    id_list = manual_id_list
else: # Initialize first batch of subjects
    id_list = get_unprocessed_subjects(bids_root, output_dir, batch_size, exclude_id_list)

if not id_list: 
    raise ValueError("No unprocessed subjects found in the BIDS directory or no IDs provided manually.")

batch = 1 # Initialize batch counter

# Create output_dir (if it does not exist)
os.makedirs(output_dir, exist_ok=True)

# Loop Docker as long as there are subjects to process (sequentially for each batch)
while id_list:
    print(f"=== Running fMRIPrep for current batch: {id_list} ===")

    docker_command = [ # no freesurfer, t1w & MNI output space
        "fmriprep-docker",
        bids_root,
        output_dir,
        "participant",
        "--participant-label", *map(str, id_list),
        "--skull-strip-t1w", "force",
        "--fs-license-file", freesurfer_txt,
        "--fs-no-reconall",
        "--output-spaces", "MNI152NLin2009cAsym:res-native", "anat",
        "--nthreads", str(nthreads),
        "--mem_mb", str(memory_mb),
        "--low-mem",
        "-w", work_dir
    ]

    try:
        subprocess.run(docker_command, check=True) # Run the Docker command
        time.sleep(20) # Sleep for 20 seconds to ensure Docker has time to clean up

        # Remove work directory if specified
        if rm_work_dir: 
            shutil.rmtree(work_dir)

        # Update id_list for the next batch or exit script
        if manual_id_list is None:
            id_list = get_unprocessed_subjects(bids_root, output_dir, batch_size, exclude_id_list) 
        else:
            id_list = [] # Exit the loop if manual_id_list is provided

        batch += 1 # Update batch counter

    except subprocess.CalledProcessError as e:
        
        print(f"*** Error running fMRIPrep for batch {batch} ({id_list}): {e} ***")
        run_time = get_runtime()
        sys.exit(1)

run_time = get_runtime()
body = f"Successfully completed fMRIPrep processing for all batches.\nTotal batches processed: {batch - 1}\nRuntime: {run_time}\nOutput directory: {output_dir}"

# Print final message to console
print("=== ", body, " ===") 