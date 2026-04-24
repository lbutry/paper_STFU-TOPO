########################################################
########################################################
# Author: Lionel Butry (@lbutry)
# Purpose: Performs nuisance regression using 'nilearn' on fmriprep-preprocessed data
#   - Detrending
#   - High & low pass filtering
#   - Scrubbing
#   - Removal of confounds (based on Wang et al. (2024): 'scrubbing.5')
#   - Normalization (zscore_sample)
# Requirements:
#   - Output path of fmriprep
########################################################
########################################################

# Import modules
import nibabel as nib, os, multiprocessing
from bids import BIDSLayout
from nilearn.interfaces.fmriprep import load_confounds_strategy
from nilearn.image import clean_img

###################### USER INPUT ######################
path_fmriprep   = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/f01_fmriprep" # fmriprep's output folder
output_dir      = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/f02_nilearn_nuireg" # output directory
n_proc = os.cpu_count() - 1
total_volumes = 190
########################################################

# Create output_dir
if not os.path.exists(output_dir):
  os.mkdir(output_dir)

# Import data
layout = BIDSLayout(path_fmriprep, validate=False)

subs          = layout.get(return_type="id", target="subject")
paths_bold    = layout.get(suffix="bold", extension=".nii.gz", space="T1w", return_type="file", task="rest")
paths_mask    = layout.get(suffix="mask", extension=".nii.gz", space="T1w", return_type="file", task="rest")

# Nuisance regression
## Confounds strategy corresponds to 'scrubbing.5' in Wang et al. (2024)
confounds, sample_mask = load_confounds_strategy(
  img_files = paths_bold,
  denoise_strategy = "scrubbing",
  fd_threshold = 0.5,
  std_dvars_threshold = None
)

## Create function for nuisance regression & saving results
def nuisance_regression(subj, bold, mask, confounds, sample_mask):

  # clean BOLD
  clean_bold = clean_img(
    imgs = bold,
    detrend = False, # included in confounds strategy (high_pass; cosines) # https://github.com/SIMEXP/load_confounds
    standardize = "zscore_sample", # (new) default
    confounds = confounds,
    low_pass= 0.1,
    high_pass = None, # already included in confounds strategy at 0.008 Hz  (high_pass; cosines)
    t_r = 2.5,
    mask_img = mask,
    clean__sample_mask = sample_mask # passed to signal.clean; used for scrubbing
  )

  # check no. of volumes < 70 % of total volumes & save as nifti
  min_volumes = 0.7 * total_volumes
  volumes = clean_bold.shape[3]

  if volumes >= min_volumes:
    nib.save(clean_bold, f"{output_dir}/sub-{subj}_scrubbing5_bold.nii.gz")
  else:
    nib.save(clean_bold, f"{output_dir}/sub-{subj}_scrubbing5_bold_less70%VolumesRemain.nii.gz")

  print("Finished preprocessing for:", subj)

# Run nuisance regression
if __name__ == "__main__":  
  with multiprocessing.Pool(processes=n_proc) as pool:
      pool.starmap(nuisance_regression, zip(subs, paths_bold, paths_mask, confounds, sample_mask))
