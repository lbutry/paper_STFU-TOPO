# Computes functional connectivity matrices for each subject (Fisher-z transformed pearson's r)
# Author: Lionel Butry (@lbutry)

from bids import BIDSLayout
from nilearn.maskers import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
import pandas as pd
import numpy as np
import os

###################### USER INPUT ######################
path_fmriprep = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/f01_fmriprep"
path_nuireg     = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/f02_nilearn_nuireg"
path_coreg_atlas = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/f03_coreg_atlas"
labels = pd.read_csv("/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/utils/BN_Atlas_246_LUT.tsv", sep="\t")["ROI.Name"]
output_dir = "/Volumes/DataDrive2/00_PREDICT-LBP/bids_root/derivatives/f04_connectome_partial"
n_proc = os.cpu_count() - 1
########################################################

# Create output_dir
os.makedirs(output_dir, exist_ok=True)

# Define subjects
layout = BIDSLayout(path_nuireg, validate=False)
id_list = layout.get_subjects()

for sub in id_list:

    # Create masker to extract BOLD data within (individualized) atlas parcels
    path = f"{path_coreg_atlas}/sub-{sub}/atlas/sub-{sub}_atlas_coreg_masked.nii.gz"
    masker = NiftiLabelsMasker(labels_img=path, standardize="zscore_sample")

    # Extract time series
    path_bold = f"{path_nuireg}/sub-{sub}_scrubbing5_bold.nii.gz"
    time_series = masker.fit_transform(path_bold)

    # Define the connectivity measure (here: partial correlation)
    conn = ConnectivityMeasure(kind="partial correlation", standardize=False)
    mat = conn.fit_transform([time_series])[0]
    np.fill_diagonal(mat, 0) # set diagonal to 0

    # Save as .csv
    mat = pd.DataFrame(mat, index=labels, columns=labels)
    mat.to_csv(f"{output_dir}/sub-{sub}_connectome_func_partial.csv")
    print(f"The FC of {sub} has been stored in: {output_dir}")
