########################################################
# Author: Lionel Butry (@lbutry)
#
# Purpose: Iterate through all participants & perform the following workflows:
#  (1) Workflow: DWI preprocessing (MRtrix3)
#      - (A) mrconvert, (B) dwidenoise, (C) mrdegibbs, 
#      - (D) dwifslpreproc, (E) dwibiascorrect
#
# Requirements:
#   - dwi & t1w data in BIDS format
#   - Installation of 'mrtrix3', 'FSL', 'ANTs'
########################################################

# Import modules
import multiprocessing as mp
import shutil
from nipype import Node, Workflow, Merge
from nipype.interfaces import mrtrix3 as mrt
from nipype.interfaces import DataSink, IdentityInterface, SelectFiles
from os.path import join

########################################################
# Create Main Workflow
########################################################

def run_dmripipe(id_list, bids_root, derivatives_dir, work_dir, name_of_wf, n_procs, memory_gb, n_gpu, ngpuproc, delete_workdir=True):

    wf_main = Workflow(name=name_of_wf, base_dir=work_dir)

    # Iterate over all subs
    idsource = Node(IdentityInterface(fields=['id']), name="idsource")
    idsource.iterables = [('id', id_list)]

    # Select files
    ## External input for wf_preproc 
    path_dwi_json = join("sub-{id}", "dwi", "sub-{id}_acq-105_dir-AP_dwi.json")
    path_dwi      = join("sub-{id}", "dwi", "sub-{id}_acq-105_dir-AP_dwi.nii.gz")
    path_bval     = join("sub-{id}", "dwi", "sub-{id}_acq-105_dir-AP_dwi.bval")
    path_bvec     = join("sub-{id}", "dwi", "sub-{id}_acq-105_dir-AP_dwi.bvec")
    path_dwi_rev  = join("sub-{id}", "dwi", "sub-{id}_acq-rev_dir-PA_dwi.nii.gz")
    path_bval_rev = join("sub-{id}", "dwi", "sub-{id}_acq-rev_dir-PA_dwi.bval")
    path_bvec_rev = join("sub-{id}", "dwi", "sub-{id}_acq-rev_dir-PA_dwi.bvec")

    path_dict = {"dwi": path_dwi,
                "dwi_rev": path_dwi_rev,
                "dwi_json": path_dwi_json,
                "bval": path_bval,
                "bvec": path_bvec,
                "bval_rev": path_bval_rev,
                "bvec_rev": path_bvec_rev}

    selectfiles = Node(SelectFiles(templates=path_dict, base_directory=bids_root), name="selectfiles")

    # Define output naming structure
    datasink = Node(DataSink(), name="datasink")
    datasink.inputs.base_directory = join(derivatives_dir, name_of_wf)
    datasink.inputs.substitutions = [("_id_%s/" % sub, "sub-%s" % sub) for sub in id_list]
    datasink.inputs.remove_dest_dir = True

    wf_main.connect([
        (idsource, selectfiles, [("id", "id")]),
        (idsource, datasink, [("id", "container")])])

    ########################################################
    # (1) Create Preprocessing Workflow
    ########################################################

    wf_preproc = Workflow(name="wf_dwi_preproc")

    # (A) Convert .nii.gz to .mif files
    mrconvert_dwi = Node(mrt.MRConvert(), name="01_mrconvert_dwi") 
    mrconvert_dwi.inputs.out_file="_dwi.mif"

    mrconvert_rev = Node(mrt.MRConvert(), name="02_mrconvert_rev") 
    mrconvert_rev.inputs.out_file="_dwi_rev.mif"

    # (B) Denoising
    dwidenoise = Node(mrt.DWIDenoise(), name="03_dwidenoise")
    dwidenoise.inputs.noise="_noise.mif"

    wf_preproc.connect([
        (mrconvert_dwi, dwidenoise, [("out_file", "in_file")])])

    # (C) Gibbs ringing correction
    mrdegibbs = Node(mrt.MRDeGibbs(), name="04_mrdegibbs")
    mrdegibbs.inputs.axes=[0, 1]

    wf_preproc.connect([
        (dwidenoise, mrdegibbs, [("out_file", "in_file")])])

    # (D) Motion, Eddy & inhomogeneity correction 
    # (D.1) Create b0 pair for eddy
    # (D.1.1) Extract b0 images from AP & PA dwi image
    dwiextract_AP = Node(mrt.DWIExtract(), name="05_dwiextract_AP")
    dwiextract_AP.inputs.bzero=True
    dwiextract_AP.inputs.out_file="_b0_AP.mif"

    dwiextract_PA = Node(mrt.DWIExtract(), name="06_dwiextract_PA")
    dwiextract_PA.inputs.bzero=True
    dwiextract_PA.inputs.out_file="_b0_PA.mif"

    wf_preproc.connect([
        (mrdegibbs, dwiextract_AP, [("out_file", "in_file")]),
        (mrconvert_rev, dwiextract_PA, [("out_file", "in_file")])])

    # (D.1.2) Estimate b0 mean image for AP & PA dwi images separately
    b0_mean_AP = Node(mrt.MRMath(), name="07_b0_mean_AP")
    b0_mean_AP.inputs.axis=3
    b0_mean_AP.inputs.operation="mean"
    b0_mean_AP.inputs.out_file="_b0_AP_mean.mif"

    b0_mean_PA = Node(mrt.MRMath(), name="08_b0_mean_PA")
    b0_mean_PA.inputs.axis=3
    b0_mean_PA.inputs.operation="mean"
    b0_mean_PA.inputs.out_file="_b0_PA_mean.mif"

    wf_preproc.connect([
        (dwiextract_AP, b0_mean_AP, [("out_file", "in_file")]),
        (dwiextract_PA, b0_mean_PA, [("out_file", "in_file")])])

    # (D.1.3) Concatenate b0 AP/PA mean images into one file
    merge = Node(Merge(2), name="09_merge") # mrcat needs the input merged

    mrcat = Node(mrt.MRCat(), name="10_mrcat")
    mrcat.inputs.axis=3
    mrcat.inputs.out_file="_b0_pair.mif"

    wf_preproc.connect([
        (b0_mean_AP, merge, [("out_file", "in1")]),
        (b0_mean_PA, merge, [("out_file", "in2")]),
        (merge, mrcat, [("out", "in_files")])])

    # (D.2) Topup & Eddy
    dwipreproc = Node(mrt.DWIPreproc(), name="11_dwipreproc")
    dwipreproc.inputs.rpe_options="pair"
    dwipreproc.inputs.align_seepi=True
    dwipreproc.inputs.eddy_options=" --repol --data_is_shelled --slm=linear"
    dwipreproc.inputs.pe_dir="AP"
    dwipreproc.inputs.out_grad_mrtrix="_grad.b"
    dwipreproc.inputs.out_file="_topup_eddy.mif"

    wf_preproc.connect([
        (mrdegibbs, dwipreproc, [("out_file", "in_file")]),
        (mrcat, dwipreproc, [("out_file", "in_epi")])])

    # (E) Bias field correction
    dwibiascorrect = Node(mrt.DWIBiasCorrect(), name="12_dwibiascorrect")
    dwibiascorrect.inputs.use_ants=True
    dwibiascorrect.inputs.out_file="_preproc_dwi.mif"

    wf_preproc.connect([
        (dwipreproc, dwibiascorrect, [("out_file", "in_file")])])


    # connect wf_preproc to wf_main
    wf_main.connect([ 
        (selectfiles, wf_preproc, [("dwi", "01_mrconvert_dwi.in_file"),
                                ("bvec", "01_mrconvert_dwi.in_bvec"),
                                ("bval", "01_mrconvert_dwi.in_bval"),
                                ("dwi_rev", "02_mrconvert_rev.in_file"),
                                ("bvec_rev", "02_mrconvert_rev.in_bvec"),
                                ("bval_rev", "02_mrconvert_rev.in_bval"),
                                ("dwi_json", "11_dwipreproc.json_import")]), 
        (wf_preproc, datasink, [("12_dwibiascorrect.out_file", "dwi.@final_preproc")])])

    ########################################################
    # (2) Run the Main Workflow
    ########################################################

    # Visualize the Main Workflow
    wf_main.write_graph(dotfilename=join(derivatives_dir, name_of_wf, f"graphs/flow_{name_of_wf}_color.dot"),
                        graph2use="colored")
    wf_main.write_graph(dotfilename=join(derivatives_dir, name_of_wf, f"graphs/flow_{name_of_wf}.dot"),
                        graph2use="flat")

    # Set parallel computing (set only when it has not been set before; important for sequential batches)
    if mp.get_start_method(allow_none=True) is None: 
        mp.set_start_method("fork")

    plugin_args = {'n_procs': n_procs,
                'memory_gb': memory_gb,
                'n_gpu': n_gpu,
                'ngpuproc': ngpuproc}

    # Run the combined workflow
    wf_main.run(plugin='MultiProc', plugin_args=plugin_args)

    # Delete working directory if 'delete_workdir == True'
    if delete_workdir: 
        shutil.rmtree(work_dir)
