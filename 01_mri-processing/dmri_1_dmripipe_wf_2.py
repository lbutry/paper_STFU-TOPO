########################################################
# Author: Lionel Butry (@lbutry)
#
# Purpose: Iterate through all participants & perform the following workflows:
#  (2) Workflow: Tensor Reconstruction (MRtrix3)
#      - (A) dwi2response, (B) dwi2mask, (C) dwi2fod, 
#      - (D) dwi2tensor, tensor2metric, (E) mtnormalise
# 
# Requirements:
#   - dwi & t1w data in BIDS format
#   - Installation of 'mrtrix3', 'FSL', 'ANTs'
########################################################

# Import modules
import multiprocessing as mp, shutil
from nipype import Node, Workflow
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
    ## External input from wf_preproc 
    path_dwi      = join(derivatives_dir, name_of_wf, "{id}", "dwi", "sub-{id}_preproc_dwi_dstriped.mif")

    path_dict = {"dwi_preproc": path_dwi}

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
    # WF (2) Create Reconstruction Workflow
    ########################################################

    wf_recon = Workflow(name="wf_tensor_recon")

    # (0) Masking 
    dwi2mask = Node(mrt.BrainMask(), name="dwi2mask")
    dwi2mask.inputs.out_file="_dwi_mask.nii.gz"

    wf_recon.connect([
        (selectfiles, dwi2mask, [("dwi_preproc", "in_file")])])

    # (A) Estimate response functions for CSD
    dwi2response = Node(mrt.ResponseSD(), name="dwi2response")
    dwi2response.inputs.algorithm="dhollander"
    dwi2response.inputs.args=" -voxels voxels.mif"
    dwi2response.inputs.csf_file="_csf_respf.txt"
    dwi2response.inputs.gm_file="_gm_respf.txt"
    dwi2response.inputs.wm_file="_wm_respf.txt"

    wf_recon.connect([
        (selectfiles, dwi2response, [("dwi_preproc", "in_file")])
    ])

    # (B) Masking
    dilate_mask = Node(mrt.MaskFilter(), name="dilate_mask")
    dilate_mask.inputs.filter="dilate"
    dilate_mask.inputs.npass=4
    dilate_mask.inputs.out_file="_mask_dilate.mif"

    wf_recon.connect([
        (dwi2mask, dilate_mask, [("out_file", "in_file")])])

    # (C) Estimate CSD-tensor
    dwi2fod = Node(mrt.ConstrainedSphericalDeconvolution(), name="dwi2fod")
    dwi2fod.inputs.algorithm="msmt_csd"
    dwi2fod.inputs.wm_odf="_wm_fod.mif"
    dwi2fod.inputs.gm_odf="_gm_fod.mif"
    dwi2fod.inputs.csf_odf="_csf_fod.mif"

    wf_recon.connect([
        (selectfiles, dwi2fod, [("dwi_preproc", "in_file")]),
        (dilate_mask, dwi2fod, [("out_file", "mask_file")]),
        (dwi2response, dwi2fod, [("wm_file", "wm_txt"),
                                ("gm_file", "gm_txt"),
                                ("csf_file", "csf_txt")])])

    # (D) Obtain FA, RD, AD, MD-maps
    dwi2tensor = Node(mrt.FitTensor(), name="dwi2tensor")
    dwi2tensor.inputs.out_file="_dti.mif"

    tensor2metric = Node(mrt.TensorMetrics(), name="tensor2metric")
    tensor2metric.inputs.out_fa="_fa_map.mif"
    tensor2metric.inputs.out_ad="_ad_map.mif"
    tensor2metric.inputs.out_rd="_rd_map.mif"
    tensor2metric.inputs.out_adc="_md_map.mif"

    wf_recon.connect([
        (selectfiles, dwi2tensor, [("dwi_preproc", "in_file")]),
        (dilate_mask, dwi2tensor, [("out_file", "in_mask")]),
        (dwi2tensor, tensor2metric, [("out_file", "in_file")])])

    # (E) Intensity normalization
    mtnormalise = Node(mrt.MTNormalise(), name = "mtnormalise")
    mtnormalise.inputs.out_file_wm="_wm_fod_norm.mif"
    mtnormalise.inputs.out_file_gm="_gm_fod_norm.mif"
    mtnormalise.inputs.out_file_csf="_csf_fod_norm.mif"

    wf_recon.connect([
        (dwi2fod, mtnormalise, [("wm_odf", "wm_fod"),
                                ("gm_odf", "gm_fod"),
                                ("csf_odf", "csf_fod")]),
        (dwi2mask, mtnormalise, [("out_file", "mask")])])

    # connect wf_recon to wf_main
    wf_main.connect([      
        (wf_recon, datasink, [("dwi2mask.out_file", "dwi.@mask"),
                            ("dwi2tensor.out_file", "recon.@dti"),
                            ("tensor2metric.out_fa", "recon.@fa_map"),
                            ("tensor2metric.out_ad", "recon.@ad_map"),
                            ("tensor2metric.out_rd", "recon.@rd_map"),
                            ("tensor2metric.out_adc", "recon.@md_map"),
                            ("mtnormalise.out_file_wm", "recon.@wm_norm"),
                            ("mtnormalise.out_file_gm", "recon.@gm_norm"),
                            ("mtnormalise.out_file_csf", "recon.@csf_norm")])])

    ########################################################
    # (5) Run the Main Workflow
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
