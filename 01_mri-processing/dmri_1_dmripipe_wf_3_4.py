########################################################
# Author: Lionel Butry (@lbutry)
#
# Purpose: Iterate through all participants & perform the following workflows:
#  (3) Workflow: Coregistration (MRtrix3, FSL, ANTs)
#      - (A) fod2dec, mrmath, FLIRT, robustfov, n4bias, bet, antsRegistrationSynQuick, applytransforms,
#      - (B) 5ttgen, 5tt2gmwmi,
#      - (C) antsRegistrationSynQuick, applytransforms
#   (4) Node: Tractography (MRtrix3)
#      - tckgen
# 
# Requirements:
#   - dwi & t1w data in BIDS format
#   - Installation of 'mrtrix3', 'FSL', 'ANTs'
########################################################

# Import modules
import multiprocessing as mp, shutil
from nipype import Node, Workflow, Merge
from nipype.interfaces import mrtrix3 as mrt
from nipype.interfaces import DataSink, IdentityInterface, SelectFiles, fsl, ants
from os.path import join

########################################################
# Create Main Workflow
########################################################

def run_dmripipe(id_list, bids_root, derivatives_dir, work_dir, name_of_wf, n_procs, memory_gb, n_gpu, ngpuproc, delete_workdir=True, output_transforms=False):
        
    wf_main = Workflow(name=name_of_wf, base_dir=work_dir)

    # Iterate over all subs
    idsource = Node(IdentityInterface(fields=['id']), name="idsource")
    idsource.iterables = [('id', id_list)]

    # Select files
    ## External input for wf_coreg
    path_t1w      = join("sub-{id}", "anat", "sub-{id}_acq-anat_T1w.nii.gz")
    path_mni      = join("derivatives", "utils", "spm152.nii.gz")
    path_atlas    = join("derivatives", "utils", "BN_Atlas_246_1mm.nii")
    path_fod_wm   = join(derivatives_dir, f"{name_of_wf}", "{id}", "recon", "sub-{id}_wm_fod_norm.mif")
    path_mask     = join(derivatives_dir, f"{name_of_wf}", "{id}", "dwi", "sub-{id}_dwi_mask.nii.gz")

    ## External input for wf_tracto
    path_dict = {"fod_wm": path_fod_wm,
                "mask": path_mask,
                "t1w": path_t1w,
                "mni": path_mni,
                "atlas": path_atlas}

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
    # (3) Create Coregistration Workflow
    ########################################################

    wf_coreg = Workflow(name="wf_coreg")

    ## (A) Coregister T1w -> DWI
    ## (A.1) Create DEC-map (resampled to 1x1x1 & cropped to brain-mask)
    fod2dec = Node(mrt.FOD2Dec(), name="fod2dec")
    fod2dec.inputs.out_file="_dec.mif"

    dec_mean = Node(mrt.MRMath(), name="dec_mean")
    dec_mean.inputs.operation="mean"
    dec_mean.inputs.axis=3
    dec_mean.inputs.out_file="_dec_mean.nii.gz"

    apply_mask2dec = Node(fsl.ApplyMask(), name="apply_mask2dec")
    apply_mask2dec.inputs.out_file = "_dec_mean_masked.nii.gz"

    flirt_resamp = Node(fsl.FLIRT(), name="flirt_resamp")
    flirt_resamp.inputs.apply_isoxfm=1
    flirt_resamp.inputs.out_file="_dec_mean_resamp.nii.gz"

    wf_coreg.connect([
        (selectfiles, fod2dec, [("fod_wm", "in_file")]),
        (fod2dec, dec_mean, [("out_file", "in_file")]),
        (dec_mean, apply_mask2dec, [("out_file", "in_file")]),
        (selectfiles, apply_mask2dec, [("mask", "mask_file")]),
        (apply_mask2dec, flirt_resamp, [("out_file", "in_file"),
                                        ("out_file", "reference")])])

    ## (A.2) Preproc t1w
    robustfov = Node(fsl.RobustFOV(), name="robustfov") # crop the t1w
    robustfov.inputs.out_roi="_t1w_crop.nii.gz"

    n4bias = Node(ants.N4BiasFieldCorrection(), name="n4bias")
    n4bias.inputs.output_image="_t1w_crop_bias.nii.gz"

    bet = Node(fsl.BET(), name="bet")
    bet.inputs.out_file="_t1w_crop_bias_bet.nii.gz"
    bet.inputs.robust=True

    wf_coreg.connect([
        (robustfov, n4bias, [("out_roi", "input_image")]),
        (n4bias, bet, [("output_image", "in_file")])])

    ## (A.3) ants registration (dwi2t1w)
    antsreg_dwi2t1w = Node(ants.RegistrationSynQuick(), name="antsreg_dwi2t1w")
    antsreg_dwi2t1w.inputs.transform_type="b" # rigid + affine + deformable b-spline syn
    antsreg_dwi2t1w.inputs.output_prefix="_dwi2t1w_"

    wf_coreg.connect([
        (bet, antsreg_dwi2t1w, [("out_file", "fixed_image")]),
        (flirt_resamp, antsreg_dwi2t1w, [("out_file", "moving_image")])])

    ## (A.4) apply inverse transforms (t1w2dwi)
    merge_transforms = Node(Merge(2), name="merge_transforms") # necessary for apply_invtransforms.inputs.transforms

    apply_invtransforms = Node(ants.ApplyTransforms(), name="apply_invtransforms")
    apply_invtransforms.inputs.dimension=3
    apply_invtransforms.inputs.invert_transform_flags=[False, True]
    apply_invtransforms.inputs.output_image="_space-dwi_t1w.nii.gz"

    wf_coreg.connect([
        (antsreg_dwi2t1w, merge_transforms, [("inverse_warp_field", "in1"),
                                            ("out_matrix", "in2")]),
        (merge_transforms, apply_invtransforms, [("out", "transforms")]),
        (n4bias, apply_invtransforms, [("output_image", "input_image")]),
        (flirt_resamp, apply_invtransforms, [("out_file", "reference_image")])])

    ## (B) 5TT Segmentation & GMWMI
    gen_5tt = Node(mrt.Generate5tt(), name="5ttgen")
    gen_5tt.inputs.algorithm="fsl"
    gen_5tt.inputs.sgm_amyg_hipp=True
    gen_5tt.inputs.nocrop=True
    gen_5tt.inputs.out_file="_space-dwi_5tt.mif"

    gmwmi = Node(mrt.Generate5tt2gmwmi(), name="5tt2gmwmi")
    gmwmi.inputs.mask_out="_space-dwi_gmwmi.mif"

    wf_coreg.connect([
        (apply_invtransforms, gen_5tt, [("output_image", "in_file")]),
        (gen_5tt, gmwmi, [("out_file", "in_file")])])

    # (C) Coregister Atlas -> DWI (DEC map)
    ## (C.1) ants registration (mni2t1w)
    antsreg_mni2t1w = Node(ants.RegistrationSynQuick(), name="antsreg_mni2t1w")
    antsreg_mni2t1w.inputs.transform_type="b" # rigid + affine + deformable b-spline syn
    antsreg_mni2t1w.inputs.output_prefix="_mni2t1w_"

    wf_coreg.connect([
        (bet, antsreg_mni2t1w, [("out_file", "fixed_image")])])

    ## (C.2) Apply inverse transformation (mni2t1w + t1w2dwi)
    merge_transforms_mni2dwi = Node(Merge(4), name="merge_transforms_mni2dwi") # necessary for applyreg_atlas.inputs.transforms

    applytransforms_mni2dwi = Node(ants.ApplyTransforms(), name="applytransforms_mni2dwi")
    applytransforms_mni2dwi.inputs.dimension=3
    applytransforms_mni2dwi.inputs.invert_transform_flags=[False, True, False, False]
    applytransforms_mni2dwi.inputs.interpolation="NearestNeighbor"
    applytransforms_mni2dwi.inputs.output_image="_space-dwi_atlas.nii.gz"

    wf_coreg.connect([
        (antsreg_dwi2t1w, merge_transforms_mni2dwi, [("inverse_warp_field", "in1"),
                                                    ("out_matrix", "in2")]),
        (antsreg_mni2t1w, merge_transforms_mni2dwi, [("forward_warp_field", "in3"),
                                                    ("out_matrix", "in4")]),
        (merge_transforms_mni2dwi, applytransforms_mni2dwi, [("out", "transforms")]),
        (flirt_resamp, applytransforms_mni2dwi, [("out_file", "reference_image")])])

    # connect wf_coreg to wf_main
    wf_main.connect([   
        (selectfiles, wf_coreg, [("t1w", "robustfov.in_file"),
                                ("mni", "antsreg_mni2t1w.moving_image"),
                                ("atlas", "applytransforms_mni2dwi.input_image")]),
        (wf_coreg, datasink, [("fod2dec.out_file", "recon.@dec"),
                            ("flirt_resamp.out_file", "recon.@dec_resamp"),
                            ("apply_invtransforms.output_image", "anat.@coreg"),
                            ("5ttgen.out_file", "anat.@5tt"),
                            ("5tt2gmwmi.mask_out", "anat.@gmwmi"),
                            ("applytransforms_mni2dwi.output_image", "atlas.@coreg")])])

    if output_transforms: 
        wf_main.connect([
            (wf_coreg, datasink, [("antsreg_dwi2t1w.out_matrix", "transforms.@matrix"),
                                ("antsreg_dwi2t1w.inverse_warp_field", "transforms.@invwarp"),
                                ("antsreg_mni2t1w.out_matrix", "transforms.@matrix2"),
                                ("antsreg_mni2t1w.inverse_warp_field", "transforms.@invwarp2")])
        ])

    ########################################################
    #(4) Create Tractography Node (iFOD2 + ACT + GMWMI seeding)
    ########################################################

    tckgen = Node(mrt.Tractography(), name="tckgen")
    tckgen.inputs.algorithm="iFOD2"
    tckgen.inputs.select=5000000 #5M streamlines
    tckgen.inputs.step_size=1 # 0.5 x voxelsize (default)
    tckgen.inputs.angle=45 # default iFOD2
    tckgen.inputs.min_length=4 # 2 x voxelsize for ACT
    tckgen.inputs.max_length=250 # default = 200
    tckgen.inputs.cutoff=0.05 # 0.1 * 0.5 because ACT
    tckgen.inputs.crop_at_gmwmi=True
    tckgen.inputs.backtrack=True
    tckgen.inputs.out_file="_tracks_ifod2.tck"

    wf_main.connect([  
        (selectfiles, tckgen, [("fod_wm", "in_file")]),
        (wf_coreg, tckgen, [("5ttgen.out_file", "act_file"),
                            ("5tt2gmwmi.mask_out", "seed_gmwmi")]),
        (tckgen, datasink, [("out_file", "recon.@tracks2")])
    ])

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
