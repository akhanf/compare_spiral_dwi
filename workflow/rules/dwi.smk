wildcard_constraints:
    shell='[0-9]+',
    ndirs='[0-9]+'

rule dwi2mif:
    input:
        dwi=lambda wildcards: config["input_path"]["dwi_nii"][wildcards.dataset],
        bval=lambda wildcards: re.sub(
            ".nii.gz", ".bval", config["input_path"]["dwi_nii"][wildcards.dataset]
        ),
        bvec=lambda wildcards: re.sub(
            ".nii.gz", ".bvec", config["input_path"]["dwi_nii"][wildcards.dataset]
        ),
    output:
        dwi=bids(
            root=root,
            datatype="dwi",
            desc='init',
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrconvert {input.dwi} {output.dwi} -fslgrad {input.bvec} {input.bval} -nthreads {threads}"

rule split_into_shells:
    input:
        dwi=bids(
            root=root,
            datatype="dwi",
            desc='init',
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
    output:
        dwi=bids(
            root=root,
            datatype="dwi",
            desc='shell{shell}',
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
    params:
        opts=lambda wildcards:  '-shells 0' if wildcards.shell=='0' else f'-no_bzero -shells {wildcards.shell}'
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwiextract {params.opts} {input.dwi} {output.dwi}"

rule subsample_shell:
    input:
        dwi=bids(
            root=root,
            datatype="dwi",
            desc='shell{shell}',
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
    output:
        dwi=bids(
            root=root,
            datatype="dwi",
            desc='shell{shell}subsampled{ndirs}',
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
    params:
        opts=lambda wildcards: '-coord 3 0:1:{n}'.format(n=int(wildcards.ndirs)-1)
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrconvert {input.dwi} {output.dwi}  {params.opts}"

   
rule combine_shells:
    input:
        #b0, b1000 subsampled, b2000 subsampled
        dwi_shells=lambda wildcards: expand(bids(root=root,
            datatype="dwi",
            desc='{desc}',
            suffix="dwi.mif",
            **config["subj_wildcards"],
            ),desc=config['subsampling'][wildcards.dataset],allow_missing=True)
   
    output:
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrcat {input} {output}"

         

rule get_avgb0:
    input:
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),        
    output:
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="avgb0.nii.gz",
            **config["subj_wildcards"],
        ),        
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        'dwiextract {input} - -bzero | mrmath - mean {output} -axis 3'


              
rule reg_b0_to_t1:
    input:
        avgb0=bids(
            root=root,
            datatype="dwi",
            suffix="avgb0.nii.gz",
            **config["subj_wildcards"],
        ),        
        t1=lambda wildcards: config["input_path"]["t1_nii"],
    params:
        general_opts="-d 3",
        rigid_opts="-m NMI -a -dof 6 -ia-identity -n 50x50"
    
    output:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_='b0',
            to='T1w',
            type_='ras',
            suffix="affine.txt",
            **config["subj_wildcards"],
        ),        
        warped_avgb0=bids(
            root=root,
            datatype="dwi",
            space='T1w',
            suffix="avgb0.nii.gz",
            **config["subj_wildcards"],
        ),        
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} {params.rigid_opts}"
        "    -i {input.t1} {input.avgb0} -o {output.xfm_ras} && "
        "greedy -threads {threads} {params.general_opts} -rf {input.t1} "
        "    -rm {input.avgb0} {output.warped_avgb0} -r {output.xfm_ras}"
   
 
rule warp_t1_to_dwi:
    """for overlay purposes only"""
    input:
        t1=lambda wildcards: config["input_path"]["t1_nii"],
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_='b0',
            to='T1w',
            type_='ras',
            suffix="affine.txt",
            **config["subj_wildcards"],
        ),        
        avgb0=bids(
            root=root,
            datatype="dwi",
            suffix="avgb0.nii.gz",
            **config["subj_wildcards"],
        ),        
    params:
        general_opts="-d 3",
 
    output:
        warped_t1=bids(
            root=root,
            datatype="dwi",
            space='dwi',
            suffix="T1w.nii.gz",
            **config["subj_wildcards"],
        ),
    shell:
        "greedy -threads {threads} {params.general_opts} -rf {input.avgb0} "
        "    -rm {input.t1} {output.warped_t1} -r {input.xfm_ras},-1"



rule fivettgen:
    input:
        aseg=config['input_path']['fs_aseg']
    output:
        fivett=bids(
            root=root,
            datatype="dwi",
            desc="freesurfer",
            suffix="5tt.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shadow: 'minimal'
    shell:
        "5ttgen freesurfer {input} {output} -nthreads {threads}"


rule dwi2response:
    input:
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
        fivett=rules.fivettgen.output.fivett,
        mask=lambda wildcards: config["input_path"]["dwi_mask"][wildcards.dataset],
    params:
        shells=",".join(config["dwi"]["shells"]),
        lmax=",".join(config["dwi"]["lmax"]),
    output:
        wm_rf=bids(
            root=root,
            datatype="dwi",
            desc="wm",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
        gm_rf=bids(
            root=root,
            datatype="dwi",
            desc="gm",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
        csf_rf=bids(
            root=root,
            datatype="dwi",
            desc="csf",
            suffix="response.txt",
            **config["subj_wildcards"],
        ),
        voxels=bids(
            root=root,
            datatype="dwi",
            suffix="responsevoxels.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2response msmt_5tt {input.dwi} {input.fivett} {output.wm_rf} {output.gm_rf} {output.csf_rf}  -nthreads {threads} -shells {params.shells} -lmax {params.lmax} -mask {input.mask} -voxels {output.voxels}"


rule dwi2fod:
    # Jeurissen, B; Tournier, J-D; Dhollander, T; Connelly, A & Sijbers, J. Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data. NeuroImage, 2014, 103, 411-426
    input:
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),

        mask=lambda wildcards: config["input_path"]["dwi_mask"][wildcards.dataset],
        wm_rf=rules.dwi2response.output.wm_rf,
        gm_rf=rules.dwi2response.output.gm_rf,
        csf_rf=rules.dwi2response.output.csf_rf,
    params:
        shells=",".join(config["dwi"]["shells"]),
    output:
        wm_fod=bids(
            root=root,
            datatype="dwi",
            desc="wm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids(
            root=root,
            datatype="dwi",
            desc="gm",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids(
            root=root,
            datatype="dwi",
            desc="csf",
            suffix="fod.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2fod -nthreads {threads} -mask {input.mask} -shell {params.shells} msmt_csd {input.dwi} {input.wm_rf} {output.wm_fod} {input.gm_rf} {output.gm_fod} {input.csf_rf} {output.csf_fod} "


rule mtnormalise:
    # Raffelt, D.; Dhollander, T.; Tournier, J.-D.; Tabbara, R.; Smith, R. E.; Pierre, E. & Connelly, A. Bias Field Correction and Intensity Normalisation for Quantitative Analysis of Apparent Fibre Density. In Proc. ISMRM, 2017, 26, 3541
    # Dhollander, T.; Tabbara, R.; Rosnarho-Tornstrand, J.; Tournier, J.-D.; Raffelt, D. & Connelly, A. Multi-tissue log-domain intensity and inhomogeneity normalisation for quantitative apparent fibre density. In Proc. ISMRM, 2021, 29, 2472
    input:
        wm_fod=rules.dwi2fod.output.wm_fod,
        gm_fod=rules.dwi2fod.output.gm_fod,
        csf_fod=rules.dwi2fod.output.csf_fod,
        mask=lambda wildcards: config["input_path"]["dwi_mask"][wildcards.dataset],
    output:
        wm_fod=bids(
            root=root,
            datatype="dwi",
            desc="normalized",
            suffix="wm_fod.mif",
            **config["subj_wildcards"],
        ),
        gm_fod=bids(
            root=root,
            datatype="dwi",
            desc="normalized",
            suffix="gm_fod.mif",
            **config["subj_wildcards"],
        ),
        csf_fod=bids(
            root=root,
            datatype="dwi",
            desc="normalized",
            suffix="csf_fod.mif",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mtnormalise -nthreads {threads} -mask {input.mask} {input.wm_fod} {output.wm_fod} {input.gm_fod} {output.gm_fod} {input.csf_fod} {output.csf_fod}"


rule dwi2tensor:
    input:
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),

    output:
        tensor=bids(
            root=root,
            datatype="dwi",
            suffix="tensor.mif",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "dwi2tensor {input} {output}"


rule tensor2metrics:
    input:
        tensor=rules.dwi2tensor.output.tensor,
        mask=lambda wildcards: config["input_path"]["dwi_mask"][wildcards.dataset],
    output:
        fa=bids(
            root=root,
            datatype="dwi",
            suffix="fa.mif",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tensor2metric -fa {output.fa} -mask {input.mask} {input.tensor}"


# -------------- MRTRIX PREPROC END ----------------#


# ----------- MRTRIX TRACTOGRAPHY BEGIN ------------#
rule create_seed:
    input:
        rules.tensor2metrics.output.fa,
    params:
        threshold=0.15,
    output:
        seed=bids(
            root=root,
            datatype="dwi",
            suffix="seed.mif",
            **config["subj_wildcards"],
        ),
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrthreshold {input} -abs {params.threshold} {output}"


def get_tckgen_runtime(wildcards,threads):
    streamlines=config["dwi"]["sl_count"]
    minimum_minutes=10
    minutes_per_streamline=0.0005
    
    return max(int(minutes_per_streamline * streamlines / threads),minimum_minutes)
 

rule tckgen:
    input:
        wm_fod=rules.mtnormalise.output.wm_fod,
        dwi=bids(
            root=root,
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),


        mask=lambda wildcards: config["input_path"]["dwi_mask"][wildcards.dataset],
        seed=rules.create_seed.output.seed,
    params:
        streamlines=config["dwi"]["sl_count"],
        cutoff=lambda wildcards: float(wildcards.cutoff),
        seed_strategy=lambda wildcards, input: f"-seed_image {input.seed}",
    output:
        tck=bids(
            root=root,
            datatype="dwi",
            desc="iFOD2",
            cutoff="{cutoff}",
            suffix="tractography.tck",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        runtime=get_tckgen_runtime
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tckgen -nthreads {threads} -algorithm iFOD2 -mask {input.mask} {params.seed_strategy} -cutoff {params.cutoff} -select {params.streamlines} {input.wm_fod} {output.tck}"


rule tcksift2:
    # Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. SIFT2: Enabling dense quantitative assessment of brain white matter connectivity using streamlines tractography. NeuroImage, 2015, 119, 338-351
    input:
        wm_fod=rules.mtnormalise.output.wm_fod,
        tck=rules.tckgen.output.tck,
    output:
        tckweights=bids(
            root=root,
            datatype="dwi",
            desc="sift2",
            cutoff="{cutoff}",
            suffix="tckweights.txt",
            **config["subj_wildcards"],
        ),
        mu=bids(
            root=root,
            datatype="dwi",
            desc="sift2",
            cutoff="{cutoff}",
            suffix="mu.txt",
            **config["subj_wildcards"],
        ),
    threads: 8
    resources:
        runtime=get_tckgen_runtime
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tcksift2 -nthreads {threads} -out_mu {output.mu} {input.tck} {input.wm_fod} {output.tckweights}"


rule tck2connectome:
    # Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015, 104, 253-265
    input:
        tckweights=rules.tcksift2.output.tckweights,
        tck=rules.tckgen.output.tck,
        parcellation=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        sl_assignment=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="streamassign.txt",
            **config["subj_wildcards"],
        ),
        conn=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="struc.conn.csv",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tck2connectome -nthreads {threads} -tck_weights_in {input.tckweights} -out_assignments {output.sl_assignment} -zero_diagonal -symmetric {input.tck} {input.parcellation} {output.conn}"


rule tck2connectome_fa:
    input:
        tck=rules.tckgen.output.tck,
        parcellation=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
        fa=bids(
            root=root,
            datatype="dwi",
            suffix="fa.mif",
            **config["subj_wildcards"],
        ),
    output:
        conn=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="strucFA.conn.csv",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    shadow: 'minimal'
    container:
        config["singularity"]["mrtrix"]
    shell:
        "tcksample -nthreads {threads} {input.tck} {input.fa} mean_FA_per_streamline.csv -stat_tck mean && "
        "tck2connectome -nthreads {threads}  -zero_diagonal -symmetric {input.tck} {input.parcellation} {output.conn} -scale_file mean_FA_per_streamline.csv -stat_edge mean "



rule vis_struct_connectome:
    """this rule doesn't depend on the cifti's from the func pipeline"""
    input:
        csv=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="{struc}.conn.csv",
            **config["subj_wildcards"],
        ),
    params:
        title='Conn {dataset}',
        vmin=lambda wildcards: 0 if wildcards.struc == 'strucFA' else 0.0000001,
        vmax=lambda wildcards: 1 if wildcards.struc == 'strucFA' else 500
    output:
        png=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="{struc}.conn.matrix.png",
            **config["subj_wildcards"],
        ),
    script:
        '../scripts/vis_struct_connectome.py'

rule connectome2tck:
    # Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015, 104, 253-265
    input:
        tck=rules.tckgen.output.tck,
        sl_assignment=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="streamassign.txt",
            **config["subj_wildcards"],
        ),
    output:
        tck_dir=directory(bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="struc.bundles",
            **config["subj_wildcards"],
        )),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mkdir -p {output} && connectome2tck -nthreads {threads} {input.tck} {input.sl_assignment} {output.tck_dir}/bundle_"


rule connectome2exemplartck:
    # Smith, R. E.; Tournier, J.-D.; Calamante, F. & Connelly, A. The effects of SIFT on the reproducibility and biological accuracy of the structural connectome. NeuroImage, 2015, 104, 253-265
    input:
        tck=rules.tckgen.output.tck,
        sl_assignment=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="streamassign.txt",
            **config["subj_wildcards"],
        ),
        parcellation=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),

    output:
        tck=bids(
            root=root,
            datatype="dwi",
            atlas="{atlas}",
            cutoff="{cutoff}",
            suffix="struc.exemplars.tck",
            **config["subj_wildcards"],
        ),
    threads: 8
    group:
        "grouped_subject"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "connectome2tck -nthreads {threads} {input.tck} {input.sl_assignment} {output.tck} -files single -exemplars {input.parcellation}"


rule mif2nii:
    input: '{root}/{{prefix}}.mif'.format(root=root)
    output: '{root}/{{prefix}}.nii.gz'.format(root=root)
    shell:
        'mrconvert {input} {output}'
