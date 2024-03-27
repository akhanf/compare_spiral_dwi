
rule plot_pconn_png:
    """generic rule for plotting pconn cifti files"""
    input:
        cifti_pconn="{prefix}.pconn.nii",
    output:
        png="{prefix}.pconn.matrix.png",
    script:
        "../scripts/plot_pconn_png.py"


rule plot_sfc_markers_png:
    """plot sfc using markers on glass brain"""
    input:
        pscalar_sfc=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="32k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="sfc.pscalar.nii",
            **config["subj_wildcards"]
        ),
        pscalar_markers="resources/atlas/atlas-{atlas}_surf-{surf}_markers.pscalar.nii",
        surfs=lambda wildcards: expand(
            config["template_surf"], surf=wildcards.surf, hemi=["L", "R"]
        ),
    output:
        png=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="32k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            surf="{surf}",
            suffix="sfc.markerplot.png",
            **config["subj_wildcards"]
        ),
    script:
        "../scripts/plot_sfc_markers_png.py"


rule plot_func_network:
    input:
        pconn=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="32k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.nii",
            **config["subj_wildcards"]
        ),
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    params:
        template="MNI152NLin6Asym",  #closest to fsLR surfs
        opts=config["netplotbrain"]["func"],
    output:
        png=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="32k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.network.png",
            **config["subj_wildcards"]
        ),
    script:
        "../scripts/plot_network.py"


rule plot_struc_network:
    input:
        pconn=bids(
            root=root,
            datatype="dwi",
            den="32k",
            atlas="{atlas}",
            suffix="struc.pconn.nii",
            **config["subj_wildcards"],
        ),
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    params:
        template="MNI152NLin6Asym",  #closest to fsLR surfs
        opts=config["netplotbrain"]["struc"],
    output:
        png=bids(
            root=root,
            datatype="dwi",
            den="32k",
            atlas="{atlas}",
            suffix="struc.pconn.network.png",
            **config["subj_wildcards"],
        ),
    script:
        "../scripts/plot_network.py"



rule plot_func_chord:
    input:
        pconn=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="32k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.nii",
            **config["subj_wildcards"]
        ),
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    params:
        opts=config['nichord']['func']
    output:
        png=bids(
            root=root,
            datatype="func",
            desc="preproc",
            den="32k",
            task="{task}",
            denoise="{denoise}",
            fwhm="{fwhm}",
            atlas="{atlas}",
            suffix="bold.pconn.chord.png",
            **config["subj_wildcards"]
        ),
    script:
        "../scripts/plot_chord.py"



rule plot_struc_chord:
    input:
        pconn=bids(
            root=root,
            datatype="dwi",
            den="32k",
            atlas="{atlas}",
            suffix="struc.pconn.nii",
            **config["subj_wildcards"],
        ),
        label_tsv="resources/atlas/atlas-{atlas}_dseg.tsv",
    params:
        opts=config['nichord']['struc']
    output:
        png=bids(
            root=root,
            datatype="dwi",
            den="32k",
            atlas="{atlas}",
            suffix="struc.pconn.chord.png",
            **config["subj_wildcards"],
        ),
    script:
        "../scripts/plot_chord.py"


