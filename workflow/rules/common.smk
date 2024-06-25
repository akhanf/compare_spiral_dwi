def get_dwi_conn_csv_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    atlas="{atlas}",
                    suffix="struc.conn.matrix.png",
                    **config["subj_wildcards"],
                ),
                subject=subjects,
                dataset=dataset,
                atlas=config["atlases"],
            )
        )
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    space='dwi',
                    suffix="T1w.nii.gz",
                    **config["subj_wildcards"],
                ),
                subject=subjects,
                dataset=dataset,
            )
        )
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    suffix="fa.nii.gz",
                    **config["subj_wildcards"],
                ),
                subject=subjects,
                dataset=dataset,
            )
        )


    return targets




def get_dwi_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="dwi",
                    den="32k",
                    atlas="{atlas}",
                    suffix="struc.pconn.{plottype}.png",
                    **config["subj_wildcards"],
                ),
                subject=subjects,
                dataset=dataset,
                plottype=["matrix", "chord"],
                atlas=config["atlases"],
            )
        )
    return targets


def get_func_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="func",
                    desc="preproc",
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="bold.pconn.{plottype}.png",
                    **config["subj_wildcards"],
                ),
                subject=subjects,
                dataset=dataset,
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                fwhm=config["func"]["fwhm"],
                atlas=config["atlases"],
                plottype=["matrix", "chord"],
            )
        )
    return targets


def get_sfc_targets():
    targets = []
    for dataset in config["datasets"]:
        targets.extend(
            expand(
                bids(
                    root=root,
                    datatype="func",
                    desc="preproc",
                    den="32k",
                    task="{task}",
                    denoise="{denoise}",
                    fwhm="{fwhm}",
                    atlas="{atlas}",
                    suffix="sfc.pscalar.nii",
                    **config["subj_wildcards"],
                ),
                subject=subjects,
                dataset=dataset,
                task=config["func"]["task"],
                denoise=config["func"]["denoise"].keys(),
                fwhm=config["func"]["fwhm"],
                atlas=config["atlases"],
            )
        )

    return targets

