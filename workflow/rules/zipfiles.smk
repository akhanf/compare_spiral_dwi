wildcard_constraints:
    subject="[a-zA-Z0-9]+",
    app="|".join(config['in_agg_zip'].keys()),
    dataset="[a-zA-Z0-9]+"


def get_zip_file(wildcards):
    if wildcards.app in config["in_agg_zip"]:
        return os.path.join(in_root,config["in_agg_zip"][wildcards.app][wildcards.dataset])
    else:
        print(
            f"ERROR: cannot find zip file for {wildcards.app}/{wildcards.dataset} in config['in_agg_zip']"
        )
        return None


rule get_from_zip:
    """ This is a generic rule to make any file within the {app} subfolder, 
        by unzipping it from a corresponding zip file"""
    input:
        zip=get_zip_file,
    output:
        "{app}/{dataset}/{file}",  # you could add temp() around this to extract on the fly and not store it
    shell:
        "unzip -d {wildcards.app} {input.zip} {wildcards.dataset}/{wildcards.file}"


rule process_file:
    """ fake example of how you can make use of files that live inside the zipfile"""
    input:
        pial="snakebatch/LOBE/derivatives/fmriprep/sub-{subject}/anat/sub-{subject}_acq-MPRvNavAvgEcho12_run-1_hemi-{hemi}_pial.surf.gii",
    output:
        postproc=bids(subject="{subject}", hemi="{hemi}", suffix="pial.surf.gii"),
    shell:
        "cp {input} {output}"
