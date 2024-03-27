import netplotbrain
import pandas as pd
import nibabel as nib
import numpy as np


df = pd.read_csv(snakemake.input.label_tsv,sep='\t')

conn_matrix = nib.load(snakemake.input.pconn).get_fdata()

kwargs_opts = snakemake.params.opts

if 'networks' in df.columns:
    kwargs_network={'node_color':'networks'}
else:
    kwargs_network={}

netplotbrain.plot(nodes=df,
                  edges=conn_matrix,
                  template=snakemake.params.template,
                  savename=snakemake.output.png,
                  **kwargs_opts,
                  **kwargs_network)

