import netplotbrain
import pandas as pd
import nibabel as nib
import numpy as np


df = pd.read_csv(snakemake.input.label_tsv,sep='\t')

conn_matrix = nib.load(snakemake.input.pconn).get_fdata()


from nichord.chord import plot_chord
from nichord.convert import convert_matrix

if 'networks' in df.columns:
    idx_to_label = {idx:label for idx,label in zip(df.index,df.networks)}
else:
    print(f'error: no networks defined for {snakemake.input.label_tsv}')
    


edges,edge_weights = convert_matrix(conn_matrix)

plot_chord(idx_to_label=idx_to_label,edges=edges,
           edge_weights=edge_weights,
           fp_chord=snakemake.output.png)



