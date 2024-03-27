import nibabel as nib
import pandas as pd
import numpy as np
import re

#add in parcel centroids/medoids, along with network column
df = pd.read_csv(snakemake.input.label_tsv,sep='\t')


coords = nib.load(snakemake.input.coords_pscalar).get_fdata().T
df_coords = pd.DataFrame(coords,columns=['x','y','z'])

df = pd.concat((df,df_coords),axis=1)

#get network label
network_pattern=snakemake.params.network_pattern
if network_pattern is not None:
    networks=[]
    for name in df.name:
        match = re.search(network_pattern,name)
        if match is not None:
            networks.append(match[1])
        else:
            networks.append('n/a')

    df['networks'] = networks



df.to_csv(snakemake.output.label_tsv,sep='\t',index=False,float_format='%.3f')

