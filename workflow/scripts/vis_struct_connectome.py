import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


mpl.use('Agg')


in_csv = snakemake.input.csv
title = snakemake.params.title
out_png = snakemake.output.png

df = pd.read_csv(in_csv)


cmap = mpl.cm.get_cmap("viridis")
cmap.set_under(color='grey')

imshow_opts={'vmin':snakemake.params.vmin,'vmax':snakemake.params.vmax,'interpolation': 'none','cmap':cmap}



fig = plt.figure(figsize=(12,8))


plt.imshow(df.to_numpy(),**imshow_opts)
plt.title(title)

plt.colorbar()
plt.savefig(out_png)
