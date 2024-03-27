import pandas as pd
import matplotlib.colors as mcolors


# Function to convert RGB to hex
def rgb_to_hex(row):
    rgb = [row['R'], row['G'], row['B']]
    return mcolors.to_hex([x / 255.0 for x in rgb])


# Define the file path
file_path = snakemake.input.label_txt


# Initialize lists to store data
data = []
labels = []

#make Background label first
#data.append([0,0,0,0,0])
#labels.append('Background')

with open(file_path, 'r') as file:
    for idx, line in enumerate(file):
        if idx % 2 == 0:
            labels.append(line.strip())
        else:
            data.append([int(x) for x in line.strip().split()])



# Create a DataFrame
df = pd.DataFrame(data)

df.columns = ['index', 'R', 'G', 'B', 'A']

# Assign labels to the DataFrame
df['name'] = labels

df['color'] = df.apply(rgb_to_hex, axis=1)


df_out = df[['index','name','color']]

df_out.to_csv(snakemake.output.label_tsv,sep='\t',index=False)


