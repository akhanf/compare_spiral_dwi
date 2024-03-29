import pandas as pd

df_subjects = pd.read_csv(snakemake.input.subjects_tsv,sep='\t',dtype={"participant_label": str})

df_demog = pd.read_table(snakemake.input.ndar_txt,header=0,skiprows=[1])

df = pd.merge(left=df_subjects,right=df_demog,left_on='participant_label',right_on='src_subject_id')


df['age'] = df['interview_age'] / 12.0 #months to years

df[['participant_label','age','sex']].to_csv(snakemake.output.tsv,sep='\t',index=False)
