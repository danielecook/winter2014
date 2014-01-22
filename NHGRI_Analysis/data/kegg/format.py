#! /usr/local/bin/Python
import pandas as pd
import os

os.chdir(os.path.dirname(__file__))

df = pd.read_csv('entrez_kegg.tmp',sep='\t')

# Filter by tabulation - Only interested in terms that are more widely used (3+)
#df = df.groupby('kegg_ID').filter(lambda x: len(x) > 50)
print df
df = df.pivot(index='Gene_ID', columns='kegg',values="pathway_desc")
print len(df.columns.tolist())
print df
df.to_csv('kegg_reshaped.txt',sep='\t')