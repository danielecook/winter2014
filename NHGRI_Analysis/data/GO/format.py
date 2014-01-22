#! /usr/local/bin/Python
import pandas as pd
from pandas.core.reshape import melt
import os

os.chdir(os.path.dirname(__file__))

df = pd.read_csv('go_annotations.tmp',sep='\t')
print df
df = df.drop_duplicates()
# Filter by tabulation - Only interested in terms that are more widely used (3+)
print df
print len(df.columns.tolist())
df = df.groupby('GO_ID').filter(lambda x: len(x) > 50)
print df
print len(df.columns.tolist())
df = df.pivot(index='Gene_ID', columns='GO_ID',values="GO_term")
print len(df.columns.tolist())
df.to_csv('GO_reshaped.txt',sep='\t')