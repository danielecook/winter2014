#! /usr/local/bin/Python
import pandas as pd
from pandas.core.reshape import melt
import os

os.chdir(os.path.dirname(__file__))

df = pd.read_csv('gene2go.human.txt',sep='\t')

df = df.drop_duplicates()


# Filter by tabulation - Only interested in terms that are more widely used (3+)

print df
df = df.groupby('GO_ID')
print df
df = df.pivot(index='GeneID', columns='GO_ID',values="GO_term")
#df.to_csv('GO_reshaped.csv')