#! /usr/local/bin/Python
import os
import sqlite3
import glob
import time
import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker,query
import datetime
import pandas as pd

os.chdir(os.path.dirname(__file__))


engine = create_engine('sqlite:///hapmap.db')
conn = engine.connect()
metadata = MetaData()

freq = Table('freq', metadata, autoload=True, autoload_with=engine)
"""
f = file('../gwas_catalog_rs_list.txt','r')
o = file('hapmap_allele_freq.txt','wr')


# Write the header line
o.write ('\t'.join([c.name for c in freq.columns]) + '\n')

for line in f.readlines():
	print line
	row_set = {}
	s = select([freq]).where(freq.c.rs == "%s" % (line.replace('\n','')))
	for row in conn.execute(s):
		o.write('\t'.join(map(str,row)) + '\n')
"""
# Reshape in Pandas (Pivot)
df = pd.read_csv('hapmap_allele_freq.txt', sep='\t')
df = df.drop_duplicates() # Not sure why - but a handful of exact line duplicates are appearing.
# Drop unneeded columns
print df
del df['chrom']
del df['pos']
del df['id']

df = df.pivot(index='rs', columns='population')
df.to_csv('hapmap_allele_freq_reshaped.csv')
