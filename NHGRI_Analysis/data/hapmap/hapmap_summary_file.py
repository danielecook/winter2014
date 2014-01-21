#! /usr/local/bin/Python
import os
import sqlite3
import glob
import time
import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker,query
import datetime


os.chdir(os.path.dirname(__file__))

conn = sqlite3.connect('hapmap.db')
c = conn.cursor()
f = file('../gwas_catalog_rs_list.txt','r')
o = file('hapmap_allele_freq.txt','wr')
for line in f.readlines()[0:10]:
    row_set = {}
    for row in c.execute("SELECT * FROM freq WHERE rs == '%s';" % (line.replace('\n',''))):
        row_set[row[1] + '_a1'] = row[5]
        row_set[row[1] + '_a2'] = row[8]
        row_set[row[1] + '_a1_freq'] = row[6]
        row_set[row[1] + '_a2_freq'] = row[7]
        
        o.write('\t'.join([str(x) for x in (row[2],'\n')]))


