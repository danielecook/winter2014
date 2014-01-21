#! /usr/local/bin/Python
import sqlite3
import os
import glob
import time
import sqlalchemy
from sqlalchemy import Table, Column, Index, Integer, String, Float, MetaData, ForeignKey
from sqlalchemy import create_engine
import datetime


os.chdir(os.path.dirname(__file__))

if os.path.isfile('hapmap.db'):
	os.remove('hapmap.db')

engine = create_engine('sqlite:///hapmap.db')
conn = engine.connect()

metadata = MetaData()

freq = Table('freq', metadata,
	Column('id', Integer, primary_key=True),
    Column('population', String(3)),
    Column('rs', Integer),
    Column('chrom', String(5)),
    Column('pos', Integer),
    #Column('strand',String(1)), # always '+''
    #Column('build',String(100)),
    #Column('center',String(100)),
    #Column('protLSID',String(100)),
    #Column('assayLSID',String(100)),
    #Column('panelLSID',String(100)),
    #Column('QC_code',String(3)),
    #
    Column('refallele',String(3)),
    Column('refallele_freq',Float),
    Column('refallele_count',Integer),
	#
    Column('otherallele',String(3)),
    Column('otherallele_freq',Float),
    Column('otherallele_count',Integer),
    #
    Column('totalcount',Integer),
    sqlite_autoincrement=True,
)



metadata.create_all(engine)

for allele_file in glob.glob("allele*"):
	f = file(allele_file,'r')
	print f
	print datetime.datetime.now()
	pop = allele_file[allele_file.find('_',allele_file.find('chr')+1)+1:allele_file.find('_',allele_file.find('chr')+1)+4]
	h = f.readline().replace('#','').replace('\n','')
	inserts = []
	c = 0
	for line in f.readlines():
		k = dict(zip(h.split(' '), line.split(' ')))
		k['population'] = pop
		c += 1
		inserts.append(k)
		if c == 1000:
			conn.execute(freq.insert(),inserts)
			inserts = []
			c = 0
	conn.execute(freq.insert(),inserts)

# Add indices
Index('population', freq.c.population).create(engine)
Index('rs', freq.c.rs).create(engine)
Index('chrom', freq.c.chrom).create(engine)
Index('pos', freq.c.pos).create(engine)
