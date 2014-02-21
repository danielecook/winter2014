from pandas import pandas as pd
import os

os.chdir(os.path.dirname(__file__))

f = open('nhgri_1kg.recode.vcf')

df = pd.DataFrame()

for l in f.readlines()[1:1000]:
	if l.startswith("#"):
		pass
	else:
		l = l.split('\t')
		l = dict({'rs':l[0],'ref_1kg':l[1],'oth_1kg':l[2]}.items() +  {x.split('=')[0]:x.split('=')[1].strip('\n') for x in l[-1].split(';')}.items())
		df = df.append(pd.DataFrame([l]))

print df
df.to_csv('1kg_formatted.txt')
