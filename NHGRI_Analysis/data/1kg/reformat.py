from pandas import pandas as pd
import os

os.chdir(os.path.dirname(__file__))

f = open('nhgri_1kg.recode.vcf')

df = pd.DataFrame()
i = 0
for l in f.readlines():
	if l.startswith("#"):
		pass
	else:
		l = l.split('\t')
		l = dict({'rs':l[2],'ref_1kg':l[3],'oth_1kg':l[4],'INFO':l[6]}.items() +  {x.split('=')[0]:x.split('=')[1].strip('\n') for x in l[-1].split(';')}.items())
		df = df.append(pd.DataFrame([l]))
		i+=1
		print i
print df
df.to_csv('1kg_formatted.txt')
