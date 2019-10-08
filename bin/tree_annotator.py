import pandas as pd  
import os
import sys
import shutil

df = pd.read_csv(sys.argv[1],sep='\t')
tf = sys.argv[2]
tf_out = sys.argv[3]

shutil.copy(tf,tf+".backup")

for index, row in df.iterrows():
	orig = row[df.columns[0]].strip()
	replace = row[df.columns[1]].strip()
	#print orig + " -- " + replace
	sedcmd = "sed -i \'s/%s/%s/g\' %s" %(orig,replace+"("+orig+")",tf)
	#print sedcmd
	os.system(sedcmd)

shutil.move(tf,tf_out)
shutil.move(tf+".backup",tf)
