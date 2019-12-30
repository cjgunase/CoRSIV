import sys
import os

## Open the file with read only permit
f = open(sys.argv[1])
## If the file is not empty keep reading line one at a time
## till the file is empty
for line in iter(f):
	cmd = "python ../full_annotation.py CpG_in_"+line.split("\t")[0]+".csv "+line.split("\t")[0]+".gene "+line.split("\t")[0]+" "+line.split("\t")[1]
	print(cmd)
	os.system(cmd)
f.close()
