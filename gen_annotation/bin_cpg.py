import pandas as pd
from operator import add
import numpy as np
import os
import sys
from collections import defaultdict
from operator import add
import numpy as np
import os
import sys
from collections import defaultdict
import math

filename = sys.argv[1]
cpg_position = pd.read_table(filename,header=None)
cpg_position.columns = ["Chr","Position"]
bin_size = 100

bin_coverage = defaultdict(lambda:0)

for index, line in cpg_position.iterrows():
	chrm_name = str(line["Chr"])
	start = int(line["Position"])
	bin_number = (start-1)//bin_size # Using integer division is the key here
	bin_label = (bin_number*bin_size) + bin_size
	bin_coverage[bin_label] += 1

i=0
for k,v in bin_coverage.items():
	i=i+1
	if (i == 50):
		break
	print ("%s - %s" % (str(k), str(v)))
