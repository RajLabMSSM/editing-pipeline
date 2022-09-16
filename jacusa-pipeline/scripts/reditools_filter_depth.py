import numpy as np
import re, os, sys

input_file = sys.argv[1]
change_depth = int(sys.argv[2])
min_depth = int(sys.argv[3])


with open(sys.argv[1]) as infile:
	for line in infile:
		linebu = line.rstrip()
		line = line.rstrip().split("\t")
		l = np.array(line[-2].strip("[]").split(",")).astype(int)
		change = line[-1]
		depth = dict(zip(["A","C","G","T"],l))
		ref_d = int(depth[change[0]])
		alt_d = int(depth[change[1]])
		if alt_d >=change_depth and (ref_d+alt_d)>min_depth:
			print(linebu)
