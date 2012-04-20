#!/usr/bin/python


import os
import math
import sys



inname = sys.argv[1]
outname = sys.argv[2]


f = open(inname,"r")
o = open(outname,"w")

line = f.readline()
while len(line) != 0:
	c=line.split()
	for j in range(0,len(c)):
		o.write(chr(int(c[j])))
		
	line = f.readline()

f.close()
o.close()