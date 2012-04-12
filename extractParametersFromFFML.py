#!/usr/bin/python


import os
import math
import sys
import xmlhandler

atoms = {}

ff=sys.argv[1]
out=sys.argv[2]

c = xmlhandler.ConvertXmlToDict(ff)

t = c["file"]["topology"]
for i in range(0,len(t)):
	if t[i]["name"]=="Mass":
		break
	
b = t[i]["class"]["interaction"]
for i in range(0,len(b)):
	atom = b[i]["atom"].upper()
	mass = b[i]["parameter"]
	atoms[atom]={}
	atoms[atom]["mass"]=mass



for i in range(0,len(t)):
	if t[i]["name"]=="TypeTwoVDWInteraction":
		break
	
b = t[i]["class"]["interaction"]
for i in range(0,len(b)):
	atom = b[i]["atom"].upper()
	epsilon = b[i]["parameter"][0]
	sigma = b[i]["parameter"][1]
	
	atoms[atom]["epsilon"]=epsilon
	atoms[atom]["sigma"]=sigma
	
o = open(out,"w")
	
for key in atoms.iterkeys():
	try:
		o.write(key+" "+atoms[key]["mass"]+" "+atoms[key]["epsilon"]+" "+atoms[key]["sigma"]+"\n")
	except:
		print(key+" will be ignored")

o.close()
