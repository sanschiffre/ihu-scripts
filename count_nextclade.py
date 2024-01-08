import os
import math 
import sys
import re

SNP = {}
DEL = {}
INS = {}



cols = {}
seqdata = {}
with open( "test.tsv" ) as f :
	for l in f :
		l = l.strip("\n") 
		if l == "" : continue

		tok = l.split("\t")
		if "seqName" in tok :
			cols["name"] = tok.index("seqName")
			cols["snps"] = tok.index("substitutions")
			cols["dels"] = tok.index("deletions")
			cols["ins"]  = tok.index("insertions")						
			cols["missing"] = tok.index("missing")
			continue
		elif cols == {} :
			sys.exit("Header missing")

		name = tok[cols["name"]]

		seqdata[name] = {"snp":tok[cols["snps"]].split(","), "del":tok[cols["dels"]].split(","), "ins":tok[cols["ins"]].split(","), "gap":[] }
		for g in tok[cols["missing"]].split(",") :
			if g != "" :
				for x in range(int(g.split("-")[0]), int(g.split("-")[-1])+1 ) :

					seqdata[name]["gap"].append(x)


		for s in seqdata[name]["snp"] :
			if s != "" :
				SNP[s] = int(re.search(r'\d+', s).group())

		for d in seqdata[name]["del"] :
			if d != "" :
				DEL[d] = [x for x in range(int(d.split("-")[0]), int(d.split("-")[-1])+1 )]
			
		for i in seqdata[name]["ins"] :
			if i != "" :
				INS[i] = int(i.split(":")[0])


r = []
for y in seqdata :
	r.append(y)
print ( "\t".join(["","","Present", "Absent", "NSP",""]+r) )	


for x in sorted(SNP,key=lambda x:SNP[x]) :
	r = []
	nbOK = 0
	nbMissing = 0
	nbNON = 0
	
	for y in seqdata :
		ing = 0
		if SNP[x] in seqdata[y]["gap"] :
			ing = 1
	
		if x in seqdata[y]["snp"] :
			nbOK += 1
			r.append("1")
		elif  ing :
			nbMissing += 1
			r.append("")
		else :
			nbNON += 1
			r.append("0")
	print ( "\t".join([x,"", str(nbOK), str(nbNON), str(nbMissing),""]+r) )
	
	


for x in sorted(DEL,key=lambda x:DEL[x][0]) :
	r = []
	nbOK = 0
	nbNON = 0
	nbMissing = 0

	for y in seqdata :
		ing = 0
		for z in DEL[x] :
			if z in seqdata[y]["gap"] :
				ing = 1
	
		if x in seqdata[y]["del"] :
			nbOK += 1
			r.append("1")
		elif  ing :
			nbMissing += 1
			r.append("")
		else :
			nbNON += 1
			r.append("0")
	print ( "\t".join([x,"", str(nbOK), str(nbNON), str(nbMissing),""]+r) )
	
	

for x in sorted(INS,key=lambda x:INS[x]) :
	r = []
	nbOK = 0
	nbNON = 0
	nbMissing = 0

	for y in seqdata :
		ing = 0
		for z in [INS[x], INS[x]+1 ]:
			if z in seqdata[y]["gap"] :
				ing = 1
	
		if x in seqdata[y]["ins"] :
			nbOK += 1
			r.append("1")
		elif  ing :
			nbMissing += 1
			r.append("")
		else :
			nbNON += 1
			r.append("0")
	print ( "\t".join([x,"", str(nbOK), str(nbNON), str(nbMissing),""]+r) )
	
