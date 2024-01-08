import os
import math 
import sys
import re

SNP = {}
DEL = {}
INS = {}

all_mutation = {"snp":{}, "del":{},"ins":{} }



countAA= 0
cols = {}
seqdata = {}
with open( "test.tsv" ) as f :
	for l in f :
		l = l.strip("\n") 
		if l == "" : continue

		tok = l.split("\t")
		if "seqName" in tok and countAA :
			cols["name"] = tok.index("seqName")
			cols["snps"] = tok.index("aaSubstitutions")
			cols["dels"] = tok.index("aaDeletions")
			cols["ins"]  = tok.index("aaInsertions")						
			continue
		elif "seqName" in tok :
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
		if name == "21K ( BA.1 )"  :
			print ( seqdata[name] )
			
			
		if countAA :
			for s in seqdata[name]["snp"] :
				if s != "" :
					all_mutation["snp"][s] = [ (s.split(":")[0], int(re.search(r'\d+', s.split(":")[1]).group())) ]

			for d in seqdata[name]["del"] :
				if d != "" :
					all_mutation["del"][d] = [ (d.split(":")[0], int(re.search(r'\d+', d.split(":")[1]).group())) ]
				
			for i in seqdata[name]["ins"] :
				if i != "" :
					all_mutation["ins"][i] = [ (i.split(":")[0], int(re.search(r'\d+', i.split(":")[1]).group())) ]
		
		
		else :
			for g in tok[cols["missing"]].split(",") :
				if g != "" :
					for x in range(int(g.split("-")[0]), int(g.split("-")[-1])+1 ) :

						seqdata[name]["gap"].append(x)

			for s in seqdata[name]["snp"] :
				if s != "" :
					all_mutation["snp"][s] = [ int(re.search(r'\d+', s).group()) ]

			for d in seqdata[name]["del"] :
				if d != "" :
					all_mutation["del"][d] = [ x for x in range(int(d.split("-")[0]), int(d.split("-")[-1])+1 ) ]
				
			for i in seqdata[name]["ins"] :
				if i != "" :
					all_mutation["ins"][i] = [ int(i.split(":")[0]) ]


r = []
for y in seqdata :
	r.append(y)
print ( "\t".join(["","","Present", "Absent", "NSP",""]+r) )	


for w in ["snp","del","ins"] :
	dic = {k: all_mutation[w][k] for k in sorted(all_mutation[w], key=lambda x:all_mutation[w][x] ) } 
	for x,pos in dic.items() :
		r = []
		nbOK = 0
		nbMissing = 0
		nbNON = 0
		
		for y in seqdata :
			ing = 0
							
			for z in pos :
				if z in seqdata[y]["gap"] :
					ing = 1
		
			if x in seqdata[y][w] :
				nbOK += 1
				r.append("1")
			elif  ing :
				nbMissing += 1
				r.append("")
			else :
				nbNON += 1
				r.append("0")
				
		print ( "\t".join([x,"", str(nbOK), str(nbNON), str(nbMissing),""]+r) )
		

