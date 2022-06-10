import os
from Bio.SeqIO import parse

def parseFile(fname):
	seqs = []
	headers = []
	for record in parse(fname, "fasta"):
		seqs.append(str(record.seq))
		headers.append(str(record.id))
	return seqs, headers	

def writeOut(l1, l2, cat, group, heads):
	outfile = open(cat + group + "/" + group + ".fasta", 'w')
	appended = []
	for i in range(len(l1)):
		outfile.write(">" + heads[i])
		outfile.write(l1[i] + "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" + l2[i])
	outfile.close()
		
	
def processCat(cat):
	for dir in os.listdir(cat):
		g2 = dir.split("_")[1] +  "_" + dir.split("_")[0]
		g1 = dir
		seqs1, heads  = parseFile(cat + dir + "/" + g1)
		seqs2, heads = parseFile(cat + dir + "/" + g2)
		writeOut(seqs1, seqs2, cat, dir, heads)
		
	




processCat("known/")
processCat("unknown/")

