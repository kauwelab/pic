import pandas as pd
import os
from scipy.stats import ttest_ind

def evaluateGroup(group):

	allFastas = os.listdir("output/"+group + "_results")
	count = 0	
	scores1 = []
	scores2 = []
	for fasta in allFastas:
		loc = group + "_results/" + fasta + "/"
		try:
			df = pd.read_csv(loc + "PF00000-NONREDUNDANT-5DD-dim0-table.txt", sep="\t")
		except:
			continue
		positions = df["# list_of_positions"]
		if len(positions) == 0:
			scores1.append(0)
			scores2.append(0)
			continue
		score1 = df["Ssymm"]
		score2 = df["Senv"]
		infile = open(group + "_concat/" + fasta + ".fasta", 'r')
		next(infile)
		seq = next(infile)
		divideStart = seq.find("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX") 
		divideEnd = divideStart + 50
		gene1 = False
		gene2 = False
		ixScores1 = []
		ixScores2 = []
		for row in range(len(positions)):
			interactors = positions[row].split("/")
			for pos in interactors:
				if int(pos) < divideStart:
					gene1 = True
				if int(pos) > divideEnd:
					gene2 = True
			if gene1 and gene2:
				count += 1
				break
				ixScores1.append(score1[row])
				ixScores2.append(score2[row])
		if len(ixScores1) == 0:
			scores1.append(0)
			scores2.append(0)
			continue
		scores1.append(max(ixScores1))
		scores2.append(max(ixScores2))
	return(count)


results = evaluateGroup("unknown")
print("Interactions found in unknown-to-interact pairs:" ,results)
results = evaluateGroup("known")
print("Interactions found in known-to-interact pais:", results

