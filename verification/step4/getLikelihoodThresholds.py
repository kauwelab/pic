import math
import numpy as np
import sys
import argparse


def setArgs(parser):

	parser.add_argument("-i", "--case", help="The tsv of the high scores from the case group with the best hyper parameter combination as determined/chosen from steps 2 and 3 ", required=True)
	parser.add_argument("-n", "--control", help="The tsv of the high scores from the control group with the best hyper parameter combination as determined/chosen from steps 2 and 3", required=True)
	parser.add_argument("-o", "--output", help="The desired location and name of the text file that will contain the python-formatted dictionary of thresholds to likelihood of interaction", required=True)
	args = parser.parse_args()
	
	return args

def getScores(filename):
	infile = open(filename, 'r')
	next(infile)
	scores = []
	for line in infile:
		line = line.strip().split("\t")
		scores.append(float(line[2]))
	infile.close()

	return scores

parser = argparse.ArgumentParser()
args = setArgs(parser)
case = args.case
control = args.control

caseScores = getScores(case)
controlScores = getScores(control)

caseScores = [ x for x in caseScores if not np.isnan(x)]
controlScores = [y for y in controlScores if not np.isnan(y)]

allScores = caseScores + controlScores

allScores.sort(reverse=True)
caseScores.sort(reverse=True)
controlScores.sort(reverse=True)

randomLikelihood = float(len(caseScores))/len(allScores)
scalar = float(len(caseScores))/len(controlScores)

thresholds = [caseScores[int(len(caseScores)/20) * i] for i in range(1,21)]
lastControlIndex = 0
thresholdToProb = {}
for i in range(len(thresholds)):
	if i == 0:
		truePos = len(caseScores) - caseScores.index(thresholds[i])
	else:
		truePos = caseScores.index(thresholds[i]) - caseScores.index(thresholds[i-1])
	falsePos = 0
	currentControl = controlScores[lastControlIndex]
	while currentControl >= thresholds[i]:
		lastControlIndex += 1
		falsePos += 1
		currentControl = controlScores[lastControlIndex]
	denom = (scalar*falsePos + truePos)
	likelihood = truePos/(scalar*falsePos + truePos)
	if likelihood < randomLikelihood:
		break
	thresholdToProb[thresholds[i]] = likelihood

outfile = open(args.output, 'w')
outfile.write(str(thresholdToProb))
outfile.close()
	
	

