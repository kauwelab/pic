import argparse
from scipy.stats import ttest_ind, iqr
from statistics import mean, stdev, median
from math import log
from scipy.stats import chi2_contingency
from math import sqrt


def getThresholdP(threshold, case, control):
	
	caseAbove = 0
	caseBelow = 0
	for score in case:
		if score > threshold:
			caseAbove += 1
		else:
			caseBelow += 1
	
	
	controlAbove = 0
	controlBelow = 0
	for score in control:
		if score > threshold:
			controlAbove += 1
		else:
			controlBelow += 1
	
	
	percentAboveCase = caseAbove/len(case)
	percentAboveControl = controlAbove/len(control)
	data = [[controlAbove, controlBelow],[caseAbove, caseBelow]]
	if controlAbove == 0 or controlBelow == 0 or caseAbove == 0 or caseBelow == 0:
		return(False)
	if percentAboveCase < percentAboveControl:
		return(False)
	
	results = chi2_contingency(data)
	return(results[1], percentAboveCase, percentAboveControl)

def bestThreshold(case, control):
	
	maxVal = max(case)
	thresholds = [i * maxVal/300 for i in range(0,300)]
	pvalues = []
	maxVal = 0
	maxIndex = 0
	currentPcAboveHigh = 0
	currentPcAboveLow = 0
	for t in thresholds:
		results = getThresholdP(t, case, control )
		if not results:
			pvalues.append(0)
			continue
		value = -1 * log(results[0], 10)
		pvalues.append(value)
		if value > maxVal:
			maxVal = value
			maxIndex = t
			currentPcAboveHigh = results[1]
			currentPcAboveLow = results[2]
	
	return  maxIndex, 10**(-1*maxVal), currentPcAboveHigh, currentPcAboveLow


def removeOutliers(data):
	IQR = iqr(data)
	med = median(data)
	upper = 2*IQR + med
	lower = med - 2*IQR
	
	noOutliers = [x for x in data if (x < upper and x > lower)]
	return noOutliers

def ttest(case, control):

	p = ttest_ind(case, control).pvalue

	cohensD = (mean(control) - mean(case)) / (sqrt((stdev(control) ** 2 + stdev(case) ** 2) / 2))	

	return p, cohensD

def getMaxes(fName):
	infile = open(fName, 'r')
	maxes = []
	next(infile)
	for line in infile:
		line = line.strip().split("\t")
		if line[2] == 'nan':
			continue
		maxes.append(float(line[2]))
	infile.close()
	return maxes



def analyze(caseFile, controlFile):
	caseMaxes = getMaxes(caseFile)
	controlMaxes = getMaxes(controlFile)
	
	caseMaxes = removeOutliers(caseMaxes)
	controlMaxes = removeOutliers(controlMaxes)
	
	ttestP, cohensD = ttest(caseMaxes, controlMaxes)	

	bestT, bestP, percentAboveCase, percentAboveControl = bestThreshold(caseMaxes, controlMaxes)	
	
	meanCases = mean(caseMaxes)
	meanControls = mean(controlMaxes)

	return [ttestP, cohensD, bestT, bestP, percentAboveCase, percentAboveControl, meanCases, meanControls]


def setArgs(parser):

	parser.add_argument("-i", "--case", help="The location of all the case or known-to-interact tsv's", required=True)
	parser.add_argument("-n", "--control", help="The location of all the control or not-known-to-interact tsv's", required=True)
	parser.add_argument("-o", "--output", help="The desired location and name of the output csv", required=True)
	args = parser.parse_args()
	
	return args

parser = argparse.ArgumentParser()
args = setArgs(parser)

caseDir = args.case
controlDir = args.control
if not caseDir.endswith("/"):
	caseDir += "/"
if not controlDir.endswith("/"):
	controlDir += "/"

results = []
aboveRandom = [i for i in range(10, 50)]
minPx = [j for j in range(0,36)]
for i in aboveRandom:
	for j in minPx:
		fileName = str(i) + "aboveR_" + str(j) + "minPx.tsv"
		results.append([i,j] + analyze(caseDir+fileName, controlDir+fileName))


outfile = open(args.output, 'w')
outfile.write(",".join(['percentAboveRandom', 'minPx', 'ttestP', 'cohensD', 'bestThreshold', 'thresholdP', 'percentAboveCase', 'percentAboveControl', 'caseMean', 'controlMean']))
outfile.write("\n")
for r in results:
	r = [str(num) for num in r]
	outfile.write(",".join(r) + "\n")

outfile.close()
		
