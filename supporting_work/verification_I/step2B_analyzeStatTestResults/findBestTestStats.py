import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import argparse



def setArgs(parser):

	parser.add_argument("-i", "--input", help="The location of the csv created with analyzeParameters.py", required=True)
	parser.add_argument("-o", "--output", help="The desired location of all the heat maps and output text file", required=False, default="")
	args = parser.parse_args()
	
	return args

parser = argparse.ArgumentParser()
args = setArgs(parser)

data = pd.read_csv(args.input)
outLoc = args.output

if not outLoc.endswith("/") and outLoc != "":
	outloc += "/"


statsOfInterest = data.columns.tolist()[3:5]

caseMeans = data['caseMean'].tolist()
controlMeans = data['controlMean'].tolist()
meanQuotient = [caseMeans[i]/controlMeans[i] for i in range(len(caseMeans))]
data['meanQuotient'] = meanQuotient


paCase = data['percentAboveCase'].tolist()
paControl = data['percentAboveControl'].tolist()
diffAboveT = [paCase[i] - paControl[i] for i in range(len(paCase))]
data["diffAboveThreshold"] = diffAboveT

ttestP = data['ttestP'].tolist()
negLogTtestP = []
for v in ttestP:
	if v == 0.0:
		negLogTtestP.append(np.nan)
	else:
		negLogTtestP.append(math.log(v, 10) * -1)
data['negLogTtestP'] = negLogTtestP

negLogThresholdP = []
thresholdP = data['thresholdP'].tolist()
for v in thresholdP:
	if v == 0.0:
		negLogThresholdP.append(np.nan)
	else:
		negLogThresholdP.append(math.log(v,10) *-1)
data['negLogThresholdP'] = negLogThresholdP

cd = data['cohensD'].tolist()
cd = [abs(x) for x in cd]
data['cohensD'] = cd

statsOfInterest += ["meanQuotient", "diffAboveThreshold", "negLogTtestP", "negLogThresholdP"]
outfile = open(outLoc + "bestSpecsForEachStat.txt", 'w')
for stat in statsOfInterest:
	
	pivot = data.pivot(index='percentAboveRandom', columns='minPx', values=stat)
	mask = pivot.isnull()
	if "Threshold" in stat or 'threshold' in stat:
		pivot2 = data.pivot(index='percentAboveRandom', columns='minPx', values='bestThreshold')
		tValues = pivot2.values.tolist()
		annot = False
		hm = sns.heatmap(pivot, annot=annot, cbar=True, mask=mask, annot_kws={"fontsize":7})
	else:
		annot=False
		hm = sns.heatmap(pivot, annot=annot, cbar=True, mask=mask, annot_kws={'fontsize':7})
	hm.set_title(stat)
	fig = hm.get_figure()
	fig.savefig(outLoc + stat + ".png")
	plt.close()

	statList = data[stat].tolist()
	maxVal = max(statList)
	bestAboveR = data["percentAboveRandom"].tolist()[statList.index(maxVal)]
	bestMinPx = data["minPx"].tolist()[statList.index(maxVal)]
	outfile.write("Best specifications for " + stat + "\n")
	outfile.write("Value: " +  maxVal + "\n")
	outfile.write("Percent above random: " +  bestAboveR + "\n")
	outfile.write("Minimum p(x)/p(y) value: " +  bestMinPx + "\n\n")

