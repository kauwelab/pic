import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse



def getMaxes(fname):

	infile = open(fname, 'r')
	scores = []	
	next(infile)
	for line in infile:
		if not np.isnan(float(line.strip().split('\t')[2])):
			scores.append(float(line.strip().split('\t')[2]))
	return scores

def getStats(threshold, cases, controls):
	tp = 0	
	fn = 0
	for score in cases:
		if score >= threshold:
			tp += 1
		else:
			fn += 1
		
	fp = 0
	tn = 0
	for score in controls:
		if score < threshold:
			tn += 1
		else:
			fp += 1
	precision = float(tp)/(tp+fp)
	recall = float(tp)/len(cases)
	tpr = tp/(len(cases))
	fpr = fp/(len(controls))
	return precision, recall, tpr, fpr



#RUN A SINGLE COMPARISON AND GET PRECISION v RECALL GRAPH and AUC GRAPH
def runOne(par, mpx, caseLoc, controlLoc, outLoc):
	
	if not outLoc.endswith("/"):
		outLoc += "/"
	fname = str(par) + "aboveR_" + str(mpx) + "minPx.tsv"
	cases = getMaxes(caseLoc + fname)
	controls = getMaxes(controlLoc + fname)
	allScores = cases + controls
	allScores.sort()
	precision_v_recall = []
	falsePos_v_truePos = []
	maxP = 0
	bestThreshold = 0
	corrRecall = 0
	for score in allScores:
		p,r,t,f = getStats(score, cases, controls)
		if p > maxP and r >= 0.2:
			bestThreshold = score
			maxP = p
			corrRecall = r
		precision_v_recall.append((r,p))
		falsePos_v_truePos.append((t,f))
	
	
	precision_v_recall.sort()
	precision_v_recall = precision_v_recall[20:]
	precisions = [x[1] for x in precision_v_recall]
	recalls = [x[0] for x in precision_v_recall]
	ycoord = len(cases)/(len(controls) + len(cases))
	plt.plot(recalls, precisions, label="Data")
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.suptitle("Precision vs Recall at Different Thresholds")
	plt.plot([0,1], [ycoord,ycoord], linestyle="--", label="Baseline")
	plt.legend()
	plt.savefig(outLoc + "precision_v_recall.png")
	plt.close()
	
	falsePos_v_truePos.sort()
	falsePosRates = [x[1] for x in falsePos_v_truePos]
	truePosRates = [x[0] for x in falsePos_v_truePos]
	plt.plot(falsePosRates, truePosRates, label="Data")
	plt.xlabel("False Positive Rate")
	plt.ylabel("True Positive Rate")
	plt.suptitle("ROC Curve")
	plt.plot([0,1],[0,1], linestyle="--", label="Baseline")
	plt.legend()
	auc = np.trapz(truePosRates, falsePosRates)
	plt.text(0.6, 0.1, "AUC=" + str(auc))
	plt.savefig(outLoc + "AUC.png")
	plt.close()


	outfile = open(outLoc + "precision_stats_" + fname[0:-4] + ".txt", 'w')
	outfile.write("Highest Precision: " + str(maxP) + "\n")
	outfile.write("Corresponding Recall: " + str(corrRecall) + "\n")
	outfile.write("Corresponding Threshold: " + str(bestThreshold) + "\n")




#RUN ALL POSSIBLE PARAMTER COMBINATIONS AND GET HEATMAPS OF THE AUC AND PRECISION OVER THE DIFFERENT HYPER PARAMETERS 
def runAll(caseLoc, controlLoc, outLoc):
	if outLoc is None:
		outLoc = ""
	bestAUC = 0
	aucs = []
	precision = []
	PARs = []
	MPXs = []
	maxP = 0
	bestThreshold = 0
	bestParamComboPrecision = ""
	bestParamComboAUC = ""
	corrRecall = 0
		
	for i in range(10,50):
		for j in range(0,36):
			PARs.append(i)
			MPXs.append(j) 
			if i == 34 and j == 0:
				precision.append(np.nan)
				aucs.append(np.nan)
				continue
			fname = str(i) + "aboveR_" + str(j) + "minPx.tsv"	
			cases = getMaxes(caseLoc + fname)
			controls = getMaxes(controlLoc + fname)
			localMaxP = 0
			falsePos_v_truePos = []
			allScores = cases + controls
			allScores.sort()
			for score in allScores:
				p,r,t,f = getStats(score, cases, controls)
				if p > localMaxP and r >= 0.2:
					localMaxP = p
				if p > maxP and r >0.2:
					bestThreshold = score
					maxP = p
					corrRecall = r
					bestParamComboPrecision = fname[0:-4]
				falsePos_v_truePos.append((t,f))
			precision.append(localMaxP)
			falsePos_v_truePos.sort()
			fpr = [x[1] for x in falsePos_v_truePos]		#For graphing purposes, I need the points to be in ascending order by true positive rate (t), which they're not guarateed to be in, but I need 
			tpr = [x[0] for x in falsePos_v_truePos]		# each false positive rate to stay coupled with its corresponding true positive rate and it was more steps to use a dictionary, so that's what this is about
			auc = np.trapz(tpr, fpr)
			aucs.append(auc)
			if auc > bestAUC:
				bestAUC = auc
				bestParamComboAUC = fname[0:-4]
	
	columns=["percentAboveRandom", "minPx", "Precision", "AUC"] 
	df = pd.DataFrame(list(zip(PARs, MPXs, precision, aucs)), columns=columns) 
	pivot = df.pivot(index="percentAboveRandom", columns="minPx", values="Precision") 
	mask = pivot.isnull()
	hm = sns.heatmap(pivot, annot=False, cbar=True, mask=mask)
	hm.set_title("Best Precisions")
	fig = hm.get_figure()
	
	fig.savefig(outLoc + "Precision.png")
	plt.close()

	pivot = df.pivot(index="percentAboveRandom", columns="minPx", values="AUC")
	mask = pivot.isnull()
	hm = sns.heatmap(pivot, annot=False, cbar=True, mask=mask)
	hm.set_title("Best AUCs")
	fig = hm.get_figure()
	fig.savefig(outLoc + "AUCs.png")
	plt.close()

	outfile = open(outLoc + "precision_stats_all.txt", 'w')
	outfile.write("Highest Precision: " + str(maxP) + "\n")
	outfile.write("Corresponding Recall: " + str(corrRecall) + "\n")
	outfile.write("Best Threshold: " + str(bestThreshold) + "\n")
	outfile.write("Corresponding Hyper Parameters: " + str(bestParamComboPrecision) + "\n")
	outfile.write("Best AUC: "+ str(bestAUC) + "\n")
	outfile.write("Corresponding Hyper Parameters: " + str(bestParamComboAUC) + "\n")
	outfile.close()



def setArgs(parser):
	"""
	Sets up the argparse
	
	Args:
		parser (ArgumentParser): the argument parser
	Return:
		args: the arguments
	"""
	parser.add_argument("-i", "--interacting", help="The directory containing all the tsv's with info about the high scores from the known-to-interact data created by testSingleHyperParamCombo.py", required=True)
	parser.add_argument("-n", "--unknown", help="The directory containing all the tsv's with info about the high scores from the not-known-to-interact data created by testSingleHyperParamCombo.py", required=True)
	parser.add_argument("-d", "--destination", help="The desired location for output files. Default current directory\n", default="", required=False, type=str)
	parser.add_argument("-a", "--all", help="Flag to compare scores for all hyper parameter pairs and create heat maps for AUC scores and best precisions at each pair\n", required=False, action="store_true")
	parser.add_argument("-o", "--one", help="Run only on a certain hyper parameter pair and create their precision vs recall and ROC graphs. List the percent above random followed by the minpx (e.g. -o 10 15)\n", required=False, nargs="+")
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	args = setArgs(parser)

	if args.all:
		runAll(args.inputCase, args.inputControl, args.destination)
	if args.one:
		runOne(args.one[0], args.one[1], args.inputCase, args.inputControl, args.destination)

	

