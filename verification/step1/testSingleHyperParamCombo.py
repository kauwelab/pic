import os
import subprocess
import sys
import argparse

def setArgs(parser):
	parser.add_argument("-r", "--percentAboveRandom", help="The particular percent above random being testsed", required=True)
	parser.add_argument("-p", "--minPx", help="The minimum p(x) or p(y) value being tested", required=True)
	parser.add_argument("-a", "--inputCases", help="The location of the directories containing the case or interacting fastas",  required=True)
	parser.add_argument("-b", "--outputCases", help="The location of the directory for the output file for the case/interacting data. Please note that both case and control files will have the same name so make this different from your output control location", required=True)
	parser.add_argument("-c", "--inputControls", help="The location of the directories containing the control or not-known-to-interact fastas", required=True)
	parser.add_argument("-d", "--outputControls", help="The location of the directory for the output file for the control data. Please note that both case and control files will have the same name so make this different from your output case location", required=True)
	args = parser.parse_args()
	
	if not args.output.endswith("/"):
		args.output = args.output + "/"
	return args

parser = argparse.ArgumentParser()
args = setArgs(parser)
subprocess.call("./computeMiEntireCategory.sh %s %s %s %s" %(args.inputCases, args.outputCases, args.percentAboveRandom, args.minPx), shell=True)
subprocess.call("./computeMiEntireCategory.sh %s %s %s %s" %(args.inputControls, args.outputControls, args.percentAboveRandom, args.minPx), shell=True)





