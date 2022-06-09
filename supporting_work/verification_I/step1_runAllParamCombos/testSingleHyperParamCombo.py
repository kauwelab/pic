import os
import subprocess
import sys
import argparse

def setArgs(parser):
	parser.add_argument("-r", "--percentAboveRandom", help="The particular percent above random being testsed", required=True)
	parser.add_argument("-p", "--minPx", help="The minimum p(x) or p(y) value being tested", required=True)
	parser.add_argument("-ii", "--interactingInput", help="The location of the directories containing the interacting fastas",  required=True)
	parser.add_argument("-io", "--interactingOutput", help="The location of the directory for the output file for the interacting data. Please note that both interacting and unkown to interact files will have the same name so make this different from your output 'unknown' location", required=True)
	parser.add_argument("-ui", "--unknownInput", help="The location of the directories containing the fastas not known to interact.", required=True)
	parser.add_argument("-uo", "--unknwonOutput", help="The location of the directory for the output file for the data not known to interact. Please note that both interacting and unknown to interact files will have the same name so make this different from your output 'interacting' location", required=True)
	args = parser.parse_args()
	
	if not args.output.endswith("/"):
		args.output = args.output + "/"
	return args

parser = argparse.ArgumentParser()
args = setArgs(parser)
subprocess.call("./computeMiEntireCategory.sh %s %s %s %s" %(args.unknownInput, args.unknownOutput, args.percentAboveRandom, args.minPx), shell=True)
subprocess.call("./computeMiEntireCategory.sh %s %s %s %s" %(args.interactingInput, args.interactingOutput, args.percentAboveRandom, args.minPx), shell=True)





