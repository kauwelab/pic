
from multiprocessing import Pool,freeze_support
from multiprocessing.pool import ThreadPool
import re 
import traceback
import sys
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import math
import argparse
import time

filtered1 = []
filtered2 = []


class MsaException(Exception):
	message = ''
	def __init__(self):
		Exception.__init__(self)
class AminoException(MsaException):
	def __init__(self, char, geneNum):
		self.message = "Invalid character in protein sequence (" + char + ") found on species #" + str(geneNum) 
class FastaException(MsaException):
	def __init__(self):
		self.message = "Two fasta headers (lines beginning with '>') in a row. Please enter correctly formatted fastas." 
class LengthException(MsaException):
	def __init__(self, length1, length2, geneNum):
		self.message = "MSA formatted incorrectly. \nLength of model gene: " + str(length2) + " \nLength of gene #" + str(geneNum) + ": "  + str(length1)
class NumSpeciesException(MsaException):
	def __init__(self):
		self.message = "ERROR: The two MSAs have an unequal number of species."
class DiffSpeciesException(MsaException):
	def __init__(self):
		self.message = "ERROR: Species list from the first file does not match that of the second.\nPlease make sure two files contain exactly the same list of species"



def makeFileMatrix(filename):	
	"""
	Turns an MSA into a list of strings with newlines removed 
	
	Args:
		fileName (string): 	the name of the fasta file with the MSA
	
	Return:
		matrix (list of string): A list representing the MSA with fasta headers as the even entries and the protein 
			sequences as the odd entries 
	"""
	fullPath = "app/static/uploads/" + filename
	infile = open(fullPath, "r")
	matrix = [next(infile).strip()]
	gene = ''
	for line in infile:
		line = line.strip()
		if line[0] == '>':
			if len(gene) == 0:
				raise FastaException()
			else:
				matrix.append(gene)
				matrix.append(line)
				gene = ''
		else:
			gene += line
	matrix.append(gene)	
	infile.close()
	return matrix

def orderMultiArray(resultTuples, outfile):
        """ 
        Orders an array of tuples (containing an integer and a list) by a single integer in the tuple, 
        and then places all corresponding lists in a multi-array, removing the single integers.
        
        Args:
                resultTuples (list of tuple of int and list): a list of tuples containing the mutual information values
                        of a single position on the first gene versus all of the positions on the second gene, 
                        and a single integer indicating the specific position on gene 1 to which the information corresponds.
                outfile (writable file): The csv file that will contain the data in the multi-array for user use. 
        
        Return:
                multiArray (numpy array): an NxM numpy multi-array where N is the length of the first translated gene 
                        and M is the length of the second translated gene
        """
        multiArray = np.empty(shape = (len(resultTuples), len(resultTuples[0][1])))
        resultTuples.sort()
        for index in range(len(resultTuples)):
                multiArray[index] = np.array(resultTuples[index][1])
                for num in multiArray[index]:
                        outfile.write(str(num) + ",") 
                outfile.write('\n')
        return multiArray
			
		
def errorCheck(matrix1, matrix2):
	"""
	Checks the input files for errors
	
	Args:
		matrix1 (list of string): A list representing the lines in a file. Odd lines are protein sequences, even are headers
		matrix2 (list of string): The same as matrix1 but for the other file
	"""
	aminoAcids = "- X \" A C D E F G H I K L M N P Q R S T V W Y"

	if len(matrix1) != len(matrix2):
		raise NumSpeciesException()
		
	lenGene1 = len(matrix1[1])
	lenGene2 = len(matrix2[1])	
	species1 = []
	species2 = []
	i = 0
	while i < len(matrix1):
		species1.append(matrix1[i])
		species2.append(matrix2[i])
		i += 1
		if len(matrix1[i]) != lenGene1:
			raise LengthException(len(matrix1[i]), lenGene1, i/2 + 1)
		if len(matrix2[i]) != lenGene2:
			raise LengthException(len(matrix2[i]), lenGene2, i/2 + 1)
		for char in matrix1[i]:
			if char not in aminoAcids:
				raise AminoException(char, i/2 + 1)	
		for char in matrix2[i]:
			if char not in aminoAcids:
				raise AminoException(char, i/2 + 1)
		i += 1
				
	species1.sort()
	species2.sort()
	
	if species1 != species2:
		raise DiffSpeciesException()
	
	return	
	
	

def setArgs(parser):
	"""
	Sets up the argparse

	Args:
		parser (ArgumentParser): the argument parser
	Return:
		args: the arguments
	"""
	parser.add_argument("fasta1", help="The multiple sequence alignment of the protein sequence of the first gene")
	parser.add_argument("gene1", help="the name of the first gene as you want it to appear on the heat map")
	parser.add_argument("fasta2", help="the multiple sequence alignment of the protein sequence of the second gene")
	parser.add_argument("gene2", help="the name of the second gene as you want it to appear on the heat map")
	parser.add_argument("-o", "--output", help="Base for the output file(s). File extension(s) will be added", nargs="?", required=False, type=str)
	parser.add_argument("-c","--color", help="the color scheme of the heatmap", nargs="?", type=str, default="viridis")	
	parser.add_argument("-m", "--noheatmap", action="store_true", help="Flag to not create a heatmap", required=False)
	parser.add_argument("-t", "--threads", type=int, help="The number of threads, default all", default=40, required=False) 
	parser.add_argument("-f", "--filter", help="The minimum hamming distance allowed for a sequence to be considered part of the alignment. \n\tDefault is 0.02", nargs='?', type=float,  default=0.02, required=False)	
	args = parser.parse_args()

	colors = "viridis plasma inferno magma Greys Purples Blues Greens Oranges Reds YlOrBr YlOrRd OrRd PuRd RdPu BuPu GnBu PuBu YlGnBu PuBuGn BuGn YlGn binary gist_yarg gist_gray gray bone pink spring summer autumn winter cool Wistia hot afmhot gist_heat copper"
	if args.color not in colors: 
		print("Please enter a valid color scheme (viridis,plasma,inferno,magma,Greys,Purples,Blues,Greens,Oranges,Reds,YlORBr,YlOrRd,OrRd,PuRd,RdPu,BuPu,GNBU,PuBU,YlGnBU,PuBUGn,BuGn,YlGn,binary,gist_yarg,gist_gray,gray,bone,pink,spring,summer,autumn,winter,cool,Wistia,hot,afmhot,gist_heat,copper)")
		exit()	
	if not(args.output is None)  and "." in args.output:
		print ("Please do not include file extensions or \".\" in file names. Extensions will be added")
		exit()	
	return args
	


def plotHeat(args, multiArray):
	"""
	Creates and saves the heat map of the data
	
	Args:
		args : arguments from the command line. Contains the user's specified color scheme for the heat map
		multiArray (list of list of float): the multi-array whose values will be turned into the heat map
	"""
	color = args.color
	plt.suptitle("  Mutual Information " + args.gene1 + " vs " + args.gene2 + "   ", fontsize=16)
	plt.ylabel(args.gene1)
	plt.xlabel(args.gene2)
	plt.ylim(0, len(multiArray))
	plt.imshow(multiArray, cmap=color, interpolation="nearest")
	plt.colorbar()
	plt.savefig(args.gene1 + "_" + args.gene2 + "2.png")



def findHammDist(g1s1, g2s1, g1s2, g2s2):
	normalize = 0
	if len(g1s1) > len(g2s1):
		longer1 = g1s1
		longer2 = g1s2
		shorter1 = g2s1
		shorter2 = g2s2
	else:
		longer1 = g2s1
		longer2 = g2s2
		shorter1 = g1s1
		shorter2 = g1s2
	normalize = len(longer1)/len(shorter1)
	diffs = 0
	i = 0
	while i < len(shorter1):
		ch1 = shorter1[i]
		ch2 = shorter2[i]
		if ch1 != ch2:
			diffs += 1
		i += 1
	j = 0
	while j < len(longer1):
		ch1 = longer1[j]
		ch2 = longer2[j]
		if (ch1 != ch2):
			diffs += normalize
		j += 1
	hd = diffs/(len(shorter1) * 2)
	return hd



def filterSequences(matrix1, matrix2, threshold):
	filtered1 = []
	filtered2 = []
	filtered1.append(matrix1[1])
	filtered2.append(matrix2[1])
	i = 3
	j = 3
	while i < len(matrix1) -1:
		j = i + 2
		valid = True
		while j < len(matrix1):
			hd = findHammDist(matrix1[i], matrix2[i], matrix1[j], matrix2[j])
			if hd < threshold:
				valid = False
			j += 2
			if valid: 
				filtered1.append(matrix1[i])
				filtered2.append(matrix2[i])
			i += 2
	
	return (filtered1, filtered2)
	
				

def transposeAlignment(matrix):
	"""
	Transposes the MSA and finds and removes non-relevant positions (i.e. positions in the model organism with a '-') as well as checks the file for errors

	Args:
		matrix(list of strings): A list representing the MSA, where every even entry is the fasta header, and every odd entry is the protein sequence
		
	Return: 
		transpose(list of lists of char): A multi-array representing the MSA, but tranposed so that columns are now rows and vise versa and with non-relevant positions removed. 
	"""
	positions = []
	transpose = []	
	modelGene = matrix[0]
	for i in range(len(modelGene)):
		if modelGene[i] != '-':
			positions.append(i)	
			column = [modelGene[i]]
			transpose.append(column)

	i = 1
	while i < len(matrix):
		for j in range(len(positions)):
			transpose[j].append(matrix[i][positions[j]])
		i += 1
			
	return transpose



def breakdown(pxy, px, py):
	"""
	Breaks down a complex calculation for mutual information. I couldn't make the equation work originally as pxy(math.log((pxy/px*py),10)), so I broke it down into
	parts to see which was causing the issue, and it ran perfectly, so I haven't messed with it since. 

	Args:
		pxy (float): the frequency of a pairing of two certain positions on the two genes
		px (float): the frequency of a residue at a certain position on one of the genes
		py (float): the frequency of a residue at a certain position on the other gene. 

	Return:
		d (float): a float that is the result of pxy(log(pxy/px*py)).
	"""
	a = px*py
	b = pxy/a
	c = math.log(b,10)
	d = pxy*c
	return d



def makeDictionary(column):
	"""
	Makes a dictionary with the different protein residue types in a given position as the keys and a list of the indicies at which those residues ocur as the values. 

	Args:
		column (list of char): A list (with a length equal to the number of species in the MSA's) of all residues found at a single certain position on one of the genes. 
			Each index in the list corresponds to the species on which the residue is found. 
	
	Return:
		dicitonary (dictionary of char to list of int): A dictionary with all possible residues in a certain position as the keys and a list of the indicies from "column"
			at which those residues occur as the values. 
	"""
	dictionary = {}
	index = 0
	for residue in column[:-1]:
		if residue in dictionary:
			dictionary[residue].append(index)
		else:
			dictionary[residue] = [index]
		index += 1
	return dictionary



def calcPxy(indices1, indices2):
	"""
	Calculates the frequency of a certain pairing by looking at the lists of the indicies(i.e. the list of ennumerated species)  at which each of the two residues in the pair occur, 
	and seeing where the lists overlap. 

	Args:
		indicies1 (list of int): a list of all the indicies (with each index corresponding to a certain species in the MSA) at which a certain amino acid is present on a certain position in one of the genes
		indicies2 (list of int): the same as indicies one, but for another amino acid at another position on the other gene. 

	Returns:
		total (int): the total number of times that the pairing occurs.  
	"""
	total = 0.0
	for entry in indices1:
		if entry in indices2:
			total += 1.0 
	return total
	
	
	
def calcMutualInfo(X, Y):
	"""
	Calculates the mutual information of each position pair. px is the probability of x being found at that particular position, py is the same for y, pxy is the probability of the pairing
	The formula is Σ[for all x in X] Σ[for all y in Y] pxy(log(pxy/px*py))
	
	Args:
		X (list of char): a list (with a lenght equal to the number of species in the MSA's) of all residues found at a single certain position on one of the genes. 
			Each index in the list corresponds to the species on which the residue is found.
		Y (list of char): same as X but for the other gene

	Returns:
		mutualInfo: The mutual information between the two positions represented by X and Y. The mutual information will be a number between 0 and 1.
	"""
	Xresidues = makeDictionary(X)
	Yresidues = makeDictionary(Y)	
	px = 0.0
	py = 0.0
	pxy = 0.0
	mutualInfo = 0.0
	for residue1 in Xresidues:
		residue1positions = Xresidues[residue1]
		px = float(len(residue1positions))/float(len(X))
		for residue2 in Yresidues:
			residue2positions = Yresidues[residue2]
			py = float(len(residue2positions))/float(len(Y))
			pxy = calcPxy(residue1positions, residue2positions)/len(Y)
			if(pxy != 0):
				bd = breakdown(pxy, px, py)
				mutualInfo += bd
			
	return mutualInfo



def makeMiMultiArray(tuple):
	"""
	When paired with multi-threading, Makes an MxN matrix where M is the number of positions on gene 2 
	and N is the number of positions on gene 1. Each entry of the matrix represents the mutual information score 
	between two of the positions.
	
	Args:
		transpose1 (list of char): list of the different amino acids that appear in a single position on one of the 
			genes across all the differnt species in the alignment
		transpose2 (list of char): same as transpose1 but for the other gene
		otufile: the file that the matrix will be written to in a csv format
		
	Return:
			tuple (tuple of int and list of float): a tuple containing the list of calculated MI values for a single
			position on the first gene when compared with every position on the second gene. The int is an indication
			of which position on the first gene. 
	"""
	x = tuple[2] 
	transpose2 = tuple[1]
	rowNum = tuple[0]
	size = len(transpose2)	
	row = np.empty(size)
	for y in range(size):	
		row[y] = calcMutualInfo(x, transpose2[y])
	tuple = (rowNum, row)
	return tuple	


if __name__ == '__main__':

	start_time = time.time()	
	freeze_support()
	print("freeze support: %s seconds" %(time.time() - start_time))	
	parser_start = time.time()	
	parser = argparse.ArgumentParser()
	args = setArgs(parser)
	print("argParse: %s seconds" %(time.time() - parser_start))

	matrix_start = time.time()
	matrix1 = makeFileMatrix(args.fasta1)
	matrix2 = makeFileMatrix(args.fasta2)
	print("make matrix: %s seconds" %(time.time() - matrix_start))

	error_start = time.time()
	try:
		errorCheck(matrix1, matrix2)
	except MsaException as m:
		traceback.print_exc()
		print(m.message)
		exit()
	print("check errors: %s seconds" %(time.time() - error_start))
	
	filter_start = time.time()
	threshold = args.filter	
	tuple = filterSequences(matrix1, matrix2, threshold)
	filtered1 = tuple[0]
	filtered2 = tuple[1]
	print("filter sequences: %s seconds" %(time.time() - filter_start))
	
	transpose_start = time.time()
	transpose1 =  transposeAlignment(filtered1)
	transpose2 =  transposeAlignment(filtered2)
	print("transpose alignments: %s seconds" %(time.time() - transpose_start))

	threading_start = time.time()

	outfile = open(args.gene1+ "_" + args.gene2 + ".csv", "w")
	p = Pool(args.threads)
	argsTuples = [(i, transpose2,  transpose1[i]) for i in range(len(transpose1))]
	multiArray = orderMultiArray(p.map( makeMiMultiArray, argsTuples))
	outfile.close()
	
	print("threading: %s seconds" % (time.time() - threading_start))
	
	plot_start = time.time()
	if(not args.noheatmap):
		plotHeat(args, multiArray)
	print("plot heat map: %s seconds" % (time.time() - plot_start))
	print("total time: %s seconds" %(time.time() - start_time))



