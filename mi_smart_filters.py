from multiprocessing import Pool,freeze_support
from multiprocessing.pool import ThreadPool
import traceback
import sys
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import math
import argparse
import statistics
import time
import os



class MsaException(Exception):
	message = ''
	def __init__(self):
		Exception.__init__(self)
class AminoException(MsaException):
	def __init__(self, char, spec):
		self.message = "Invalid character in protein sequence (" + char + ") found on species " + spec 
class FastaException(MsaException):
	def __init__(self):
		self.message = "Two fasta headers (lines beginning with '>') in a row. Please enter correctly formatted fastas." 
class LengthException(MsaException):
	def __init__(self, length1, length2, speciesName):
		self.message = "MSA formatted incorrectly. \nLength of model gene: " + str(length2) + " \nLength of gene for " + speciesName + ": "  + str(length1)
class DiffSpeciesException(MsaException):
	def __init__(self):
		self.message = "ERROR: Species list from the first file does not match that of the second.\nPlease make sure two files contain exactly the same set of species"
class TooFewAlignmentsException(MsaException):
    def __init__(self, length):
        self.message = "ERROR: Too few sequences. For accurate results please use only alignments with 100 sequences or more. You only provided " + str(length)





class MiMatrix():
    def __init__(self, matrix):
        self.matrix = matrix
        self.distributions = {}
    

    def orderMatrix(self): 
        """
        Orders the output of multithreading by ordering by row number and then extracting the row data.
        Before this function is called, self.matrix is a list of tuples, where each tuple contains a single
        integer indicating the desired row number and a list of MI scores. The list is ordered by the row number
        and then the row is extracted and placed in a numpy multi array. 
        """
        multiArray = np.empty(shape = (len(self.matrix), len(self.matrix[0][1])))
        self.matrix.sort()
        for index in range(len(self.matrix)):
            line = self.matrix[index][1]
            multiArray[index] = np.array(line)
        self.matrix = multiArray


    def makeCsv(self, name):
        """
        Writes the mutual information score matrix to a .csv file

        Args:
            name (string): The name of the output file
        """
        outfile = open(name, 'w')
        for line in self.matrix: 
            for data in line:
                outfile.write(str(data)+",")
            outfile.write("\n") 
        outfile.close()
    

    def correctScores(self):
        """
        Performs a correction on all scores to eliminate phylogenetic bias and bias towards longer proteins
        by dividing each score by the mean score and then by the total number of scores in the score matrix.
        """
        numScores = len(self.matrix) * len(self.matrix[0])
        meanScore = statistics.mean([statistics.mean(row) for row in self.matrix])
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):
                self.matrix[i][j] = 10000*(self.matrix[i][j]/meanScore)/numScores
        
    
          

        



class MutualInfoCalculator():
    def __init__(self, file1, file2, minPx, percentAboveRandom):
        self.fasta1 = file1
        self.fasta2 = file2
        self.transpose1 = []
        self.transpose2 = []
        self.minPx = minPx/100.0
        self.interactionRequirement = 5
        self.percentOfRand = 1.0 + percentAboveRandom/100


    def makeFileDictionary(self, fullPath):
        """
        Reads through the file and creates a dictionary with headers as keys and sequences as values. Newlines are removed
        from sequence data. This is done in order to ensure that both files have the same set of headers and to order their
        sequences accordingly. Also checks to make sure there are at least 100 different species in the file. 
        Args:
            fullPath (string): The path to the file
        Returns:
            headerToGene (dictionary of string:string): The above mentioned dictionary
        """
        try: 
            infile = open(fullPath, "r")
        except FileNotFoundError:
            print(fullPath) 
            print("files missing")
            exit()
        headerToGene = {}
        key = next(infile).strip()
        self.modelOrg = key
        value = next(infile).strip()
        if value.startswith(">"):
            raise FastaException()
        for line in infile:
            line = line.strip()
            if line[0] == '>':
                if len(value) == 0:
                    raise FastaException()
                else:
                    headerToGene[key] = value
                    value = ''
                    key = line
            else:
                value += line
        headerToGene[key] = value
        if len(headerToGene) < 100:
            raise TooFewAlignmentsException(len(headerToGene))
        infile.close()
        return headerToGene


    def prepareAlignments(self):
        """
        Calls the necessary functions to read the files into a workable data structure, check the files for errors, 
        and transpose the alignments for mutual info calculations. 
        """
        try:
            self.alignment1 = self.makeFileDictionary(self.fasta1)
        except MsaException as m:
            print(m.message)
            exit()
        try:
            self.alignment2 = self.makeFileDictionary(self.fasta2)
        except MsaException as m:
            print(m.message)
            exit()
       
        try:
            self.errorCheck()
        except MsaException as m:
            print(m.message)
            exit()

        self.transpose1 = self.transpose(self.alignment1)
        self.transpose2 = self.transpose(self.alignment2)
       

    def errorCheck(self):
        """
        Checks the input files for further errors. Checks that the set of species between the two alignments
        are identical, that the sequences are fully aligned, and that there are no invalid characters in the
        sequence data. 
        """
        
        aminoAcids = "- X J Z B A C D E F G H I K L M N P Q R S T V W Y"
        
        species1 = list(self.alignment1.keys())
        species2 = list(self.alignment2.keys())
        lenGene1 = len(self.alignment1[species1[0]])
        lenGene2 = len(self.alignment2[species2[0]])	

        species1.sort()
        species2.sort()
        
        if species1 != species2:
            raise DiffSpeciesException()

        self.species = species1
        
        for key in species1:
            if len(self.alignment1[key]) != lenGene1:
                raise LengthException(len(self.alignment1[key]), lenGene1, key)
            if len(self.alignment2[key]) != lenGene2:
                raise LengthException(len(self.alignment2[key]), lenGene2, key)
            for char in self.alignment1[key]:
                if char not in aminoAcids:
                    raise AminoException(char, key)
            for char in self.alignment2[key]:
                if char not in aminoAcids:
                    raise AminoException(char, key)
    

    def transpose(self, alignment):
        """
        Transposes the alignments. Takes the sequence data, removes all positions not present in the model sequence
        (represented by a '-') and then creates a list of lists of char, with each inner list being a list of the different
        amino acids present at that position across all the species in the dataset. Each index corresponds to a specific species.
        Args:
            alignment (dict): The header:sequence dictionary that needs to be converted to a transposed list of lists
        Returns:
            transpose (list of list of char): The transposed alignment mentioned above
        """
        positions = []
        transpose = []
        modelSequence = alignment[self.modelOrg]
        for i in range(len(modelSequence)):
            if modelSequence[i] != '-':
                positions.append(i)	
                column = [modelSequence[i]]
                transpose.append(column)

        for key in self.species:
            if key == self.modelOrg:
                continue 
            sequence = alignment[key]
            for i in range(len(positions)):
                transpose[i].append(sequence[positions[i]])
        
        return transpose

    
    def makeMiMatrix(self):
        """
        Calls the necessary functions to calculate the mutual information between the two proteins. 
        Uses mulit-processing to optimize runtime, but as multi-processing returns results in an 
        order other than the one it was called in, the returned data also needs to be correctly orderd.
        Saves results to a csv.  
        """
        p = Pool()

        basic_mi_start = time.time()
        argList = [i for i in range(len(self.transpose1))]
        self.miMatrix = MiMatrix(p.map(self.makeRow, argList)) 
        p.close()

        self.miMatrix.orderMatrix()
        self.miMatrix.correctScores()
        
        print("Basic MI time:", time.time()-basic_mi_start)

    
    def recordHighestValue(self, fname, p1, p2, taxGroup):
        """
        For mega-analyses involving more than one protein-protein pair, in which the user wishes to only store the 
        highest score from each pairing in a tsv file as a table with the following columns: the name of the first protein,
        the name of the second protein, the score, the position on the first protein corresponding to the max score, the 
        position on the second protein corresponding to the max score, and the likelihood of interaction. Note that the 
        likelihood of interaction is not known if the minPx and percentAboveRandom values are set by the user. Also note
        that this function appends to a existent file rather than overwriting it if that file does exist. 

        Args:
            fname (string): The name of the output file
            p1 (string): The name of the first protein
            p2 (string): The name of the second protein
            taxGroup (string): The taxonomic group that the data is from
        """
        getLikelihood = True
        if taxGroup == 'v' or taxGroup == 'vertebrates':
            thresholds = {480.76335120970106: 0.9967920187830541, 151.51649042416005: 0.8161532125205931, 79.28735697652426: 0.8443779720816076, 45.61458962695341: 0.8161532125205931, 30.18496425069254: 0.6355705515962872, 20.9852341587332: 0.7146303729268318, 15.044388316202244: 0.631461211223173, 10.925487048296537: 0.6233999112853246, 7.851299350820726: 0.5722610958215775, 6.166235010210569: 0.5790468732193523, 4.620538882440486: 0.5466372054916955}
        elif taxGroup == 'b' or taxGroup == 'bacteria':
            thresholds = {44.399839098866174: 0.9580806529465971, 32.01057839583139: 0.7620271888427133, 23.930940257329848: 0.7558249205295701, 20.136639306157356: 0.8155672823218998, 18.264811466880378: 0.885564605247002, 17.03460914998562: 0.845261348542855, 15.001936391295194: 0.8014891094977034, 13.194554167079934: 0.7683320904797415, 11.065147470225472: 0.6989372355202378, 10.43074661430499: 0.885564605247002, 9.460142689310178: 0.8084669132757912, 8.339216723912237: 0.7206328059950042, 7.089151582107779: 0.6639560574444582, 6.40719242029028: 0.7262687969924811, 5.546289684960158: 0.6454952267303103, 4.693976626377164: 0.5532061771323379, 4.1379838378616345: 0.599645262311892, 3.748001781578582: 0.6410393150238497, 3.206328273973718: 0.6280332056194126, 2.316631254998976: 0.3839479007701317}
        else:
            getLikelihood = False

        maxScore = 0
        position1 = 0
        position2 = 0
        for i in range(len(self.miMatrix.matrix)):
            for j in range(len(self.miMatrix.matrix[0])):
                if self.miMatrix.matrix[i][j] > maxScore:
                        maxScore = self.miMatrix.matrix[i][j]
                        position1 = i
                        position2 = j
        if getLikelihood:
            keys = list(thresholds.keys())
            keys.sort()
            likelihood = 0.5
            for key in keys:
                if maxScore >= key:
                    likelihood = thresholds[key]

        else:
            likelihood = np.nan

        writeHeader = False
        if not os.path.exists(fname):
            writeHeader = True
        outfile = open(fname, 'a+')
        if writeHeader:
            headers = ['Protein1', 'Protein2', 'MaxScore', 'Position1', 'Position2', 'LikelihoodOfInteraction']
            outfile.write("\t".join(headers) + "\n")
        output = [p1, p2, str(maxScore), str(position1), str(position2), str(likelihood)]
        outfile.write("\t".join(output) + "\n")
        outfile.close()



    def makeRow(self, index):
        """
        Makes a single row of the mutual information matrix by comparing a single position on protein 1 
        to every position on gene 2. 
        Args:
            index (int): The index of the position on protein 1 that will be compared to every position on gene 2
        Returns:
            tuple (tuple of (int, list of float)): The calculated row of mutual information coupled with an int indicating 
                                                    where this row should fall in the eventual MI matrix, used because of multi-threading
        """
        X = self.transpose1[index] 
        rowNum = index
        size = len(self.transpose2)	
        row = np.empty(size)
        for Y in range(size):	
            row[Y] = self.calcMutualInfo(X, self.transpose2[Y])     
        tuple = (rowNum, row)
        return tuple

        
    def calcMutualInfo(self,X,Y):
        """
        Calculates the mutual information between two protein positions. Mutual information with our method is 
        calculated as the sum of partial mutual information scores, where each type of amino acid pairing between 
        the two positions (i.e. where 'A' occurs at the same indices in X as 'M' occurs in Y) has its own partial MI score. 
        Each partial mutual information score is calculated as p(xy) * log(p(xy)/(px*py)) where px is the percentage 
        of position X made up by the first amino acid of interest, py is the percentage of position Y made up by the
        second amino acid of interest, and pxy is the percentage of instances in which both occur at the same indices.
        Partial MI scores are not considered if their components fall below certain thresholds that were either the 
        thresholds we found to be best for precision, or those set by the user. px and py must meet a certain threshold (minPx), and
        pxy must be better than what may randomly occur given px and py by a certain amount. Further, pairings involving
        ambiguous amino acids or a deletion (-) are not considered. 

        Args:
            X (list of char): A list of char representing the amino acids manifest across different species at a certain position on the first protein
            Y (list of char): Same as X, but for a certain position on the second protein
        
        Returns:
            mutualInfo (float): The mutual information score calculated as per the above method.
        """                               
        ambiguousAminoAcids = "X J Z B -"
        Xresidues = self.makeDictionary(X)
        Yresidues = self.makeDictionary(Y)	
        px = 0.0 
        py = 0.0
        pxy = 0.0
        mutualInfo = 0.0

        for residue1 in Xresidues:
            if residue1 in ambiguousAminoAcids:
                continue

            residue1positions = Xresidues[residue1]
            px = float(len(residue1positions))/float(len(X))
            if px < self.minPx:                                                             
                continue

            for residue2 in Yresidues:
                if residue2 in ambiguousAminoAcids:
                    continue

                residue2positions = Yresidues[residue2]
                py = float(len(residue2positions))/float(len(Y))
                if py < self.minPx:
                    continue                                                       

                pxy = len([value for value in residue1positions if value in residue2positions])/len(Y)
                betterThanRandom = max([self.percentOfRand * px * py, px*py + 0.05])
                if pxy != 0 and pxy >= self.interactionRequirement/len(Y) and pxy >= betterThanRandom: 
                    partialMi = pxy*math.log(pxy/(px*py), 10)
                    mutualInfo += partialMi
                
                
        return mutualInfo
    

    def makeDictionary(self, column):
        """
        Makes a dictionary of the lists of the list of chars that represents a position on a protein as manifest across 
        different species. This dictionary will have the different manifested amino acids as keys and a list of indices
        at which the ocurr as values. This makes it easy to see, between the two different positions, where amino acids
        co-occur. 

        Args:
            column (list of char): The list of char representing a position on a protein as manifest across different species
        
        Returns:
            dictionary (dict): The dictionary described above
        """
        dictionary = {}
        keys = list(set(column))
        for key in keys:
                dictionary[key] = [i for i, x in enumerate(column) if x == key]
        return dictionary


    def removeNeg(self):
        """
        Removes all negative values by changing them to 0 for heatmap creation, allowing users to more easily see 
        areas of interest. Does not affect the final CSV or other calculations
        """
        for i in range(len(self.miMatrix.matrix)):	    
            line = self.miMatrix.matrix[i]
            for j in range(len(line)):
                if j < 0:
                    self.miMatrix.matrix[i][j] = 0
	
    		
    def makeHeatMap(self, removeNeg, color, p1, p2, imagename):
        """
        Makes a heatmap of the mutual information matrix. 

        Args
            removeNeg (bool): whether or not to have negative values changed to 0 to see high scores better
            color (string): the matplotlib color scheme to use for the heatmap
        """
        if removeNeg:
            self.removeNeg()
        
        multiArray = self.miMatrix.matrix
        plt.suptitle("  Mutual Information" + p1 + " vs " + p2) 
        plt.ylabel(p1)
        plt.xlabel(p2)
        plt.ylim(0, len(multiArray))
        plt.imshow(multiArray, cmap=color, interpolation="nearest")
        plt.colorbar() 
        plt.savefig(imagename)
        

def setArgs(parser):
	"""
	Sets up the argparse and checks for errors in user input.

	Args:
		parser (ArgumentParser): the argument parser

	Return:
		args: the arguments
	"""
	parser.add_argument("-f1", "--fasta1", help="The multiple sequence alignment of the protein sequence of the first gene", required=True)
	parser.add_argument("-f2", "--fasta2", help="The multiple sequence alignment of the protein sequence of the second gene", required=True)
	parser.add_argument("-d", "--dataname", help="The filepath to the output csv. Default name made from protein names and written to the current directory", required=False)
	parser.add_argument("-i", "--imagename", help="The filepath to the heat map png. Default name made from the protein names and written to the current directory", required=False)
	parser.add_argument("-c","--color", help="The matplotlib color scheme of the heatmap", nargs="?", type=str, default="viridis")	
	parser.add_argument("-m", "--noheatmap", action="store_true", help="Flag to not create a heatmap", required=False)
	parser.add_argument("-t", "--threads", type=int, help="The number of threads, default all", default=40, required=False)
	parser.add_argument("-n", "--noneg", help="Change all negative values to zero exclusively for the creation of the heatmap so that high values stand out more", action="store_true", required=False)
	parser.add_argument("-p", "--minPx", help="Minimum percentage of the domain that a residue needs to make up to count towards the MI score.", type=float, required=False)
	parser.add_argument("-r", "--percentAboveRandom", help="Minimum percentage above random that a pairing should occur to be factored into MI scores. Default 35", type=float, required=False)
	parser.add_argument("-s", "--highestscore", help="Flag to record only the highest mutual information score and position of said score in a tsv. Also provide the filename (e.g. -s allMaxes.tsv)", required=False)
	parser.add_argument("-p1", "--protein1", help="The name of the first protein being compared", required=True)
	parser.add_argument("-p2", "--protein2", help="The name of the second protein being compared", required=True)
	parser.add_argument("-x", "--taxonomic", help="The taxonomic group of the input files. Flag to use the best values for percentAboveRandom and minPx as determined by our research. Input 'b' or 'bacteria', 'v' or 'vertebrates", required=False)
	args = parser.parse_args()

	colors = "viridis plasma inferno magma Greys Purples Blues Greens Oranges Reds YlOrBr YlOrRd OrRd PuRd RdPu BuPu GnBu PuBu YlGnBu PuBuGn BuGn YlGn binary gist_yarg gist_gray gray bone pink spring summer autumn winter cool Wistia hot afmhot gist_heat copper"
	taxGroups = ["v", 'vertebrates', 'b', 'bacteria']
	if args.color not in colors: 
		print("Please enter a valid matplotlib color scheme (viridis,plasma,inferno,magma,Greys,Purples,Blues,Greens,Oranges,Reds,YlORBr,YlOrRd,OrRd,PuRd,RdPu,BuPu,GNBU,PuBU,YlGnBU,PuBUGn,BuGn,YlGn,binary,gist_yarg,gist_gray,gray,bone,pink,spring,summer,autumn,winter,cool,Wistia,hot,afmhot,gist_heat,copper)")
		exit()
    
	if args.dataname:
		try:
			open(args.dataname, 'w')
		except OSError:
			print("No such file path for --dataname. Please enter a valid file path")
			exit()
		#outfile.close()
	else:
		args.dataname = args.protein1 + "_" + args.protein2 + ".csv"
	
	if args.imagename:
		try:
			open(args.imagename, 'w')
		except OSError:
			print("No such file path for --imagename. Please enter a valid file path")
			exit()
		#outfile.close()
	else:
		args.imagename = args.protein1 + "_" + args.protein2 + ".png"
    
	if args.highestscore and "/" in args.highestscore:
		fname = args.highestscore.split("/")[-1]
		directory = args.highestscore[0:args.highestscore.find(fname)]
		if not os.path.exists(directory):
			print("No such file path for --highestscore. Please enter a valid file path")
			exit()

	if args.percentAboveRandom and not args.minPx:
		print("If using --percentAboveRandom please also include --minPx")
		exit()
    
	if args.minPx and not args.percentAboveRandom:
		print("If using --minPx please also include --percentAboveRandom")
		exit()

	if args.taxonomic and args.percentAboveRandom:
		args.taxonomic = None

	if not args.percentAboveRandom and (args.taxonomic == 'v' or args.taxonomic == 'vertebrates'):
		args.percentAboveRandom = 35.0
		args.minPx = 17.0
    
	if not args.percentAboveRandom and (args.taxonomic == 'b' or args.taxonomic == 'bacteria'):
		args.percentAboveRandom = 22.0
		args.minPx = 29.0

	if args.taxonomic and args.taxonomic not in taxGroups:
		print("Invalid option entered for --taxonomic. Please enter a valid option (v, b, vetebrates, bacteria")
		exit()
    
	if not args.taxonomic and not args.minPx:
		print("Please either choose a taxonomic group or set the minPx and percentAboveRandom values yourself.")
		exit()
	
	return args



if __name__ == '__main__':
    start_time = time.time()
    freeze_support()
    print("freeze support: %s seconds" %(time.time() - start_time))
    parser_start = time.time()
    parser = argparse.ArgumentParser()
    args = setArgs(parser)
    print("argParse: %s seconds" %(time.time() - parser_start))

    calculator = MutualInfoCalculator(args.fasta1, args.fasta2, args.minPx, args.percentAboveRandom)
    prep_start = time.time()
    calculator.prepareAlignments()
    print("alignment prep: %s seconds" %(time.time() - prep_start))
    calculator.makeMiMatrix()
    
    if not (args.highestscore is None):
        calculator.recordHighestValue(args.highestscore, args.protein1, args.protein2, args.taxonomic)
    else:
        calculator.miMatrix.makeCsv(args.dataname)
        if not args.noheatmap:
            calculator.makeHeatMap(args.noneg, args.color, args.protein1, args.protein2, args.imagename) 
        print("total time: %s seconds" %(time.time() - start_time))
