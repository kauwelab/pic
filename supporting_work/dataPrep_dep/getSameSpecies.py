#! /usr/bin/env python
import sys
import ast
import subprocess
from multiprocessing import Pool
path = sys.argv[1]
if path[-1] != '/':
    path += '/'
excludedSpecies = open(path +"excludedSpecies") #allExcludedStecies/exxx
print(sys.argv[1])
count = 0
threads =16
def readFiles(line):
    gene1,gene2,excluded1,excluded2=line.strip().split("\t")
    excluded1=set(list(ast.literal_eval(excluded1)))
    excluded2=set(list(ast.literal_eval(excluded2)))
    ###fasta1 = open("fastasFromTable/"+gene1 + ".fasta")
    fasta1 = open(path+gene1 + ".fasta")
    head = fasta1.readline()
    newFasta = ""
    while head !="":
        seq = fasta1.readline()
        if head[1:].strip() not in excluded1:
            newFasta += head
            newFasta += seq
        head =fasta1.readline()
    fasta1.close()
    newFasta= unicode(newFasta,"utf-8")
    #p=subprocess.Popen(["clustalo","--force","--infile=-","--threads=1","-o","sameSpeciesFastas/"+gene1+"/"+gene1+"_"+gene2],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    p=subprocess.Popen(["clustalo","--force","--infile=-","--threads=1","-o",path +gene1+"_"+gene2],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    p.stdin.write(newFasta)
    p.communicate()[0]
    p.stdin.close()
    p.wait()
    ##fasta2 = open("fastasFromTable/"+gene2 + ".fasta")
    fasta2 = open(path+gene2 + ".fasta")
    head = fasta2.readline()
    newFasta = ""
    while head !="":
        seq = fasta2.readline()
        if head[1:].strip() not in excluded2:
            newFasta += head
            newFasta += seq
        head =fasta2.readline()
    newFasta= unicode(newFasta,"utf-8")
    #p=subprocess.Popen(["clustalo","--force","--infile=-","--threads=1","-o","sameSpeciesFastas/"+gene2+"/"+gene2+"_"+gene1],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    p=subprocess.Popen(["clustalo","--force","--infile=-","--threads=1","-o",path+gene2+"_"+gene1],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    p.stdin.write(newFasta)
    p.communicate()[0]
    p.stdin.close()
    p.wait()
    fasta2.close()

allLines = []
for line in excludedSpecies:
    allLines.append(line)
    count +=1
    if count %1000 == 0:
        print(count)
        p=Pool(threads)
        p.map(readFiles,allLines)
        allLines = []
p=Pool(threads)
p.map(readFiles,allLines)
excludedSpecies.close()
