#! /usr/bin/env python
import os
import sys
allFasta = []
path = sys.argv[1]
allFasta = os.listdir(path)
if path[-1] != '/':
    path += '/'
allInputFiles = [path +i for i in allFasta]

notIncluded = dict() #gene->gene->set of species not included
fasta = dict() #gene->species->sequence
allSpecies = dict() #gene->set of species
count = 0
for firstRead in allInputFiles:
    count +=1
    if not ".fasta" in firstRead:
        continue
    gene = firstRead.split("/")[-1][0:-6]
    inputF = open(firstRead)
    allSpecies[gene] = set()
    head = inputF.readline()
    while head != "":
        allSpecies[gene].add(head[1:].strip())
        inputF.readline()
        head = inputF.readline()
    inputF.close()
count=0
usedGenes = set()
output = open(path + "excludedSpecies",'w') #gene1\tgene2\texcludedFrom1\texcludedFrom2\n
for gene1,species1 in allSpecies.items():
    count +=1
    print(count)
    usedGenes.add(gene1)
    for gene2,species2 in allSpecies.items():
        if gene2 in usedGenes:
            continue

        excluded1 = sorted(list(species1-species2))
        #if len(species1)-len(excluded1) <100:
            #continue
        excluded2 = sorted(list(species2-species1))
        #if len(species2)-len(excluded2) <100:
            #continue
        output.write(gene1 + "\t" +gene2 + "\t" + str(excluded1) +"\t" +str(excluded2) +"\n")
output.close()

