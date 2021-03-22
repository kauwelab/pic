#${1} is the path to the directory with the fastas
python2 findOverlaps.py ${1}
python2 getSameSpecies.py ${1}
