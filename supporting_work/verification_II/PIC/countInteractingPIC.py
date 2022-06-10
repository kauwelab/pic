import pandas as pd
import sys

threshold = float(sys.argv[1])

df = pd.read_csv("output/PIC_output_unknown.tsv", sep="\t")
aboveT = df[df["LikelihoodOfInteraction"] > threshold] 
print("Interactions found in unknown-to-interact pairs:",len(aboveT.index))

df = pd.read_csv("output/PIC_output_known.tsv", sep="\t")
aboveT = df[df["LikelihoodOfInteraction"] > threshold]
print("Interactions found in known-to-interact pairs:", len(aboveT.index))
