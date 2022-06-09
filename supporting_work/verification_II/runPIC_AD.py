import os 

infiles_unknown = os.listdir("unknown")
outfile_unknown = "PIC_output_unknown.tsv"
infiles_known = os.listdir("known")
outfile_known = "PIC_output_known.tsv"
for i in infiles_unknown:
	p1 = i.split("_")[0]
	p2 = i.split("_")[1]
	f1 = "comparing_scores_AD/unknown/" + i + "/" + i
	f2 = "comparing_scores_AD/unknown/" + i + "/" + p2 + "_" + p1
	os.system("python mi_smart_filters.py -x v -f1 " + f1 + " -f2 " + f2 + " -p1 " + p1 + " -p2 " + p2 + " -s " + outfile_unknown)
for i in infiles_known:
	p1 = i.split("_")[0]
	p2 = i.split("_")[1]
	f1 = "comparing_scores_AD/known/" + i + "/" + i
	f2 = "comparing_scores_AD/known/" + i + "/" + p2 + "_" + p1
	os.system("python mi_smart_filters.py -x v -f1 " + f1 + " -f2 " + f2 + " -p1 " + p1 + " -p2 " + p2 + " -s " + outfile_known)



	
