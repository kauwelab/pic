import os 

infiles_unknown = os.listdir("inputData/unknown")
outfile_unknown = "output/PIC_output_unknown.tsv"
infiles_known = os.listdir("inputData/known")
outfile_known = "output/PIC_output_known.tsv"
for i in infiles_unknown:
	p1 = i.split("_")[0]
	p2 = i.split("_")[1]
	f1 = "inputData/unknown/" + i + "/" + i
	f2 = "inputData/unknown/" + i + "/" + p2 + "_" + p1
	os.system("python ../../../mi_smart_filters.py -x v -f1 " + f1 + " -f2 " + f2 + " -p1 " + p1 + " -p2 " + p2 + " -s " + outfile_unknown)
for i in infiles_known:
	p1 = i.split("_")[0]
	p2 = i.split("_")[1]
	f1 = "inputData/known/" + i + "/" + i
	f2 = "inputData/known/" + i + "/" + p2 + "_" + p1
	os.system("python ../../../mi_smart_filters.py -g v -xd -xi -f1 " + f1 + " -f2 " + f2 + " -p1 " + p1 + " -p2 " + p2 + " -s " + outfile_known)



	
