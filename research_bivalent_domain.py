# Bivalent domain identification
# H3K4me3/H3K27me3 analysis
# Bongsoo Park, Ph.D
# National Institute on Aging
# ESCA/TGB/NIA

#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import operator
import math

f = open("./summary_reports/summary_h3k4me3_annotation.All.final.csv")
cnt = 0
# H3K4me3 genes associated with TSS-proximal sites
# The key: NCBI RefSeq ID
# The value: GeneSymbol
gene1 = {}
gene1_young = {}
gene1_old = {}
gene1_log2FC = {}
gene1_peak = {}
gene1_peak2 = {}
x3 = {} #PeakOrigin H3K4me3
y3 = {} #PeakOrigin H3K27me3
z3 = {} #PeakOrigin for union set
for line in f:
	line = line.strip()
	data = line.split(",")
	# Excluding header
	# Line number
	if cnt > 0:
		# TSS-proximal sites were defined 1Kbp from TSS
		if data[3] == "TSS-proximal":
			new_line = data[3]
			the_key = data[4]
			gene1.update({the_key:data[5]})
			gene1_young.update({the_key:str((float(data[8])+float(data[9]))/2)})
			gene1_old.update({the_key:str((float(data[10])+float(data[11]))/2)})
			gene1_log2FC.update({the_key:data[7]})
			gene1_peak.update({the_key:data[0]})
		x3.update({the_key:data[12]})
		z3.update({the_key:data[12]})
	cnt += 1
f.close()
f = open("./summary_reports/summary_h3k27me3_annotation.All.final.csv")
cnt = 0
# H3K27me3 genes associated with TSS-proximal sites
# The key: NCBI RefSeq ID
# The value: GeneSymbol
gene2 = {}
gene2_young = {}
gene2_old = {}
gene2_log2FC = {}
gene2_peak = {}
gene2_peak2 = {}
for line in f:
	line = line.strip()
	data = line.split(",")
	if cnt > 0:
		if data[3] == "TSS-proximal":
			new_line = data[3]
			the_key = data[4]
			gene2.update({the_key:data[5]})
			gene2_young.update({the_key:str((float(data[8])+float(data[9]))/2)})
			gene2_old.update({the_key:str((float(data[10])+float(data[11]))/2)})
			gene2_log2FC.update({the_key:data[7]})
			gene2_peak.update({the_key:data[0]})
		y3.update({the_key:data[12]})
		z3.update({the_key:data[12]})
	cnt += 1
f.close()

#Finding union geneset by merging H3K4me3/H3K27me3 datasets
gene = {}
for the_key in gene1:
	gene.update({the_key:gene1[the_key]})
for the_key in gene2:
	gene.update({the_key:gene2[the_key]})

x = {} #Annotation H3K4me3
z = {} #Annotation for union set
x2 = {} #Fold change H3K4me3
z2 = {} #Fold change for union set
f = open("./summary_reports/summary_h3k4me3_annotationFC1.5All.final.csv")
fw = open("./summary_reports/summary_h3k4me3_annotationFC1.5_TSS.final.csv","w")
cnt = 0
for line in f:
	line = line.strip()
	data = line.split(",")
	if cnt > 0:
		if data[3] == "TSS-proximal":
			fw.write(line + "\n")
			x.update({data[4]:data[5]})
			z.update({data[4]:data[5]})
			x2.update({data[4]:data[6]})
			z2.update({data[4]:data[6]})
			x3.update({data[4]:data[12]})
			z3.update({data[4]:data[12]})
	else:
		fw.write(line + "\n")
	cnt += 1
f.close()

y = {} #Annotation H3K27me3
y2 = {} #Fold change H3K27me3
f = open("./summary_reports/summary_h3k27me3_annotationFC1.5All.final.csv")
fw = open("./summary_reports/summary_h3k27me3_annotationFC1.5_TSS.final.csv","w")
cnt = 0
for line in f:
	line = line.strip()
	data = line.split(",")
	if cnt > 0:
		if data[3] == "TSS-proximal":
			fw.write(line + "\n")
			y.update({data[4]:data[5]})
			z.update({data[4]:data[5]})
			y2.update({data[4]:data[6]})
			z2.update({data[4]:data[6]})
			y3.update({data[4]:data[12]})
			z3.update({data[4]:data[12]})
	else:
		fw.write(line + "\n")
	cnt += 1
f.close()

f = open("DEG2_none_ercc_YHSC-OHSC_all.txt")
deg_up = {}
deg_dn = {}
deg_no = {}
cnt = 0
for line in f:
	line = line.strip()
	data = line.split("\t")
	if cnt > 0:
		try:
			if float(data[6]) < 0.05 and float(data[2]) > 0.263:
				deg_up.update({data[0]:data[2]+"\t"+data[8]})	
			elif float(data[6]) < 0.05 and float(data[2]) < -0.263:
				deg_dn.update({data[0]:data[2]+"\t"+data[8]})	
			else:
				deg_no.update({data[0]:data[2]+"\t"+data[8]})	
		except:
			pass
	cnt += 1
f.close()

print "Annotation TSS-proximal H3K4me3, H3K27me3, Union"
print len(x), len(y), len(z)
print "FoldChange TSS-proximal H3K4me3, H3K27me3, Union"
print len(x2), len(y2), len(z2)
print "Annotation All H3K4me3, H3K27me3, Union"
print len(x3), len(y3), len(z3)

