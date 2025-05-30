# Annotation version 2
# Report all sites including TSS/TES/Genebody/Intergenic
import operator
import math
import re

# If one peak shared the same gene within close distance, we will use only bigger size transcript.
# Okay.
# Rgmb gene
# Tas2r115 (too small), check genebody first then 2K

narrowPeak = ["h3k4me3","h3k27ac","h3k27me3"]
for file_name in narrowPeak:
	f1 = open("./summary_reports/summary_"+file_name+"_annotation.All.csv")
	f2 = open("./summary_reports/summary_"+file_name+"_annotation.TSS.csv")
	f3 = open("./summary_reports/summary_"+file_name+"_annotationFC1.5All.csv")
	fw = open("./summary_reports/summary_"+file_name+"_annotation.All.final.csv","w")
	fw2 = open("./summary_reports/summary_"+file_name+"_annotation.TSS.final.csv","w")
	fw3 = open("./summary_reports/summary_"+file_name+"_annotationFC1.5All.final.csv","w")
	fw4 = open("./summary_reports/summary_"+file_name+"_annotationFC1.5All.final.removeDUP.csv","w")
	fw.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
	fw2.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
	fw3.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
	fw4.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
	
	cnt = 0
	cnt2 = 0
	duplicate_peak = {}
	duplicate_gene_all = {}
	for line in f1:
		line = line.strip()
		data = line.split(",")
		if cnt > 0:
			the_key = data[0] + "\t" + data[4] + "\t" + data[5]
			if duplicate_gene_all.has_key(the_key):
				cnt2 += 1
				#print cnt2, file_name, the_key
				pass
			else:
				duplicate_gene_all.update({the_key:1})
				fw.write(line + "\n")	
		cnt += 1
	
	cnt = 0
	cnt2 = 0
	duplicate_gene_tss = {}
	for line in f2:
		line = line.strip()
		data = line.split(",")
		if cnt > 0:
			the_key = data[0] + "\t" + data[4] + "\t" + data[5]
			if duplicate_gene_tss.has_key(the_key):
				cnt2 += 1
				#print cnt2, file_name, the_key
				pass
			else:
				duplicate_gene_tss.update({the_key:1})
				fw2.write(line + "\n")	
		cnt += 1
	cnt = 0
	cnt2 = 0
	duplicate_gene_all = {}
	for line in f3:
		line = line.strip()
		data = line.split(",")
		if cnt > 0:
			the_key = data[0] + "\t" + data[4] + "\t" + data[5]
			if duplicate_gene_all.has_key(the_key):
				cnt2 += 1
				#print cnt2, file_name, the_key
				pass
			else:
				duplicate_gene_all.update({the_key:1})
				fw3.write(line + "\n")	
				if duplicate_peak.has_key(the_key):
					pass
				else:
					duplicate_peak.update({the_key:1})
					fw4.write(line + "\n")	
				if duplicate_peak.has_key(the_key):
		cnt += 1
	
	cnt = 0
		
	f1.close()
	f2.close()
	f3.close()
	fw4.close()
	fw3.close()
	fw2.close()
	fw.close()

