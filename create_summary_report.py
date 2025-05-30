# Import libraries
import operator
import math
import re

# 1/24/2022 TSS/Genebody/TES/Intergenic matching update
# Bongsoo Park, Ph.D, National Institute on Aging
# H3K4me3/ATAC/H3K27ac - narrowPeak
# H3K27me3 - IDEAS-repressive mark + extension with low/broad peaks, comparison with SICER dataset later
# H3K27me3 - merged overlapped peaks (4k) and exclude small singleton peaks (<2k size)

# Fix error belows
# Gm13483
# Cacnb4
# Hoxd1
# Wipf1

# H3K27me3 bivalent
# Rgmb 
# Lsr

# 2/2/2022 Update: excluding chrUn chromosomes

# NCBI RefSeq Jan 2020 version
# We use the unique TSS site as a reference due to the duplicated outcomes with the original dataset
# When two transcripts shared the same TSS, we will keep the transcript that has a larger size (TES-TSS).
# 2/7/2022 Update: excluding duplicated TSS sites

f = open("../00_reference_files/UCSC-RefSeq.unique.TSS.bed")
fw = open("../00_reference_files/UCSC-RefSeq.TSS2kb.upstream.bed","w")
tss = {}
tss_size = {}
tss_strand = {}
tes = {}
genebody = {}
for line in f:
	line = line.strip()
	data = line.split("\t")
	chrom = data[0].split("_")[0]
	if chrom != "chrUn":
		the_size = int(data[2]) - int(data[1])
		data[3] = data[3] + "|" + str(the_size)
		if data[5] == "+":
			# Searching TSS/TES/Genebody
			for i in range(int(data[1])-500, int(data[1])+500, 20):
				the_key = data[0] + ":" + str((int(data[1])/20)*20)
				tss.update({the_key:data[3]})
				tss_size.update({the_key:the_size})
			the_key = data[0] + ":" + str((int(data[2])/20)*20)
			tes.update({the_key:data[3]})
			tss_strand.update({the_key:"+"})
			for i in range(int(data[1])+500, int(data[2]), 50):
				the_key3 = data[0] + ":" + str((i/50)*50)
				genebody.update({the_key3:data[3]})
			fw.write(data[0] + "\t" + str(int(data[1])-2000) + "\t" + str(int(data[1])) + "\t" + data[3] + "\t" + data[4] + "\t" + data[5] + "\n")
		else:
			# Searching TSS/TES/Genebody
			for i in range(int(data[2])-500, int(data[2])+500, 20):
				the_key = data[0] + ":" + str((int(data[2])/20)*20)
				tss.update({the_key:data[3]})
			tss_size.update({the_key:the_size})
			the_key = data[0] + ":" + str((int(data[1])/20)*20)
			tes.update({the_key:data[3]})
			tss_strand.update({the_key:"-"})
			for i in range(int(data[1]), int(data[2])-500, 50):
				the_key3 = data[0] + ":" + str((i/50)*50)
				genebody.update({the_key3:data[3]})
			fw.write(data[0] + "\t" + str(int(data[2])) + "\t" + str(int(data[2])+2000) + "\t" + data[3] + "\t" + data[4] + "\t" + data[5] + "\n")
fw.close()

f = open("../00_reference_files/UCSC-RefSeq.unique.TSS.bed")
tss_2kb = {}
tss_2kb_strand = {}
for line in f:
	line = line.strip()
	data = line.split("\t")
	chrom = data[0].split("_")[0]
	if chrom != "chrUn":
		the_size = int(data[2]) - int(data[1])
		data[3] = data[3] + "|" + str(the_size)
		if data[5] == "+":
		# Tas2r115 (too small), check genebody first then 2K
			the_key = data[0] + ":" + str((int(data[1])/50)*50)
			for i in range(int(data[1]), int(data[1])+2000, 50):
				the_key2 = data[0] + ":" + str((i/50)*50)
				if genebody.has_key(the_key2):
					if genebody[the_key2] == data[3]:
						tss_2kb.update({the_key2:data[3]})
						tss_2kb_strand.update({the_key2:"+"})
		else:
			the_key = data[0] + ":" + str((int(data[2])/50)*50)
			for i in range(int(data[2])-2000, int(data[2]), 50):
				the_key2 = data[0] + ":" + str((i/50)*50)
				if genebody.has_key(the_key2):
					if genebody[the_key2] == data[3]:
						tss_2kb.update({the_key2:data[3]})
						tss_2kb_strand.update({the_key2:"-"})
f.close()

print "TSS number:", len(tss)
print "TES number:", len(tes)
print "Genebody:", len(genebody)
print "TSS 2kb:", len(tss_2kb)

array_name = ["h3k4me3", "h3k27ac", "h3k27me3", "atac"]
cnt = 0
for ele in array_name:
	f = open("preliminary_reports/summary_"+array_name[cnt]+".bed.annPeaks.txt")
	x = {}
	y = {}
	cnt2 = 0
	for line in f:
		line = line.strip()
		data = line.split("\t")
		if cnt2 > 0:
			the_key = data[0]
			the_pos = data[7].split(" ")[0]
			the_dis = data[9]
			try:
				the_gene = data[15]
				x.update({the_key:the_pos+","+the_gene+","+the_dis})
			except:
				x.update({the_key:the_pos+",.,."})
			flag = 0
			if ele == "h3k4me3" or ele == "h3k27ac" or ele == "atac":
				coor = data[0].split(":")
				chrom = coor[0]
				# narrowPeak, search +/- 1kb from peaks
				start = int(coor[1].split("-")[0]) - 1000
				end = int(coor[1].split("-")[1]) + 1000
				for i in range(start, end, 20):
					the_key2 = chrom+":"+str((i/20)*20)
					if tss.has_key(the_key2):
						if y.has_key(the_key):
							y[the_key] += "," + tss[the_key2]
						else:
							y.update({the_key:"TSS-proximal,"+tss[the_key2]})				
						flag = 1
			else:
				#H3K27me3: checking TSS-proximal site
				coor = data[0].split(":")
				chrom = coor[0]
				start = int(coor[1].split("-")[0])  
				end = int(coor[1].split("-")[1]) 
				for i in range(start, end, 20):
					the_key2 = chrom+":"+str((i/20)*20)
					if tss.has_key(the_key2):
						if y.has_key(the_key):
							y[the_key] += "," + tss[the_key2]
						else:
							y.update({the_key:"TSS-proximal,"+tss[the_key2]})				
						flag = 1
				if flag == 0:
					# If H3K27me3 peaks were not overlaped with TSS/ searching downstream 2Kb
					# But we should assign only one nearest TSS
					start = int(coor[1].split("-")[0])  
					end = int(coor[1].split("-")[1]) 
					for i in range(start, end, 50):
						the_key2 = chrom+":"+str((i/50)*50)
						if tss_2kb.has_key(the_key2):
							if tss_2kb_strand[the_key2] == "+":
								y.update({the_key:"TSS-proximal,"+tss_2kb[the_key2]})				
								flag = 1
								break
					if flag == 0:	
						for i in range(end, start, -50):
							the_key2 = chrom+":"+str((i/50)*50)
							if tss_2kb.has_key(the_key2):
								if tss_2kb_strand[the_key2] == "-":
									y.update({the_key:"TSS-proximal,"+tss_2kb[the_key2]})				
									flag = 1
									break
			# Checking TES/Genebody/Intergenic
			if flag == 0:
				# Checking TES
				coor = data[0].split(":")
				chrom = coor[0]
				start = int(coor[1].split("-")[0]) 
				end = int(coor[1].split("-")[1]) 
				for i in range(start, end, 20):
					the_key2 = chrom+":"+str((i/20)*20)
					if tes.has_key(the_key2):
						y.update({the_key:"TES,"+tes[the_key2]})				
						flag = 1
				if flag == 0:
					# Checking Genebody
					for i in range(start, end, 50):
						the_key2 = chrom+":"+str((i/50)*50)
						if genebody.has_key(the_key2):
							y.update({the_key:"Genebody,"+genebody[the_key2]})				
							flag = 1
					if flag == 0:
						y.update({the_key:"Intergenic,NA"})
		cnt2 += 1	
	f.close()

	fo = open("./preliminary_reports/summary_"+array_name[cnt]+".csv")
	fw = open("./summary_reports/summary_"+array_name[cnt]+"_annotation.csv","w")
	if ele == "atac":
		# ATAC 4x8 comparison
		fw.write("Coordination,Young1,Young2,Young3,Young4,Old1,Old2,Old3,Old4,Old5,Old6,Old7,Old8,FoldChange(Old/Young),log2FC,Homer-annotation,GeneName,Distance2TSS,TSS-proximal,RefSeq-annotation,MultipleTSS\n")
		for line in fo:
			line = line.strip()
			data = line.split(",")
			coor = data[0].split(":")
			chrom = coor[0].split("_")[0]
			if chrom != "chrUn":
				auc_value1 = float(data[1]) + float(data[2]) + float(data[3]) + float(data[4])
				auc_value2 = float(data[5]) + float(data[6]) + float(data[7]) + float(data[8])
				auc_value2 += float(data[9]) + float(data[10]) + float(data[11]) + float(data[12])
				auc_value = auc_value1 + auc_value2
				# Extremtly high or low signal filtering
				if auc_value > 0.2 and auc_value < 2000:
					fold_change = auc_value2 / auc_value1
					log2FC = math.log(fold_change,2)
					fw.write(line + "," + str(fold_change)+","+str(log2FC) + "," + x[data[0]]+"," + y[data[0]] + "\n")
		fw.close()
	else:
		# ChIP-seq 2x2 comparison
		fw.write("Coordination,Young1,Young2,Old1,Old2,FoldChange(Old/Young),log2FC,Homer-annotation,GeneName,Distance2TSS,TSS-proximal,RefSeq-annotation,MultipleTSS\n")
		for line in fo:
			line = line.strip()
			data = line.split(",")
			coor = data[0].split(":")
			chrom = coor[0].split("_")[0]
			if chrom != "chrUn":
				auc_value = float(data[1]) + float(data[2]) + float(data[3]) + float(data[4])
				# Extremtly high or low signal filtering
				#if auc_value > 0.3 and auc_value < 800:
				try:
					fold_change = (float(data[3]) + float(data[4])) / (float(data[1]) + float(data[2]))
				except:
					fold_change = (float(data[3]) + float(data[4])) / 0.01
				log2FC = math.log(fold_change,2)
				fw.write(line + "," + str(fold_change)+","+str(log2FC) + "," + x[data[0]]+"," + y[data[0]] + "\n")
		fw.close()
	fo.close()
	cnt += 1

# Annotation version 2
# Summary Report Step #2

f = open("../00_reference_files/mm10-blacklist.v2.bed")
blacklist = {}
for line in f:
	line = line.strip()
	data = line.split("\t")
	for i in range(int(data[1]), int(data[2]), 50):
		the_key = data[0] + ":" + str(i)
		data[3] = re.sub(" ","-",data[3])
		blacklist.update({the_key:data[3]})
f.close()

# Check the origin of peaks (YHSC or OHSC or Common)
f = open("peak_calling_MACS/ChIP_YOHSC_H3K4me3_Created_R1.dedup.uniq.bam_peaks.narrowPeak.bed")
h3k4me3_peaks = {}
for line in f:
	line = line.strip()
	data = line.split("\t")
	the_key = data[3]
	the_value = data[4]
	h3k4me3_peaks.update({the_key:the_value})
f.close()
f = open("peak_calling_MACS/ChIP_YOHSC_H3K27ac_Created_R1.dedup.uniq.bam_peaks.narrowPeak.bed")
h3k27ac_peaks = {}
for line in f:
	line = line.strip()
	data = line.split("\t")
	the_key = data[3]
	the_value = data[4]
	h3k27ac_peaks.update({the_key:the_value})
f.close()
f = open("peak_calling_IDEAS/ChIP_YOHSC_H3K27me3_Created_R1.dedup.uniq.bam_peaks.IDEAS.bed")
h3k27me3_peaks = {}
for line in f:
	line = line.strip()
	data = line.split("\t")
	the_key = data[3]
	the_value = data[4]
	h3k27me3_peaks.update({the_key:the_value})
f.close()

narrowPeak = ["h3k4me3","h3k27ac"]
narrow_peak_array = [h3k4me3_peaks,h3k27ac_peaks]
narrow_peak_cnt = 0
for file_name in narrowPeak:
	f1 = open("./summary_reports/summary_"+file_name+"_annotation.All.csv","w")
	f2 = open("./summary_reports/summary_"+file_name+"_annotation.TSS.csv","w")
	fw = open("./summary_reports/summary_"+file_name+"_annotationFC1.5All.csv","w")
	f1.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
	f2.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
	fw.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
	
	test = ["up-regulation","down-regulation"]
	for test_name in test:
		f = open("./summary_reports/summary_"+file_name+"_annotation.csv")
		cnt = 0
		duplicate_gene = {}
		duplicate_gene_all = {}
		narrow_peaks = narrow_peak_array[narrow_peak_cnt]
		for line in f:
			line = line.strip()
			data = line.split(",")
			if cnt > 0:
				new_line = data[10]
				gene = {}
				gene2 = {}
				gene2_size = {}
				for i in range(11, len(data)):
					try:
						# Associated genes
						# RefSeqID|GeneSymbol|GeneSize
						# If there is overlapped gene, then we will assign the longer gene here
						ncbi_ref = data[i].split("|")[0]
						gene_symbol = data[i].split("|")[1]
						gene_size = int(data[i].split("|")[2])
						if gene2.has_key(gene_symbol):
							if gene2_size[gene_symbol] < gene_size:
								saved_ncbi_ref = gene2[gene_symbol]	
								gene.pop(saved_ncbi_ref)
								# if the size of gene archived in gene dictionary is smaller than current one
								# we will replace it with new one
								gene.update({ncbi_ref:gene_symbol})
								gene2.update({ncbi_ref:gene_size})
							else:
								# otherwise, skip
								pass
						else:
							# The-key: RefSeq, The-value: Symbol
							gene.update({ncbi_ref:gene_symbol})
							gene2.update({gene_symbol:ncbi_ref})
							gene2_size.update({gene_symbol:gene_size})
					except:
						gene.update({data[i].split("|")[0]:"NA"})
				for the_key in gene:
					flag = "None"
					coor = data[0].split(":")	
					start = int(coor[1].split("-")[0])
					end = int(coor[1].split("-")[1])
					peak_size = end - start
					the_value = data[1] + "," + data[2] + "," + data[3] + "," + data[4]
					for i in range(start, end, 50):
						the_key_black = coor[0] + ":" + str((i/50)*50)
						if blacklist.has_key(the_key_black):
							flag = blacklist[the_key_black]
					if flag != "High-Signal-Region":
						# Weak peak but high AUC rescue
						young_value = float(data[1]) + float(data[2])
						old_value = float(data[3]) + float(data[4])
						if young_value > 0.5 and old_value > 0.5:
							pass
							#Need to test more
							#narrow_peaks[data[0]] = "YOHSC" 
						if duplicate_gene_all.has_key(data[0]+"|"+gene[the_key]):
							pass
						else:
							duplicate_gene_all.update({data[0]+"|"+gene[the_key]:1})
							f1.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + narrow_peaks[data[0]] + "\n")
							# check duplicate coordination+gene 
							if data[10] == "TSS-proximal":
								f2.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + narrow_peaks[data[0]] + "\n")
						if test_name == "up-regulation":
							if float(data[6]) > 0.58496:
								if float(data[3]) > float(data[1])*1.2 and float(data[3]) > float(data[2])*1.2 and float(data[4]) > float(data[1])*1.2 and float(data[4]) > float(data[2])*1.2:
									if duplicate_gene.has_key(data[0]+"|"+gene[the_key]):
										pass
									else:
										fw.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + narrow_peaks[data[0]] + "\n")
										duplicate_gene.update({data[0]+"|"+gene[the_key]:1})
						else:
							if float(data[6]) < -0.58496:
								if float(data[3])*1.2 < float(data[1]) and float(data[3])*1.2 < float(data[2]) and float(data[4])*1.2 < float(data[1]) and float(data[4])*1.2 < float(data[2]):
									if duplicate_gene.has_key(data[0]+"|"+gene[the_key]):
										pass
									else:
										fw.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + narrow_peaks[data[0]] + "\n")
										duplicate_gene.update({data[0]+"|"+gene[the_key]:1})
			cnt += 1
		f.close()
	narrow_peak_cnt += 1
	f2.close()
	fw.close()

f1 = open("./summary_reports/summary_h3k27me3_annotation.All.csv","w")
f2 = open("./summary_reports/summary_h3k27me3_annotation.TSS.csv","w")
fw = open("./summary_reports/summary_h3k27me3_annotationFC1.5All.csv","w")
f1.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
f2.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
fw.write("Coordination,PeakSize,Blacklist,Annotation,RefSeq,Genename,FoldChange,log2FC,Young1,Young2,Old1,Old2,PeakOrigin\n")
test = ["up-regulation","down-regulation"]
for test_name in test:
	f = open("./summary_reports/summary_h3k27me3_annotation.csv")
	cnt = 0
	duplicate_gene = {}
	duplicate_gene_all = {}
	for line in f:
		line = line.strip()
		data = line.split(",")
		if cnt > 0:
			new_line = data[10]
			gene = {}
			gene2 = {}
			for i in range(11, len(data)):
				try:
					# Associated genes
					# RefSeqID|GeneSymbol|GeneSize
					# If there is overlapped gene, then we will assign the longer gene here
					ncbi_ref = data[i].split("|")[0]
					gene_symbol = data[i].split("|")[1]
					gene_size = int(data[i].split("|")[2])
					if gene2.has_key(gene_symbol):
						if gene2_size[gene_symbol] < gene_size:
							saved_ncbi_ref = gene2[gene_symbol]	
							gene.pop(saved_ncbi_ref)
							# if the size of gene archived in gene dictionary is smaller than current one
							# we will replace it with new one
							gene.update({ncbi_ref:gene_symbol})
							gene2.update({ncbi_ref:gene_size})
						else:
							# otherwise, skip
							pass
					else:
						# The-key: RefSeq, The-value: Symbol
						gene.update({ncbi_ref:gene_symbol})
						gene2.update({gene_symbol:ncbi_ref})
						gene2_size.update({gene_symbol:gene_size})
				except:
					gene.update({data[i].split("|")[0]:"NA"})
			for the_key in gene:
				flag = "None"
				coor = data[0].split(":")	
				start = int(coor[1].split("-")[0])
				end = int(coor[1].split("-")[1])
				peak_size = end - start
				the_value = data[1] + "," + data[2] + "," + data[3] + "," + data[4]
				for i in range(start, end, 50):
					the_key_black = coor[0] + ":" + str((i/50)*50)
					if blacklist.has_key(the_key_black):
						flag = blacklist[the_key_black]
				if flag != "High-Signal-Region":
					# Weak peak but high AUC rescue
					young_value = float(data[1]) + float(data[2])
					old_value = float(data[3]) + float(data[4])
					if young_value > 0.5 and old_value > 0.5:
						pass
						#Need to test more
						#h3k27me3_peaks[data[0]] = "YOHSC" 
					if duplicate_gene_all.has_key(data[0]+"|"+gene[the_key]):
						pass
					else:
						duplicate_gene_all.update({data[0]+"|"+gene[the_key]:1})
						# check duplicate coordination+gene 
						f1.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + h3k27me3_peaks[data[0]] + "\n")
						if data[10] == "TSS-proximal":
							f2.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + h3k27me3_peaks[data[0]] + "\n")
					if test_name == "up-regulation":
						if float(data[6]) > 0.58496:
							if float(data[3]) > float(data[1])*1.2 and float(data[3]) > float(data[2])*1.2 and float(data[4]) > float(data[1])*1.2 and float(data[4]) > float(data[2])*1.2:
								if duplicate_gene.has_key(data[0]+"|"+gene[the_key]):
									pass
								else:
									duplicate_gene.update({data[0]+"|"+gene[the_key]:1})
									# check duplicate coordination+gene 
									if (data[10] == "Intergenic" or data[10] == "Genebody") and peak_size <= 2000:
										pass
									elif flag != "High-Signal-Region":
										fw.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + h3k27me3_peaks[data[0]] + "\n")
					else:
						if float(data[6]) < -0.58496:
							if float(data[3])*1.2 < float(data[1]) and float(data[3])*1.2 < float(data[2]) and float(data[4])*1.2 < float(data[1]) and float(data[4])*1.2 < float(data[2]):
								if duplicate_gene.has_key(data[0]+"|"+gene[the_key]):
									pass
								else:
									duplicate_gene.update({data[0]+"|"+gene[the_key]:1})
									# check duplicate coordination+gene 
									if (data[10] == "Intergenic" or data[10] == "Genebody") and peak_size <= 2000:
										pass
									elif flag != "High-Signal-Region":
										fw.write(data[0]+","+str(peak_size) + "," + flag +","+data[10]+","+the_key+","+gene[the_key]+","+data[5]+","+data[6]+"," + the_value + "," + h3k27me3_peaks[data[0]] + "\n")
		cnt += 1
	f.close()
f2.close()
fw.close()
