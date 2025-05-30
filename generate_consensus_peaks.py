import os
import sys

# MACS peak call for narrowPeaks

h3k4me3_peaks = [
"peak_calling_MACS/ChIP_YHSC_H3K4me3_Merged_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_OHSC_H3K4me3_Merged_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_YOHSC_H3K4me3_Merged_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_YOHSC_H3K4me3_Intersect_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_YOHSC_H3K4me3_Created_R1.dedup.uniq.bam_peaks.narrowPeak.bed"
]

h3k27ac_peaks = [
"peak_calling_MACS/ChIP_YHSC_H3K27ac_Merged_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_OHSC_H3K27ac_Merged_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_YOHSC_H3K27ac_Merged_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_YOHSC_H3K27ac_Intersect_R1.dedup.uniq.bam_peaks.narrowPeak",
"peak_calling_MACS/ChIP_YOHSC_H3K27ac_Created_R1.dedup.uniq.bam_peaks.narrowPeak.bed"
]

# Hamp
# SICER FDR<0.05
h3k27me3_peaks_stringent = [
"peak_calling_SICER/ChIP_YHSC_H3K27me3_Merged_R1.dedup.uniq.bam-W200-G1000.bed",
"peak_calling_SICER/ChIP_OHSC_H3K27me3_Merged_R1.dedup.uniq.bam-W200-G1000.bed",
"peak_calling_SICER/ChIP_YOHSC_H3K27me3_Merged_R1.dedup.uniq.bam-W200-G1000.bed",
"peak_calling_SICER/ChIP_YOHSC_H3K27me3_Intersect_R1.dedup.uniq.bam-W200-G1000.bed",
"peak_calling_SICER/ChIP_YOHSC_H3K27me3_Created_R1.dedup.uniq.bam-W200-G1000.bed"
]

# IDEAS + minimum peak hight >0.4
h3k27me3_peaks_relaxed = [
"peak_calling_IDEAS/ChIP_YHSC_H3K27me3_Merged.normalized10M.bigWig_auc_peak.bed",
"peak_calling_IDEAS/ChIP_OHSC_H3K27me3_Merged.normalized10M.bigWig_auc_peak.bed",
"peak_calling_IDEAS/ChIP_YOHSC_H3K27me3_Merged.normalized10M.bigWig_auc_peak.bed",
"peak_calling_IDEAS/ChIP_YOHSC_H3K27me3_Intersect.normalized10M.bigWig_auc_peak.bed",
"peak_calling_IDEAS/ChIP_YOHSC_H3K27me3_Created_R1.dedup.uniq.bam_peaks.IDEAS.bed"
]

file_array = [h3k4me3_peaks, h3k27ac_peaks, h3k27me3_peaks_stringent, h3k27me3_peaks_relaxed]
peak_cnt = 0
for ele in file_array:
	file_name_yhsc = ele[0]
	file_name_ohsc = ele[1]
	file_name_merged = ele[2]
	file_name_intersect = ele[3]
	file_name_created = ele[4]
	f = open(file_name_intersect)
	fw = open(file_name_created,"w")
	cnt = 0
	x = {}
	y = {}
	z = {}
	for line in f:
		line = line.strip()
		data = line.split("\t")
		the_key1 = data[0] + ":" + data[1] + "-" + data[2]
		if len(data) == 6:
			the_key2 = data[3] + ":" + data[4] + "-" + data[5]
			start = int(data[4])
			end = int(data[5])
		else:
			the_key2 = data[10] + ":" + data[11] + "-" + data[12]
			start = int(data[11])
			end = int(data[12])
		if int(data[1]) < start:
			start = int(data[1])
		if int(data[2]) > end:
			end = int(data[2])
		if peak_cnt == 3:
			fw.write(data[0] + "\t" + str(start) + "\t" + str(end-1) + "\t" + data[0] + ":" + str(start) + "-" + str(end) + "\tYOHSC\t.\n")	
		else:
			fw.write(data[0] + "\t" + str(start) + "\t" + str(end) + "\t" + data[0] + ":" + str(start) + "-" + str(end) + "\tYOHSC\t.\n")	
		x.update({the_key1:the_key2})
		y.update({the_key2:the_key1})
		for i in range(int(data[1]) + 20, int(data[2]) + 20, 20):
			the_key_only = data[0] + ":" + str((i/20)*20)
			z.update({the_key_only:1})
		if len(data) == 6:
			for i in range(int(data[4]), int(data[5]) + 20, 20):
				the_key_only = data[0] + ":" + str((i/20)*20)
				z.update({the_key_only:1})
		else:
			for i in range(int(data[11]), int(data[12]) + 20, 20):
				the_key_only = data[0] + ":" + str((i/20)*20)
				z.update({the_key_only:1})
		cnt += 1
	f.close()
	print cnt, file_name_intersect, len(x), len(y)

	f = open(file_name_yhsc)
	cnt = 0
	x_only = {}
	for line in f:
		line = line.strip()
		data = line.split("\t")
		the_key = data[0] + ":" + data[1] + "-" + data[2]
		if x.has_key(the_key):
			pass
		else:
			start = int(data[1])
			end = int(data[2])
			if peak_cnt == 3:
				fw.write(data[0] + "\t" + data[1] + "\t" + str(end-1) + "\t" + data[0] + ":" + data[1] + "-" + data[2] + "\tYHSC\t.\n")
			else:
				fw.write(data[0] + "\t" + data[1] + "\t" + str(end) + "\t" + data[0] + ":" + data[1] + "-" + data[2] + "\tYHSC\t.\n")
			for i in range(int(data[1]), int(data[2]), 10):
				the_key_only = data[0] + ":" + str((i/10)*10)
				x_only.update({the_key_only:1})
		cnt += 1
	f.close()
	print cnt, file_name_yhsc
	
	f = open(file_name_ohsc)
	cnt = 0
	y_only = {}
	for line in f:
		line = line.strip()
		data = line.split("\t")
		the_key = data[0] + ":" + data[1] + "-" + data[2]
		if y.has_key(the_key):
			pass
		else:
			start = int(data[1])
			end = int(data[2])
			if peak_cnt == 3:
				fw.write(data[0] + "\t" + data[1] + "\t" + str(end-1) + "\t" + data[0] + ":" + data[1] + "-" + data[2] + "\tOHSC\t.\n")
			else:
				fw.write(data[0] + "\t" + data[1] + "\t" + str(end) + "\t" + data[0] + ":" + data[1] + "-" + data[2] + "\tOHSC\t.\n")

			for i in range(int(data[1]), int(data[2]), 20):
				the_key_only = data[0] + ":" + str((i/20)*20)
				y_only.update({the_key_only:1})
		cnt += 1
	f.close()
	print cnt, file_name_ohsc
	fw.close()
	os.system("bedtools sort -i "+file_name_created+" > tmp.bed")
	if peak_cnt > 1:
		os.system("bedtools merge -i tmp.bed > tmp2.bed")
		#os.system("bedtools merge -d 4000 -i tmp.bed > tmp2.bed")
	else:
		os.system("bedtools merge -i tmp.bed > tmp2.bed")
	f = open("tmp2.bed")
	fw = open(file_name_created,"w")
	for line in f:
		line = line.strip()
		data = line.split("\t")
		flag = 0
		for i in range(int(data[1]), int(data[2]), 20):
			the_key_only = data[0] + ":" + str((i/20)*20)
			if z.has_key(the_key_only):
				flag = 0
				break
			elif x_only.has_key(the_key_only):
				flag = 1
			elif y_only.has_key(the_key_only):
				flag = 2
		if flag == 0:
			fw.write(data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[0] + ":" + data[1] + "-" + data[2] + "\tYOHSC\n")
		elif flag == 1:
			fw.write(data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[0] + ":" + data[1] + "-" + data[2] + "\tYHSC\n")
		elif flag == 2:
			fw.write(data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[0] + ":" + data[1] + "-" + data[2] + "\tOHSC\n")
	fw.close()
	f.close()
	os.system("rm -rf tmp*.bed")
	peak_cnt += 1
