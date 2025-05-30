# This script will generate prelimenary summary report
# The initial version didn't exclude any peaks and contains all predicted datasts
# We will follow the curation work in the next step
# The input files are the processed files from parse_auc_peak.py
# It contains the peak location and the average height per peak (normalized by peak size)
# H3K4me3, H3K27me3 => Bivalent domain detection
# H3K27ac, ATAC-seq => Enhancer detection
import operator

f = open("./preliminary_reports/file_list.txt")

# Peak dictionary
h3k4me3 = {}
h3k27me3 = {}
h3k27ac = {}
atac = {}

# Peak array
h3k4me3_files = []
h3k27me3_files = []
h3k27ac_files = []
atac_files = []
cnt = 0
for line in f:
	line = line.strip()
	data = line.split("\t")
	if cnt < 4:
		tmp = h3k4me3
		h3k4me3_files.append(line)
	elif cnt < 8:
		tmp = h3k27me3
		h3k27me3_files.append(line)
	elif cnt < 12:
		tmp = h3k27ac
		h3k27ac_files.append(line)
	else:
		tmp = atac
		atac_files.append(line)	
	fo = open("./preliminary_reports/"+line)
	for line in fo:
		line = line.strip()
		data = line.split("\t")
		if tmp.has_key(data[3]):
			tmp[data[3]] += ","+data[5]
		else:
			tmp.update({data[3]:data[5]})
	cnt += 1
	fo.close()

array1 = [h3k4me3, h3k27me3, h3k27ac, atac]
array2 = [h3k4me3_files, h3k27me3_files, h3k27ac_files, atac_files]
array3 = ["h3k4me3", "h3k27me3", "h3k27ac", "atac"]
cnt = 0
for ele in array1:
	sorted_x = sorted(array1[cnt].items(), key=operator.itemgetter(1), reverse=True)
	fw = open("./preliminary_reports/summary_"+array3[cnt]+".csv","w")
	fw2 = open("./preliminary_reports/summary_"+array3[cnt]+".bed","w")
	cnt3 = 0
	for the_key, the_value in sorted_x:
		data = the_value.split(",")
		cnt2 = 0
		tot2 = 0.0
		for i in range(0, len(data)):
			tot2 += float(data[i])
			cnt2 += 1
		if tot2/cnt2 > 0.3:
			pass
		else:
			cnt3 += 1
		fw.write(the_key + "," + the_value + "," + str(tot2/cnt2) + "\n")
		coor = the_key.split(":")
		chrom = coor[0]
		start = coor[1].split("-")[0]
		end = coor[1].split("-")[1]
		fw2.write(chrom+"\t"+start+"\t"+end+"\t"+the_key+"\n")
	print array3[cnt], cnt3
	fw2.close()
	fw.close()
	cnt += 1


