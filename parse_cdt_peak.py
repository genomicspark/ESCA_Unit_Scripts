import operator
import pyBigWig
import sys, os
import random

f = open("DAR__YHSC-OHSC_up.name.bed")
x = {}

heatmap_size = 5000
bin_size = 10
for line in f:
	line = line.strip()
	data = line.split("\t")
	chrom = data[0]
	# +/- from TSS
	mid_point = int(data[1]) + int(data[2])
	mid_point = mid_point / 2
	start = int(mid_point) - heatmap_size
	end = int(mid_point) + heatmap_size
	the_key = chrom+":"+str(start)+"-"+str(end)
	x.update({the_key:[chrom,start,end,str(end-start),chrom+":"+str(start)+"-"+str(end)]})
f.close()

print(len(x))

fn = sys.argv[1]
f = open("file_list_"+fn+".txt","w")
f.write(fn+"\n")
f.close()

f = open("file_list_"+fn+".txt")
for fn in f:
	fn = fn.strip()
	print(fn)
	bw = pyBigWig.open(fn)
	fw1 = open(fn+"_auc_dar_up.cdt","w")
	auc_list = [x]
	fw_list = [fw1]
	file_names = ["dar_up"]
	cnt = 0
	for ele in auc_list:
		print(cnt, file_names[cnt])
		x_auc = {}
		x_val = {}
		line_cnt = 0
		for the_key in ele:
			# AUC version #1
			coor = ele[the_key][0] + ":" + str(ele[the_key][1]) + "-" + str(ele[the_key][2])
			try:
				tmp = bw.values(ele[the_key][0], int(ele[the_key][1]), int(ele[the_key][2]))
				auc = 0.0
				line_cnt += 1
				if line_cnt % 100 == 0:
					print(fn,line_cnt,the_key,ele[the_key][4],len(tmp))
				ele3 = []
				for ele2 in tmp:
					if str(ele2) != "nan":
						ele2 = ele2/40
						auc += float(ele2)
						ele3.append(str(ele2).split(".")[0] + "." + str(ele2).split(".")[1][:2]) 
					else:
						ele3.append("0.00")
				if ele[the_key][4] == "-":
					ele3.reverse()
					ele4 = the_key+"\t1"
					for element in ele3:
						ele4 += "\t" + element
				else:
					ele4 = the_key+"\t1"
					for element in ele3:
						ele4 += "\t" + element
				if bin_size == 1:
					x_auc.update({the_key:ele4})
				else:
					ele5 = the_key+"\t1"
					data = ele4.split("\t")
					bin_value = 0.00
					cnt2 = 0
					for i in range(2, len(data)-2):
						cnt2 += 1
						bin_value += float(data[i])
						if cnt2 % bin_size == 0:
							ele5 += "\t" + str(bin_value)
							bin_value = 0.0
					x_auc.update({the_key:ele5})
				x_val.update({the_key:auc})
			except:
				pass
		sorted_x = sorted(x_val.items(), key=operator.itemgetter(1), reverse=True)
		header = ""
		for i in range(-1*heatmap_size+bin_size,heatmap_size,bin_size):
			header += "\t" + str(i)
		fw_list[cnt].write("Uniqe\tWeight"+header+"\n")
		for the_key, the_value in sorted_x:
			fw_list[cnt].write(x_auc[the_key] + "\n")
		cnt += 1
	fw1.close()
f.close()
