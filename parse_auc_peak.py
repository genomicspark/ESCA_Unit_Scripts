import pyBigWig
import sys, os

fn = sys.argv[1]
f = open(fn)
x = {} # +/- 1kb from TSS

for line in f:
	line = line.strip()
	data = line.split("\t")
	chrom = data[0]
	# +/- from TSS
	start = int(data[1])
	end = int(data[2])
	try:
		x.update({data[3]:[chrom,start,end,str(end-start)]})
	except:
		the_key = data[0] + ":" + data[1] + "-" + data[2]
		x.update({the_key:[chrom,start,end,str(end-start)]})
f.close()

print(len(x))

fn = sys.argv[2]
f = open("file_list_"+fn+".txt","w")
f.write(fn+"\n")
f.close()

f = open("file_list_"+fn+".txt")
for fn in f:
	fn = fn.strip()
	print(fn)
	bw = pyBigWig.open(fn)
	fw1 = open(fn+"_auc_peak.txt","w")
	auc_list = [x]
	fw_list = [fw1]
	file_names = ["peak"]
	cnt = 0
	for ele in auc_list:
		print(cnt, file_names[cnt])
		for the_key in ele:
			# AUC version #1
			coor = ele[the_key][0] + ":" + str(ele[the_key][1]) + "-" + str(ele[the_key][2])
			try:
				tmp = bw.values(ele[the_key][0], int(ele[the_key][1]), int(ele[the_key][2]))
				auc = 0.0
				for ele2 in tmp:
					if str(ele2) != "nan":
						auc += float(ele2)
			except:
				auc = 0.0
			try:
				the_value = the_key.split("|")[0] + "\t" + the_key.split("|")[1]
			except:
				the_value = the_key.split("|")[0] + "\t."
			try:
				read_per_bp = float(auc)/float(ele[the_key][3])
				fw_list[cnt].write(the_value+"\t"+ele[the_key][3]+"\t"+coor+"\t"+str(auc)+"\t"+str(read_per_bp)+"\n")
			except:
				read_per_bp = 0.0
				fw_list[cnt].write(the_value+"\t"+ele[the_key][3]+"\t"+coor+"\t"+str(auc)+"\t"+str(read_per_bp)+"\n")
		cnt += 1
	fw1.close()
f.close()
