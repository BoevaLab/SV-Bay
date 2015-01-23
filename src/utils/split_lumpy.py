f =open('/Users/dasha/Downloads/sim400_PairedEnds.pesr.exclude.bedpe','r')
new = open('lumpy_res.txt','w')
chrom_names=[]
for i in range(1,22,1):
	chrom_names.append('chr'+str(i))

chrom_names.append('chrY')
chrom_names.append('chrX')
for line in f:
	line1 = line.split()
	if line1[0] in chrom_names and line1[3] in chrom_names:
		new.write(line)
new.close()
f.close()