chromosomes =[]
for i in range(1,23):
	chromosomes.append('chr'+str(i))
chromosomes.append('chrY')
chromosomes.append('chrX')
print chromosomes

fat=open('sim400_PairedEnds_hg38.sam','r')
fchr = dict()
for chrom in chromosomes:
	fchr[chrom] = open(chrom + '.sam','w')
couple_flag=0
n = 0
k = 1
print '!!!', fchr.keys()
for line in fat:
	if line[0]=='@':
		for name in fchr.keys():
			fchr[name].write(line)
	else:
		line_sp =line.split()
		if line_sp[2] in fchr.keys():
			if line_sp[6]=='=':
				fchr[line_sp[2]].write(line)
			elif line_sp[6] in fchr.keys():
				fchr[min(line_sp[2],line_sp[6])].write(line)

	n+=1
	if n==1000000:
		k+=1
		print k*n, ' lines were processed'
		n=0
for name in fchr.keys():
	fchr[name].close()
fat.close()

