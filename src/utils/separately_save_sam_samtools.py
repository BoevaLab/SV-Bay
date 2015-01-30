import os
import optparse
chromosomes =[]
for i in range(1,23):
	chromosomes.append('chr'+str(i))
chromosomes.append('chrY')
chromosomes.append('chrX')
#chromosomes =['chr1']
parser = optparse.OptionParser()
parser.add_option('-i', '--input_file', dest='input_file', help='Path to the initial bam/sam file')
parser.add_option('-o', '--out_put_folder', dest='out_dir', help='Path to output directory were sam files by chromosome would be created')
#
(options, args) = parser.parse_args()

if options.input_file[-3:]=='bam':
	print 'bam'
	fat = os.popen("samtools view -h "+options.input_file)
elif options.input_file[-3:]=='sam':
	print 'sam'
	fat=open(options.input_file,'r')
else:
	print 'Input files could be sam/bam format only'
fchr = dict()

for chrom in chromosomes:
	fchr[chrom] = open(options.out_dir+chrom + '_ga.sam','w')
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

