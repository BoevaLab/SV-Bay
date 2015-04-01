import os
import optparse

##########################################################
# Split whole genome bam/sam file to per-chromosome sam files
# Bam files are processed using system samtools
##########################################################

parser = optparse.OptionParser()
parser.add_option('-i', '--input_file', dest = 'input_file', help = 'bam/sam file to split')
parser.add_option('-o', '--output_dir', dest = 'output_dir', help = 'output directory were per-chromosome sam files will be created')

(options, args) = parser.parse_args()

print 'Input file:', options.input_file
print 'Output dir:', options.output_dir

# Add chromosomes chr{1-23}, chrX, chrY to a list
chromosomes = []
for i in range(1,23):
	chromosomes.append('chr' + str(i))
chromosomes.append('chrY')
chromosomes.append('chrX')

if options.input_file[-3:] == 'bam':
	fin = os.popen("samtools view -h " + options.input_file)
elif options.input_file[-3:]=='sam':
	fin = open(options.input_file, 'r')
else:
	print 'Input file should be sam/bam'

# Output file handles dictionary
fchr = dict()
for chrom in chromosomes:
	fchr[chrom] = open(options.output_dir + chrom + '.sam','w')

couple_flag = 0
n = 0
for line in fin:
	# Copy header to all output files
	if line[0] == '@':
		for name in fchr.keys():
			fchr[name].write(line)
	# And data lines only to one
	else:
		line_sp = line.split()
		if line_sp[2] in fchr.keys():
			if line_sp[6] == '=':
				fchr[line_sp[2]].write(line)
			elif line_sp[6] in fchr.keys():
				fchr[min(line_sp[2],line_sp[6])].write(line)
	n += 1
	if n % 1000000 == 0:
		print 'Lines processed:', n
		
for name in fchr.keys():
	fchr[name].close()
fin.close()

