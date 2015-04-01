import os
import optparse

##########################################################
# Split whole genome gem file to per-chromosome files
##########################################################

parser = optparse.OptionParser()
parser.add_option('-i', '--input_file', dest = 'input_file', help = 'gem file to split')
parser.add_option('-o', '--output_dir', dest = 'output_dir', help = 'output directory were per-chromosome gem files will be created')

(options, args) = parser.parse_args()

print 'Input file:', options.input_file
print 'Output dir:', options.output_dir

fin = open(options.input_file)

flag = 0
for line in fin:
	if 'chr' in line:
		print 'Processing chromosome', line[1:-1]
		if flag:
			curr_file.close()
		curr_file = open(options.output_dir + line[1:-1] + '_gem.txt','w')
		flag = 1
	elif flag:
		curr_file.write(line)

curr_file.close()
fin.close()