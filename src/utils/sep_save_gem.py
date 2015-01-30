
import os
import optparse
parser = optparse.OptionParser()
parser.add_option('-i', '--input_file', dest='input_file', help='Path to the initial bam/sam file')
parser.add_option('-o', '--out_put_folder', dest='out_dir', help='Path to output directory were sam files by chromosome would be created')

(options, args) = parser.parse_args()
fat = open(options.input_file)
path = options.out_dir
print options.input_file
print path
flag=0	
for line in fat:
	if 'chr' in line:
		print 'Now for this chromosome ',line[1:-1]
		if flag:
			curr_file.close()
		curr_file = open(path+line[1:-1]+'_gem.txt','w')
		flag = 1
	elif flag:
		curr_file.write(line)
curr_file.close()