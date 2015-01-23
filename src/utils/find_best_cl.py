import glob
import logging
import yaml
import optparse
from cluster import *
from sublink import *
from numpy import arange,array,ones,linalg
from pylab import plot,show
from pylab import *
from operator import is_not
from functools import partial
f = open('/Users/dasha/PhD/data/best_cl_chr15_names.txt','r')
new = open('/Users/dasha/PhD/data/best_cl_all_info_chr15.txt','w')
# Parse command line arguments to get config file name
parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

# Load configuration
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()
names = f.readlines()
names = [i.split()[0] for i in names]
(chr1_prev,chr2_prev,direction_prev) = ('','','')
for name in names:
	name_sp = name.split('_')
	(chr1,chr2,direction) = (name_sp[1],name_sp[2],name_sp[-1])
	#if chr1!=chr1_prev or chr2!=chr2_prev or direction!=direction_prev:
	cl = open(config['working_dir']+config['clusters_files_dir']+'all_info/'+'clusters_all_info_'+chr1+'_'+chr2+'_'+direction+'.txt','r')
	fl=0
	for line in cl:
		if name in line:
			new.write(line)
			fl=1
		elif fl and 'cl' not in line:
			new.write(line)
		elif fl and 'cl' in line:
			cl.close()
			break
new.close()
f.close()

