import glob
import bisect
import optparse
import yaml


parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

def LoadGem(config, chrom):
		f = open(config['working_dir'] + config['gem_files_dir'] + chrom + '_gem.txt')
		gem_line = ''
		for line in f:
			if '~' not in line:
				if line[-1] == '\n':
					gem_line += line[:-1]
				else:
					gem_line += line
		f.close()
		return(gem_line)

def LoadNormalFragments(chrom, config):
		chr_nf_list = []
		f = open(config['working_dir'] + config['normal_fragments_dir'] + 'normal_fragments_' + chrom + '.txt')
		for line in f:
			(c, nf) = line.split()
			chr_nf_list.append(int(nf))
		f.close()
		return(chr_nf_list)

config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()

# Load serialized stats, calculated by main_clustering.py
serialized_stats_file = open(config['working_dir'] + config['serialized_stats_file'])
stats = yaml.load(serialized_stats_file)
serialized_stats_file.close()
chromosomes = config['chromosomes']

mappble_positions_all = 0
unmappable_postions_all =0
count_map_all = 0
count_unmp_pos_loc = []
all_norm_frag = 0
all_un_mapp = 0
all_mapp = 0
count_frag_all = 0 
repats = 0
####### Main process ###########
for chrom in chromosomes:
	mappble_positions_loc = 0
	count_mapp_loc = []
	count_loc = 0
	normal_fragments = []
	unmappable_postions_loc = 0
	count_frag_loc = 0
	print 'I have started to process this chromosome '+chrom
	normal_fragments = LoadNormalFragments(chrom,config)
	print 'Number of normal fragmanents in '+chrom+' = ',len(normal_fragments)
	gem_line = LoadGem(config,chrom)
	all_norm_frag+=len(normal_fragments)
	print 'I finished all loading!'
	for i in normal_fragments:
		if gem_line[i]!='!':
			if i not in count_mapp_loc:
				count_mapp_loc.append(i)
			else:
				repats+=1
			count_frag_loc+=1
	count_frag_all+=count_frag_loc
	count_map_all+=len(count_mapp_loc)
	mappble_positions_loc=gem_line.count('!')
	mappble_positions_all+=mappble_positions_loc
	unmappable_postions_loc = len(gem_line) -  mappble_positions_loc
	unmappable_postions_all += unmappable_postions_loc
	print 'Number of mappable positions in ', chrom,' is = ',mappble_positions_loc
	print 'Number of unmappable positions in ',chrom,' is = ',unmappable_postions_loc
	print 'Percentage of the normal fragmanents that were mapped on unmappable positions for chrom = ',chrom
	print 'Percentage = ', float(count_frag_loc)/len(normal_fragments)
	print 'Percentage of unmappable positions in the genome that appear to be mappable = ',len(count_mapp_loc)/float(unmappable_postions_loc)
print 'All unmappable positions = ',unmappable_postions_all
print 'All mappable positions =',mappble_positions_all
print 'Percentage of the normal fragmanents that were mapped on unmappable positions for all chromosomes'
print 'Percentage = ', float(count_frag_all)/all_norm_frag
print 'Percentage of unmappable positions in the genome that appear to be mappable for all chromosomes = ',count_map_all/float(unmappable_postions_all)
print repats, ' !!!!!!!!!!!!'




