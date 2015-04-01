import yaml
import optparse
import sys
import logging
import bisect
from joblib import Parallel, delayed
from bayesian import *
from bayesianinput import *
import utils

# Parse command line arguments to get config file name
parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

# Load configuration
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()

# Check that all used directories exist
utils.check_dirs_probabilities(config)

# Load serialized stats, calculated by main_clustering.py
serialized_stats_file = open(config['working_dir'] + config['serialized_stats_file'])
stats = yaml.load(serialized_stats_file)
serialized_stats_file.close()

# Initialize logging
# Levels of logging:
# DEBUG - all messages
# INFO - basic information about execution phases etc.
# WARNING - warnings
# ERROR - errors
# CRYTICAL - crytical errors (not actually used)
logger = logging.getLogger('main_logger')
if config['debug'] == 1:
	logger.setLevel(logging.DEBUG)
else:
	logger.setLevel(logging.INFO)
#handler = logging.StreamHandler()
handler = logging.FileHandler(config['working_dir'] + config['probabilites_log_file'])
formatter = logging.Formatter('%(filename)-25sline:%(lineno)-5d%(levelname)-8s [%(asctime)s]   %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# Helper class representing flanking regions for cluster
class links_flank_reg:
	def __init__(self, name, flank_A1, flank_A2, flank_B1, flank_B2, num_abnormal, link):
		self.name = name
		self.A1 = flank_A1
		self.A2 = flank_A2
		self.B1 = flank_B1
		self.B2 = flank_B2
		self.num_abnormal = num_abnormal
		self.link = link

	def pr(self):
		print self.name, self.A1, self.A2, self.B1, self.B2, self.num_abnormal, self.link

# Perform centromer check for flanking region
def CentromCheck(flank_reg, list_centr, flag_beg_end):
	reg = flank_reg[:]
	for cent in list_centr:
		if flank_reg[1] <= cent[1] and flank_reg[0] >= cent[0]:
			logger.debug(flank_reg)
			logger.debug(cent)
			reg = [flank_reg[1],flank_reg[1]]
			break
		if flag_beg_end=='beg':
			if flank_reg[0] < cent[1] < flank_reg[1]:
				reg[0] = cent[1]
				break
		else:
			if flank_reg[0] < cent[0] < flank_reg[1]:
				reg[1] = cent[0]
				break
	return(reg)

# Find nearest neighbours for each sublink with current elements number
def FindNeighbors(current_sub_links, all_sub_links_chr, list_centr, chrom, chrom_len, curr_num):
	logger.debug('Here I write the number of current sublinks and others ' + str(len(current_sub_links))+' ' + str(len(all_sub_links_chr)))
	sl_with_max_length = max(all_sub_links_chr, key = lambda sl: sl.safe_end - sl.safe_start)
	max_sublink_lenght = sl_with_max_length.safe_end - sl_with_max_length.safe_start
	logger.debug('max_sublink_lenght = ' + str(max_sublink_lenght))

	# All candidates
	candidates = [sublink for sublink in all_sub_links_chr if sublink.link.num_elements >= curr_num or sublink.link.gamma_alelles>=1]
	# Candidates sorted by safe_start and safe starts separately (for bisect)
	candidates_sort_beg = sorted(candidates, key = lambda link: link.safe_start)
	begins = [sub_link.safe_start for sub_link in candidates_sort_beg]
	# Candidates sorted by safe_end and safe ends separately (for bisect)
	candidates_sort_end = sorted(candidates, key = lambda link: link.safe_end)
	ends = [sub_link.safe_end for sub_link in candidates_sort_end]

	current_sub_links.sort(key = lambda link: link.safe_start)	
	
	for current_sub_link in current_sub_links:
		reg = []
		logger.debug('next sublink, safe_start = ' + str(current_sub_link.safe_start) + ' safe_end = ' + str(current_sub_link.safe_end))
		# Indexes for interseption search
		left_ind = bisect.bisect_left(begins,current_sub_link.safe_start - max_sublink_lenght)
		right_ind = bisect.bisect_left(begins,current_sub_link.safe_end - max_sublink_lenght)
		# Current indexes for the current sublink
		curr_ind_beg =  bisect.bisect_left(begins,current_sub_link.safe_start)
		curr_ind_end = bisect.bisect_left(begins,current_sub_link.safe_end)
		# Indexes for the search without interseption
		curr_non_iters_ind_left = bisect.bisect_left(ends,current_sub_link.safe_start) - 1
		curr_non_iters_ind_right = bisect.bisect_left(begins,current_sub_link.safe_end)

		logger.debug('The first link according to left_ind = ' + str(candidates_sort_beg[left_ind].name))
		if curr_ind_beg<len(candidates_sort_beg):
			logger.debug(candidates_sort_beg[curr_ind_beg].name + ' ' + str(candidates_sort_beg[curr_ind_beg].safe_start))

		if curr_ind_beg == 0:
			current_sub_link.left_neighbor_end = 0 
			current_sub_link.left_neighbor_name = 'Chromosome begin'
		else:
			for other_sub_link in candidates_sort_beg[left_ind:curr_ind_beg]:
				logger.debug(other_sub_link.name + 'other_sub_link.safe_end = ' + str(other_sub_link.safe_end) + 'current_sub_link.safe_start=' + str(current_sub_link.safe_start))
				if current_sub_link == other_sub_link:
					continue
				elif current_sub_link.safe_start <= other_sub_link.safe_end: # Left intersection
					current_sub_link.left_neighbor_end = current_sub_link.safe_start
					current_sub_link.left_neighbor_name = other_sub_link.name + ' intersection'
					break
		if current_sub_link.left_neighbor_name == '':
			current_sub_link.left_neighbor_end =candidates_sort_end[curr_non_iters_ind_left].safe_end
			current_sub_link.left_neighbor_name = candidates_sort_end[curr_non_iters_ind_left].name
		if candidates_sort_end[-1] == current_sub_link:
			current_sub_link.right_neighbor_begin = chrom_len
			current_sub_link.right_neighbor_name = 'Chromosome lenght'
		else:
			for other_sub_link in candidates_sort_beg[right_ind:min(len(candidates_sort_beg),curr_ind_end)]:
				logger.debug(other_sub_link.name + ' other_sub_link.safe_start ' + str(other_sub_link.safe_start) + ' other_sub_link.safe_end ' + str(other_sub_link.safe_end))
				if other_sub_link == current_sub_link:
					continue			
				if current_sub_link.safe_end <= other_sub_link.safe_end  : # Right intersection 
					current_sub_link.right_neighbor_begin = current_sub_link.safe_end
					current_sub_link.right_neighbor_name = other_sub_link.name + ' intersection'
					break			
			if current_sub_link.right_neighbor_name == '':				
				current_sub_link.right_neighbor_begin = candidates_sort_beg[min(len(candidates_sort_beg)-1,curr_non_iters_ind_right)].safe_start
				current_sub_link.right_neighbor_name = candidates_sort_beg[min(len(candidates_sort_beg)-1,curr_non_iters_ind_right)].name		
		reg_initial = [current_sub_link.left_neighbor_end , current_sub_link.safe_start]
		
		# Centromer part
		if reg_initial[1]-reg_initial[0] > 0:
			reg = CentromCheck(reg_initial, list_centr, 'beg')
			if reg[0] != reg_initial[0]:
				current_sub_link.left_neighbor_end = reg[0]
				current_sub_link.left_neighbor_name += '_centr'

		reg_initial=[current_sub_link.safe_end , current_sub_link.right_neighbor_begin]
		if reg_initial[1]-reg_initial[0]>0:
			reg=CentromCheck(reg_initial,list_centr,'end')
			if reg[1]!=reg_initial[1]:
				current_sub_link.right_neighbor_begin = reg[1]
				current_sub_link.right_neighbor_name+='_centr'

# Helper to initailise flanking regions
def CreateFR(name, begin, end):
	if begin == end:
		logger.debug('Flanking region ' + name + ' empty')
		return []
	else:
		logger.debug('Flanking region ' + name + ': [' + str(begin) + ';' + str(end) + ']')
		return [begin, end]

# Initialise flanking regions for all current sublinks
def GetFlankingRegions(current_sub_links, numb_elem):
	link_names = [i.link.name for i in current_sub_links]
	unique_link_names = set(link_names)

	curr_fl_lnk = []
	for name in unique_link_names:
		logger.debug('Link: ' + name)
		curr_fl_lnk.append(links_flank_reg(name,[],[],[],[],numb_elem,''))
		for csl in current_sub_links:
			if csl.link.name==name:
				#curr_fl_lnk[-1].position=[csl.link.rightmost_begin,csl.link.leftmost_end]
				curr_fl_lnk[-1].link = csl.link
				logger.debug('Sublink ' + csl.name + ' [ ' + str(csl.safe_start) + '; ' + str(csl.safe_end) + ']')
				logger.debug('Neighbor in the left: ' + csl.left_neighbor_name + ' ' + str(csl.left_neighbor_end) + ', ' + \
							'neighbor in the rigth: ' + csl.right_neighbor_name + ' ' + str(csl.right_neighbor_begin))

				if csl.name[-2:]=='_1':
					curr_fl_lnk[-1].A1 = CreateFR('A1', csl.left_neighbor_end, csl.safe_start)
					curr_fl_lnk[-1].A2 = CreateFR('A2', csl.safe_end, csl.right_neighbor_begin)
				else:
					curr_fl_lnk[-1].B1 = CreateFR('B1', csl.left_neighbor_end, csl.safe_start)
					curr_fl_lnk[-1].B2 = CreateFR('B1', csl.safe_end, csl.right_neighbor_begin)
					break
	return curr_fl_lnk

# Process intra-chromosomal clusters for one chromosome
def ProcessChromosome(input_data, curr_sublinks, curr_num, chrom, all_in_mem):
	logger.info('Processing sublinks for chromosome ' + chrom)
	# Set flanking regions to relevant sublinks
	chr_sublinks = [sl for sl in curr_sublinks if sl.link.chr1 == chrom and sl.link.chr2 == chrom]
	links_flank_regions = GetFlankingRegions(chr_sublinks, curr_num)
	# Run bayesian models
	BayesianModels(input_data, links_flank_regions, chrom, chrom, config, stats)

# Process translocations for chromsome pair
def ProcessTranslocations(input_data, curr_sublinks, curr_num, chrom1, chrom2, all_in_mem):
	logger.info('Processing translocation sublinks for chromosome pair ' + chrom1 + '-' + chrom2)
	# Set flanking regions to relevant sublinks
	tr_sublinks = [sl for sl in curr_sublinks if sl.link.chr1 == chrom1 and sl.link.chr2 == chrom2]
	links_flank_regions_tr = GetFlankingRegions(tr_sublinks,curr_num)
	# Neighbors search was already performed when processing sepaparate chromosomes
	# So it's not needed
	# Run bayesian models
	BayesianModels(input_data, links_flank_regions_tr, chrom1, chrom2, config, stats)

##########################################################################################
# Main code starts here
##########################################################################################

# Get list of chromosomes to process from config
chromosomes = config['chromosomes']
logger.info('Chromosomes to process:' + str(chromosomes))
# Load all data from files
input_data = BayesianInputData(config, stats, chromosomes)
# If we have enough memory to store all chrom an gem lines, just load it now
# It will speed up processing significantly
all_in_mem = 1
if all_in_mem:
	logger.info('Loading chrom and gem lines, patience...')
	for chrom in chromosomes:
		input_data.LoadChrom(config, chrom)

# Procceed for all sizes (number of links) from the smallest to the biggest
# Skip clusters with very small number of elements (< 4) as noise
# TODO move this to config
for curr_num in input_data.numb_elem[4:]:
	logger.info('====================================================')
	logger.info('Processing sub links with number of elements ' + str(curr_num))
	logger.info('====================================================')

	curr_sublinks = [sl for sl in input_data.sub_links if sl.num_elements == curr_num]
	logger.info('Number of sub links with ' + str(curr_num) + ' elements: ' + str(len(curr_sublinks)))
	
	# Find neifghbors for all sublinks for each chromosome
	# And process intra-chr sublimks
	logger.info('Searching neighbors and processing intra-chr sublinks for each chrom')
	for chrom in chromosomes:
		# Load chromosome fa and gem in input data
		if not all_in_mem:
			input_data.LoadChrom(config, chrom)
		# Find neighbors and process
		all_sublinks_chr = [sl for sl in input_data.sub_links if sl.chrom == chrom]
		curr_sublinks_chr = [sl for sl in all_sublinks_chr if sl.num_elements == curr_num]
		FindNeighbors(curr_sublinks_chr, all_sublinks_chr, input_data.chr_centrom_dict[chrom], \
			chrom, len(input_data.chr_line_dict[chrom]), curr_num)
		ProcessChromosome(input_data, curr_sublinks, curr_num, chrom, all_in_mem)
		# Unload fa and gem 
		if not all_in_mem:
			input_data.UnloadChrom(chrom)
	logger.info('Finished neighbors and intra-chr sublinks processing')
		
	logger.info('Proceessing tranlocation sublinks...')
	chrom_pairs = [(sl.link.chr1, sl.link.chr2) for sl in curr_sublinks if \
		sl.link.chr1 != sl.link.chr2 and sl.link.chr1 in chromosomes and sl.link.chr2 in chromosomes]
	uniq_chrom_pairs = set(chrom_pairs)
	for (chrom1, chrom2) in uniq_chrom_pairs:
		# Load chromosomes fa and gem in input data
		if not all_in_mem:
			input_data.LoadChrom(config, chrom1)
			input_data.LoadChrom(config, chrom2)
		# Process translocations
		ProcessTranslocations(input_data, curr_sublinks, curr_num, chrom1, chrom2, all_in_mem)
		# Unload fa and gem
		if not all_in_mem:
			input_data.UnloadChrom(chrom1)
			input_data.UnloadChrom(chrom2)
	logger.info('Finished proceessing tranlocation sublinks')
	
