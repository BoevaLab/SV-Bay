import yaml
import optparse
import sys
import logging
import bisect
from joblib import Parallel, delayed
from bayesian import *
from bayesianinput import *
import cluster
import utils
import resource

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
	candidates = [sublink for sublink in all_sub_links_chr if sublink.num_elements >= curr_num or sublink.gamma_alelles >= 1]
	# Candidates sorted by safe_start and safe starts separately (for bisect)
	candidates_sort_beg = sorted(candidates, key = lambda sl: sl.safe_start)
	begins = [sub_link.safe_start for sub_link in candidates_sort_beg]
	# Candidates sorted by safe_end and safe ends separately (for bisect)
	candidates_sort_end = sorted(candidates, key = lambda sl: sl.safe_end)
	ends = [sub_link.safe_end for sub_link in candidates_sort_end]

	current_sub_links.sort(key = lambda sl: sl.safe_start)	
	
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
def SetFlankingRegions(current_links):
	for l in current_links:
		logger.debug('Setting flanking regions for link ' + l.name)
		l.flanking_regions = cluster.FlankingRegions([], [], [], [])
		l.flanking_regions.A1 = CreateFR('A1', l.sublink1.left_neighbor_end, l.sublink1.safe_start)
		l.flanking_regions.A2 = CreateFR('A2', l.sublink1.safe_end, l.sublink1.right_neighbor_begin)
		l.flanking_regions.B1 = CreateFR('B1', l.sublink2.left_neighbor_end, l.sublink2.safe_start)
		l.flanking_regions.B2 = CreateFR('B2', l.sublink2.safe_end, l.sublink2.right_neighbor_begin)

# Process intra-chromosomal clusters for one chromosome
def ProcessChromosome(input_data, curr_num, chrom, all_in_mem):
	logger.info('Processing links for chromosome ' + chrom)
	# Set flanking regions to relevant links
	links_to_process = [l for l in input_data.links if \
		l.chr1 == chrom and l.chr2 == chrom and l.num_elements == curr_num]
	SetFlankingRegions(links_to_process)
	# Run bayesian models
	BayesianModels(input_data, links_to_process, chrom, chrom, config, stats)

# Process translocations for chromsome pair
def ProcessTranslocations(input_data, curr_num, chrom1, chrom2, all_in_mem):
	logger.info('Processing translocation links for chromosome pair ' + chrom1 + '-' + chrom2)
	# Set flanking regions to relevant links
	tr_links_to_process = [l for l in input_data.links if \
		l.chr1 == chrom1 and l.chr2 == chrom2 and l.num_elements == curr_num]
	SetFlankingRegions(tr_links_to_process)
	# Neighbors search was already performed when processing sepaparate chromosomes
	# So it's not needed
	# Run bayesian models
	BayesianModels(input_data, tr_links_to_process, chrom1, chrom2, config, stats)

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

logger.info('Done loading data, memory consumed: ' + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

# Procceed for all sizes (number of links) from the smallest to the biggest
# Skip clusters with very small number of elements (< 4) as noise
# TODO move this to config
for curr_num in input_data.numb_elem[4:]:
	logger.info('====================================================')
	logger.info('Processing sub links with number of elements ' + str(curr_num))
	logger.info('====================================================')
	
	# Find neifghbors for all sublinks for each chromosome
	# And process intra-chr sublimks
	logger.info('Searching neighbors and processing intra-chr sublinks for each chrom')
	for chrom in chromosomes:
		# Load chromosome fa and gem in input data
		if not all_in_mem:
			input_data.LoadChrom(config, chrom)

		# Find neighbors and process
		# First fill in arrays of all sublinks for this chr and 
		# and sublinks for this chr with current elements number
		all_sublinks_chr = []
		curr_sublinks_chr = []
		for l in input_data.links:
			if l.chr1 == chrom:
				all_sublinks_chr.append(l.sublink1)
				if l.num_elements == curr_num:
					curr_sublinks_chr.append(l.sublink1)
			if l.chr2 == chrom:
				all_sublinks_chr.append(l.sublink2)
				if l.num_elements == curr_num:
					curr_sublinks_chr.append(l.sublink2)

		logger.info('Number of sub links for chr ' + chrom + ': ' + str(len(all_sublinks_chr)))
		logger.info('Number of sub links for chr ' + chrom + ' with ' + str(curr_num) + \
			' elements: ' + str(len(curr_sublinks_chr)))

		safe_chr_len = min(len(input_data.chr_line_dict[chrom]), len(input_data.chr_gem_dict[chrom])) - 1
		FindNeighbors(curr_sublinks_chr, all_sublinks_chr, \
			input_data.chr_centrom_dict[chrom], chrom, safe_chr_len, curr_num)
		ProcessChromosome(input_data, curr_num, chrom, all_in_mem)
		# Unload fa and gem 
		if not all_in_mem:
			input_data.UnloadChrom(chrom)
	logger.info('Finished neighbors and intra-chr sublinks processing')
		
	logger.info('Proceessing tranlocation sublinks...')
	chrom_pairs = [(l.chr1, l.chr2) for l in input_data.links if \
		l.chr1 != l.chr2 and l.chr1 in chromosomes and l.chr2 in chromosomes]
	uniq_chrom_pairs = set(chrom_pairs)
	for (chrom1, chrom2) in uniq_chrom_pairs:
		# Load chromosomes fa and gem in input data
		if not all_in_mem:
			input_data.LoadChrom(config, chrom1)
			input_data.LoadChrom(config, chrom2)
		# Process translocations
		ProcessTranslocations(input_data, curr_num, chrom1, chrom2, all_in_mem)
		# Unload fa and gem
		if not all_in_mem:
			input_data.UnloadChrom(chrom1)
			input_data.UnloadChrom(chrom2)
	logger.info('Finished proceessing tranlocation sublinks')
	
