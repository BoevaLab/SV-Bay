import yaml
import optparse
import glob
import itertools
import logging
import os
from joblib import Parallel, delayed

import fragment
import chrfragments
import utils

# Parse command line arguments to get config file name
parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

# Load configuration
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()

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
handler = logging.FileHandler(config['working_dir'] + config['clustering_log_file'])
formatter = logging.Formatter('%(filename)-20sline:%(lineno)-5d%(levelname)-8s [%(asctime)s]   %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# Check that all used directories exist
utils.check_dirs_clustering(config)

# Get list of bam and sam files in sam_files_dir
sam_files = glob.glob(config['working_dir'] + config['sam_files_dir'] + '*.sam')
sam_gz_files = glob.glob(config['working_dir'] + config['sam_files_dir'] + '*.sam.gz')
bam_files = glob.glob(config['working_dir'] + config['sam_files_dir'] + '*.bam')
all_data_files = sam_files + sam_gz_files + bam_files
# Process only files with names containing one of
# chromosomes from config
data_files_to_process = []
logger.info(all_data_files)
logger.info(config['chromosomes'])
for data_file in all_data_files:
	if any((chrom) in data_file for chrom in config['chromosomes']):
		data_files_to_process.append(data_file)

logger.info('Processing files matching chromosomes from ' + config['working_dir'] + config['sam_files_dir'])
logger.info('Files: ' + str(data_files_to_process))

# Function to process one data file
def Process(data_file, config):
	# Create ChrFragments object, load fragments from sam/bam file
	cf = chrfragments.ChrFragments(data_file, config)
	# Calculate different chromosome statistics
	cf.CollectStats()
	# Print collected stats
	cf.PrintStats()
	# Main loading (collect abnormal fragments to memory,
	# write normal to file)
	cf.MainLoad()
	# Write stats for one chromosome to temporary file
	# (not use global dict to enable parallelism)
	cf.SerializeStatsToTmp()
	# Save translocations to separate temporary files 
	# (not use global array to enable parallelism)
	cf.SaveTranslocationsToTmp()
	# Cluster abnormal fragments
	cf.PerformClustering()

# Process data files in several threads
# according to configuration
npp = 1
if 'clustering_parallel_processes' in config:
	npp = config['clustering_parallel_processes']
#Parallel(n_jobs = npp)(delayed(Process)(data_file, config) for data_file in data_files_to_process)
for data_file in data_files_to_process:
	Process(data_file, config)

# Serialize stats
stats_to_serialize = dict()
stats_to_serialize['per_chr_stats'] = dict()
for chrom in config['chromosomes']:
	stats_temp_fname = '/tmp/stats_' + chrom + '.tmp'
	stats_temp_file = open(stats_temp_fname, 'r')
	chrom_stats = yaml.load(stats_temp_file)
	stats_to_serialize['per_chr_stats'][chrom] = dict()
	stats_to_serialize['per_chr_stats'][chrom]['num_all_abn'] = chrom_stats['num_all_abn']
	stats_to_serialize['per_chr_stats'][chrom]['R'] = chrom_stats['R']
	stats_to_serialize['per_chr_stats'][chrom]['smallest_normal'] = chrom_stats['smallest_normal']
	stats_to_serialize['per_chr_stats'][chrom]['biggest_normal'] = chrom_stats['biggest_normal']
	stats_to_serialize['per_chr_stats'][chrom]['median'] = chrom_stats['median']
	stats_to_serialize['per_chr_stats'][chrom]['flag_direction'] = chrom_stats['flag_direction']
	if 'a' in chrom_stats:
		stats_to_serialize['a'] = chrom_stats['a']
	if 'b' in chrom_stats:
		stats_to_serialize['b'] = chrom_stats['b']
	if 'p' in chrom_stats:
		stats_to_serialize['p'] = chrom_stats['p']		
	stats_temp_file.close()
	os.remove(stats_temp_fname)
serialized_stats_file = open(config['working_dir'] + config['serialized_stats_file'], 'w')
serialized_stats_file.write(yaml.dump(stats_to_serialize))
serialized_stats_file.close()

# Load translocations from temporary files
all_translocations = []
for chrom in config['chromosomes']:
	trans_temp_fname = '/tmp/trans_' + chrom + '.tmp'
	trans_temp_file = open(trans_temp_fname, 'r')
	for line in trans_temp_file:
		tr = fragment.Fragment()
		tr.from_string(line)
		all_translocations.append(tr)
	trans_temp_file.close()
	os.remove(trans_temp_fname)

# Separate chromosomes processed, now process translocations
combinations_chr = itertools.permutations(config['chromosomes'], 2)
for combination_chr in combinations_chr:
	logger.info('Processing chromosome pair ' + str(combination_chr))
	c1 = combination_chr[0]
	c2 = combination_chr[1]
	# We order translocation reads while loading, c1 is always less than c2
	if c1 >= c2:
		continue
	M = stats_to_serialize['per_chr_stats'][c1]['biggest_normal']
	S = stats_to_serialize['per_chr_stats'][c1]['smallest_normal']
	exp_dir = stats_to_serialize['per_chr_stats'][c1]['flag_direction']
	
	# Split translocations with each direction
	# To cluster them separately
	trans_ff = [f for f in all_translocations if \
		f.first_read_chr == c1 and f.second_read_chr == c2 and f.direction == 'ff']
	trans_fr = [f for f in all_translocations if \
		f.first_read_chr == c1 and f.second_read_chr == c2 and f.direction == 'fr']
	trans_rf = [f for f in all_translocations if \
		f.first_read_chr == c1 and f.second_read_chr == c2 and f.direction == 'rf']
	trans_rr = [f for f in all_translocations if \
		f.first_read_chr == c1 and f.second_read_chr == c2 and f.direction == 'rr']
	
	if trans_ff:
		utils.clust_translocations(trans_ff, c1, c2, M, S, config, exp_dir)
	if trans_fr:
		utils.clust_translocations(trans_fr, c1, c2, M, S, config, exp_dir)
	if trans_rf:
		utils.clust_translocations(trans_rf, c1, c2, M, S, config, exp_dir)
	if trans_rr:
		utils.clust_translocations(trans_rr, c1, c2, M, S, config, exp_dir)

