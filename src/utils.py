import os
import sys
import cluster
import logging

logger = logging.getLogger('main_logger')

# Check if directory exists. If not:
# if crytical is True - terminate, else - warn and create
def check_dir(name, path, crytical):
	if not os.path.isdir(path):
		if crytical:
			print ('Path for ' + name + ' option (' + path + ') doesn\'t exist, terminating')
			sys.exit()
		else:
			print ('Path for ' + name + ' option (' + path + ') doesn\'t exist, creating')
			os.mkdir(path)

# Check that used directories exist
# 1) working_dir obviously must exist
# 2) sam_files_dir is input, it also must exist
# 3) normal_fragments_dir, clusters_files_dir are for output
# if they don't exist, just create them
def check_dirs_clustering(config):
	working_dir = config['working_dir']
	check_dir('working_dir', working_dir, True)

	sam_files_dir = config['working_dir'] + config['sam_files_dir']
	check_dir('sam_files_dir', sam_files_dir, True)

	normal_fragments_dir = config['working_dir'] + config['normal_fragments_dir']
	check_dir('normal_fragments_dir', normal_fragments_dir, False)

	clusters_files_dir = config['working_dir'] + config['clusters_files_dir']
	check_dir('clusters_files_dir', clusters_files_dir, False)

# Check that used directories exist
# 1) working_dir obviously must exist
# 2) normal_fragments_dir, clusters_files_dir, fa_files_dir and gem_files_dir
# are input, they also must exist
# 3) valid_links is for output, if it doesn't exist, just create it	
def check_dirs_probabilities(config):
	working_dir = config['working_dir']
	check_dir('working_dir', working_dir, True)

	normal_fragments_dir = config['working_dir'] + config['normal_fragments_dir']
	check_dir('normal_fragments_dir', normal_fragments_dir, True)

	clusters_files_dir = config['working_dir'] + config['clusters_files_dir']
	check_dir('clusters_files_dir', clusters_files_dir, True)

	fa_files_dir = config['working_dir'] + config['fa_files_dir']
	if (config['fa_files_dir'][0]=='/') or (config['fa_files_dir'][0]=='\\'):
		fa_files_dir = config['fa_files_dir']	

	check_dir('fa_files_dir', fa_files_dir, True)

	gem_files_dir = config['working_dir'] + config['gem_files_dir']
	if (config['gem_files_dir'][0]=='/') or (config['gem_files_dir'][0]=='\\'):
		gem_files_dir = config['gem_files_dir']	
	check_dir('gem_files_dir', gem_files_dir, True)

	valid_links_dir = config['working_dir'] + config['valid_links_dir']
	check_dir('valid_links_dir', valid_links_dir, False)	

# Cluster fragments of one type (fr, rf, rr, ff) for one chromosome and dump clusters to file
def clust(fragments, chrm, type, M, S, config, direction):
	logger.info('Clustering fragments for chromosome ' + chrm + ', direction type ' + type + ', input fragments: ' + str(len(fragments)) + '...')
	if fragments:
		clusters = cluster.Clustering(fragments, M, S, chrm, chrm, 0, config,direction)
	else:
		clusters = []
	logger.info('Done, clusters: ' + str(len(clusters)))
	if len(clusters) > 0:
		f = open(config['working_dir'] + config['clusters_files_dir'] + 'clusters_' + chrm + '_' + type + '.txt', 'w' )
		logger.debug('Here would be results!!')
		logger.debug(config['working_dir'] + config['clusters_files_dir'] + 'clusters_' + chrm + '_' + type + '.txt')
		for cl in clusters:
			if cl.num_elements > 1:
				f.write(cl.to_string())
				#logger.info(cl.to_string())
				f.write('\n')
		f.close()
	
# Cluster translocartions for a pair of chromosomes and dump clusters to file
def clust_translocations(fragments, chrm1, chrm2, M, S, config, direction):
	logger.info('Clustering translocations for chromosomes ' + chrm1 + ' ' + chrm2)
	if fragments:
		clusters = cluster.Clustering(fragments, M, S, chrm1, chrm2, 1, config,direction)
	else:
		clusters = []
	logger.info('Done, clusters: ' + str(len(clusters)))
	if len(clusters) > 0:
		f = open(config['working_dir'] + config['clusters_files_dir'] + 'clusters_translocations_' + chrm1 + '_' + chrm2 + '.txt', 'a')
		for cl in clusters:
			if cl.num_elements > 1:
				f.write(cl.to_string())
				f.write('\n')
		f.close()
