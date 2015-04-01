import yaml
import optparse
import os

##########################################################
# Removes all intermidiate files, 
# written by main_clustering.py and used by main_probabilities.py
# Uses the same config as main scripts
##########################################################

parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

# Load configuration
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()

# Remove clusters
##########################################################
cmd = 'rm -rf ' + config['working_dir'] + config['clusters_files_dir'] + '*'
print cmd
os.system(cmd)

# Remove normal fragments
##########################################################
cmd = 'rm -f ' + config['working_dir'] + config['normal_fragments_dir'] + '*'
print cmd
os.system(cmd)

# Remove length histogram
##########################################################
# cmd = 'rm -f ' + config['working_dir'] + config['length_histogram_file']
# print cmd
# os.system(cmd)

# Remove serialized stats file
##########################################################
cmd = 'rm -f ' + config['working_dir'] + config['serialized_stats_file']
print cmd
os.system(cmd)

# Remove log
##########################################################
cmd = 'rm -f ' + config['working_dir'] + config['clustering_log_file']
print cmd
os.system(cmd)
