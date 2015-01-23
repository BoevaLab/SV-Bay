import yaml
import optparse
import os

##########################################################
# Removes results, written by main_probabilities.py
# Uses the same config as main scripts
##########################################################

parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

# Load configuration
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()

# Remove links probabilities
##########################################################
cmd = 'rm -f ' + config['working_dir'] + config['links_probabilities_file']
print cmd
os.system(cmd)

# Remove valid links
##########################################################
cmd = 'rm -rf ' + config['working_dir'] + config['valid_links_dir'] + '*'
print cmd
os.system(cmd)

# Remove log
##########################################################
cmd = 'rm -f ' + config['working_dir'] + config['probabilites_log_file']
print cmd
os.system(cmd)