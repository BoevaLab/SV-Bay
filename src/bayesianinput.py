import glob
import logging
from cluster import *
from sublink import *
from numpy import arange,array,ones,linalg
from pylab import plot,show
from pylab import *
from operator import is_not
from functools import partial
########################################################################
# Input data for bayesian algorithms
########################################################################

logger = logging.getLogger('main_logger')
def SortLambda(ll):
	return ll[0]
class BayesianInputData:
	###########################
	# Several dictionaries of data related to one chromosome
	# With chromosome name as a key
	###########################
	# value - list of normal fragments begins
	chr_normal_fragments_dict = dict()
	chr_normal_fragments_dict_unique = dict()
	# value - list of centromer (begin, end) pairs
	chr_centrom_dict = dict()
	# value - string with mappability info from gem file
	chr_gem_dict = dict()
	# value - chromosome string from fa file
	chr_line_dict = dict()
	# CN information about each chromosome
	chr_cnv_dict = dict()
	###########################
	# Chromosomes for which data is loaded
	chromosomes = []
	# List of lambda data pairs (gc share, probability)
	lambda_gen = []
	# Lambda gen values for possible percents of gc (0-100)
	#lambda_gen_index = []
	lambda_norm_index = []
	lambda_abnorm_index =[]
	lambda_abnorm = []
	lambda_norm =[]
	# List of pairs (fragment length, probability of this length)
	length_probabilities = []
	# List of pairs (fragment length, probability of length <= than this)
	cummulative_length_probabilities = []
	# All sub links for all chromosomes
	sub_links = []
	# Sorted list of element numbers observed in clusters
	numb_elem = []

	def __init__(self, config, stats, chromosomes):
		logger.info('Loading data...')
		self.chromosomes = chromosomes
		# Load centrom, normal fragments and sublinks for all chromosomes
		# Do not load fa and gem line for all chromosomes - it requires too much memory,
		# so we load it for separate chromosomes on demand
		#self.__LoadFreec(config)
		self.__LoadNormalFragments(config)
		self.__LoadCentrom(config)
		self.__LoadLinksAndConstructSubLinks(config, stats)
		# Load common data
		self.__LoadLambda(config, stats)
		self.__LoadLengthProbabilities(config, stats)
		logger.info('Finished loading data')

	# Load fa and gem data required to process chromosome
	def LoadChrom(self, config, chrom):
		self.__LoadGem(config, chrom)
		self.__LoadChromLine(config, chrom)
		self.__LoadFreec(config,chrom)

	def UnloadChrom(self, chrom):
		self.chr_gem_dict.pop(chrom)
		self.chr_line_dict.pop(chrom)

	# def LogLoadedDataInfo(self):
	# 	logger.info('Loaded data for chromosomes ' + str(self.chromosomes))
	# 	logger.info('Number of elements in lambda_gen: ' + str(len(self.lambda_gen)))
	# 	logger.info('First: ' + str(self.lambda_gen[0]), ', last: ' + str(self.lambda_gen[-1]))
	# 	logger.info('Number of elements in length_probabilities: ' + str(len(self.length_probabilities)))
	# 	logger.info('First: ' + str(self.length_probabilities[0]), ', last: ' + str(self.length_probabilities[-1]))
	# 	logger.info('Number of elements in cummulative_length_probabilities: ' + str(len(self.cummulative_length_probabilities)))
	# 	logger.info('First: ' + str(self.cummulative_length_probabilities[0]), ', last: ' + str(self.cummulative_length_probabilities[-1]))
	# 	logger.info('Number of elements in sub_links: ' + str(len(self.sub_links)))
	# 	logger.info('First: ' + str(self.sub_links[0]), ', last: ' + str(self.sub_links[-1]))
	# 	logger.info('Number of elements in numb_elem: ' + str(len(self.numb_elem)))
	# 	logger.info('First: ' + str(self.numb_elem[0]), ', last: ' + str(self.numb_elem[-1]))

	#########################################################################
	# Functions below should not be used directly from outside the class
	# So names are mungled
	#########################################################################
	def __LoadGem(self, config, chrom):
		#f = open(config['working_dir'] + config['gem_files_dir'] + chrom + '_hg38')
		f = open(config['working_dir'] + config['gem_files_dir'] + chrom + '_gem.txt')
		gem_line = ''
		for line in f:
			if '~' not in line:
				if line[-1] == '\n':
					gem_line += line[:-1]
				else:
					gem_line += line
		f.close()
		self.chr_gem_dict[chrom] = gem_line

	def __LoadChromLine(self, config, chrom):
		f = open(config['working_dir'] + config['fa_files_dir'] + chrom + '.fa')
		chrom_line = ''
		for line in f:
			chrom_line += line.rstrip()
		f.close()
		self.chr_line_dict[chrom] = chrom_line.lower()

	def __LoadNormalFragments(self, config):
		for chrom in self.chromosomes:
			chr_nf_list = []
			chr_nf_uniq =[]
			f = open(config['working_dir'] + config['normal_fragments_dir'] + 'normal_fragments_' + chrom + '.txt')
			for line in f:
				if '| chr |' in line:
					continue
				(c, nf,uq,mp,name) = line.split()
				if uq == '1':
					chr_nf_uniq.append(int(nf))
				if mp == '1':
					chr_nf_list.append(int(nf))
				
			f.close()
			self.chr_normal_fragments_dict[chrom] = chr_nf_list
			self.chr_normal_fragments_dict_unique[chrom] = chr_nf_uniq
	def __LoadCentrom(self, config):
		f = open(config['working_dir'] + config['centromic_file'])
		for line in f:
			(chrom, centrom_begin, centrom_end) = line[1:-2].split(',')
			chrom = chrom[1:-1] 
			if chrom in self.chromosomes:
				if chrom not in self.chr_centrom_dict:
					self.chr_centrom_dict[chrom] = []
				self.chr_centrom_dict[chrom].append((int(centrom_begin),int(centrom_end)))

	def __LoadLinksAndConstructSubLinks(self, config, stats):
		links = []
		file_names = glob.glob(config['working_dir'] + config['clusters_files_dir'] + '*.txt')
		# Process only files with names containing one of
		# chromosomes from config
		file_names_to_process = []
		for fname in file_names:
			if any(chrom in fname for chrom in config['chromosomes']):
				file_names_to_process.append(fname)

		for fname in file_names_to_process:
			logger.debug(fname)
			f = open(fname)
			for line in f:
				links.append(Cluster.from_string(line))
				if links[-1].num_elements not in self.numb_elem:
					self.numb_elem.append(links[-1].num_elements)
			f.close()
		# Sort numb elem and remove smallest numbers
		self.numb_elem = sorted(self.numb_elem)
		biggest_normal = stats['per_chr_stats'][config['chromosomes'][0]]['biggest_normal']
		initial_direction = stats['per_chr_stats'][config['chromosomes'][0]]['flag_direction'] 
		for link in links:
			self.sub_links.append(SubLink(link, 'beg', biggest_normal,initial_direction))
			self.sub_links.append(SubLink(link, 'end', biggest_normal,initial_direction))

	# TODO should be load separately for each chromosome
	def __LoadLengthProbabilities(self, config, stats):
		f = open(config['working_dir'] + config['length_histogram_file'])
		f.readline()
		for line in f:
			(l, p) = line.split()
			if float(l) >= stats['per_chr_stats'][self.chromosomes[0]]['smallest_normal']:
				self.length_probabilities.append((int(float(l)),float(p)))
			elif float(l) > stats['per_chr_stats'][self.chromosomes[0]]['biggest_normal']:
				break
		f.close()
		# Sort length probabilities by length (first tuple element will be used automaically)
		self.length_probabilities = sorted(self.length_probabilities)
		# Calculate cummulative probabilities
		# For each length l in this list we will calculate probability that length of fragment is less or equal than l
		# It is calculated as a sum of all previous probabilities, because length_probabilities is sorted
		# For the first length it will be it's own probablility, for the last one it will be 1.0
		cum_prob = 0.0
		for (l, p) in self.length_probabilities:
			cum_prob = cum_prob + p
			self.cummulative_length_probabilities.append((l, cum_prob))

	def __LoadFreec(self,config, chrom):
		logger.debug('I am here! Burn in the Hell!')
		#for chrom in self.chromosomes:
		f = open(config['working_dir'] + config['cnv_file'],'r')
		self.chr_cnv_dict[chrom]=[]
		for line in f:
			line_sp = line.split()
			if 'chr'+line_sp[0]==chrom:
				self.chr_cnv_dict[chrom].append([int(line_sp[1]),int(line_sp[2])])
		logger.debug('Freec for '+chrom)
		logger.debug(self.chr_cnv_dict[chrom])
		f.close()


	# Load lambda from file
	# If file doesn't esist, re-generate lambda
	def __LoadLambda(self, config, stats):
		try:
			flag_norm='n'
			f = open(config['working_dir'] + config['lambda_file'],'r')
			lines = f.readlines()
			logger.debug('Lambda load f.readlines()' + str(len(lines)))
			for line in lines:
				if line :
					if 'Normal lambda' in line:
						flag_norm = 1
					elif 'Abnormal lambda' in line:
						flag_norm = 0
					elif flag_norm==1:
						line = line.split()
						self.lambda_norm.append((float(line[0]), float(line[1])))
					elif flag_norm ==0:
						line = line.split()
						self.lambda_abnorm.append((float(line[0]), float(line[1])))
			f.close()

		except IOError:
			logger.warning('Lambda file not found, calculating lambda')
			(lambda_norm,lambda_abnorm) = self.__CalculateLambda(config, stats)
			self.lambda_norm= lambda_norm
			self.lambda_abnorm = lambda_abnorm

		# Form lambda_gen_index
		# For each percent store
		self.lambda_norm_index =[]
		self.lambda_abnorm_index=[]

		step = 0.01
		for i in xrange(0, 101):
			if float(i) / 100 < self.lambda_norm[0][0]:
				self.lambda_norm_index.append(self.lambda_norm[0][1])
			elif float(i) / 100 > self.lambda_norm[-1][0]:
				self.lambda_norm_index.append(self.lambda_norm[-1][1])
			else:
				self.lambda_norm_index.append(self.lambda_norm[int((round(float(i) / 100, 2) - self.lambda_norm[0][0]) / step)][1])
		for i in xrange(0, 101):
			if float(i) / 100 < self.lambda_abnorm[0][0]:
				self.lambda_abnorm_index.append(self.lambda_abnorm[0][1])
			elif float(i) / 100 > self.lambda_abnorm[-1][0]:
				self.lambda_abnorm_index.append(self.lambda_abnorm[-1][1])
			else:
				self.lambda_abnorm_index.append(self.lambda_abnorm[int((round(float(i) / 100, 2) - self.lambda_abnorm[0][0]) / step)][1])

	def __CalculateLambda(self, config, stats):
		logger.info('Started calculating lambda_gen')

		Fgc_Ngc_norm = dict()
		Fgc_Ngc_abnorm	= dict()	
		step = 61
		count = 0

		#if len(self.chromosomes)>=2:
		#	chrom_lambda = self.chromosomes[:2]
		#else:
		#	chrom_lambda = self.chromosomes
		#	step= step/2
		#print self.chromosomes[1]
		#chrom_lambda = [self.chromosomes[0]]
		for chrom in self.chromosomes:
			logger.debug('Processing chromosome ' + chrom+' for lambda calculation')
			logger.debug('count = '+str(count))
			self.LoadChrom(config, chrom)
			gem_val_normal = []
			gc_rate_normal = []
			gc_gem_normal = []
			
			gem_val_abnormal = []
			gc_rate_abnormal = []
			gc_gem_abnormal = []

			median = stats['per_chr_stats'][chrom]['median']
			chrom_line = self.chr_line_dict[chrom]
			logger.debug('chrom_line =  '+str(len(chrom_line)))
			logger.debug(self.chr_cnv_dict)
			cnv = self.chr_cnv_dict[chrom]
			cnv_line =[]
			centrom_line = []
			sum_cnv =0 
			gem_line = self.chr_gem_dict[chrom]
			for j in cnv:
				sum_cnv += j[1] - j[0] 
				cnv_line.append(xrange(j[0],j[1]+1))
			for j in self.chr_centrom_dict[chrom]:
				centrom_line.append(xrange(j[0],j[1]+1))
			logger.debug('sum_cnv = '+str(sum_cnv))
			logger.debug('len(gem_line)= '+str(len(gem_line)))	
			logger.debug('sum_cnv/len(gem_line) = '+str(sum_cnv/float(len(gem_line))))
			if sum_cnv/float(len(gem_line))>0.7:
				logger.info('This chromosome we skip, i.e. the average ploidy is not equal to excepted '+chrom)
				continue
			else:
				count+=1
			
			normal_fragments_normal = sorted(self.chr_normal_fragments_dict[chrom])
			normal_fragments_unique = sorted(self.chr_normal_fragments_dict_unique[chrom])
			logger.debug('Number of normal fragments = '+str(len(normal_fragments_normal)))
			logger.debug('len(normal_fragments_unique) = '+str(len(normal_fragments_unique)))
			len_lenght = len(self.length_probabilities)
			min_boundaries = int(stats['a'])
			max_boundaries = int(stats['b'])
			len_omega = float(max_boundaries - min_boundaries)
			logger.debug('min_boundaries ='+str(min_boundaries))
			logger.debug('max_boundaries ='+str(max_boundaries))
			logger.debug('len(gem) = '+str(len(gem_line)-max_boundaries))
			logger.debug('len(chr) = '+str(len(chrom_line) - max_boundaries))
			for i in xrange(0, len(gem_line) - max_boundaries, step):
				if i not in cnv_line and i not in centrom_line:
					if gem_line[i] == '!':
						gem_val_normal.append(1)
						abnormal_value = (gem_line[i+min_boundaries:i+max_boundaries].count('!'))/float(len_omega)
						if abnormal_value == 0.0:
							gem_val_abnormal.append(-1)
						else:
							gem_val_abnormal.append(abnormal_value)
					else:
						gem_val_normal.append(-1)
						gem_val_abnormal.append(-1)

			logger.debug('Finished with gem')
			gc_rate =[]
			for i in xrange(0, len(chrom_line) - max_boundaries, step):
				if i not in cnv_line and i not in centrom_line: 
					gc_rate.append(round((chrom_line[i : i + median].count('g') + chrom_line[i : i + median].count('c')) / float(median),2))
			#for i in xrange(10200,10250):
			#	logger.info(str(gc_rate[i]))
			#logger.info(len(gem_val_abnormal))
			#logger.info(len(gc_rate))
			gc_gem_normal = [gc_rate[i] * gem_val_normal[i] for i in xrange(len(gc_rate) - 1)]
			#gc_gem_abnormal = [gc_rate[i] * gem_val_abnormal[i] for i in xrange(len(gem_val_abnormal) - 1)]

			logger.debug('All windows are counted, now Fgc')
			for i in xrange(len(gc_gem_normal)):
				if gc_gem_normal[i] < 0:
					gc_gem_normal[i] = 'n'
			#for i in xrange(len(gc_gem_abnormal)):
			#	if gc_gem_abnormal[i]<0:
			#		gc_gem_abnormal[i] = 'n'
			
			for i in set(gc_gem_normal):
				if i not in Fgc_Ngc_norm.keys():
					Fgc_Ngc_norm[i] = {'Fgc':0,'Ngc':0}
				Fgc_Ngc_norm[i]['Ngc'] += gc_gem_normal.count(i)
			#for i in set(gc_gem_abnormal):
			#	if i not in Fgc_Ngc_abnorm.keys():
			#		Fgc_Ngc_abnorm[i] = {'Fgc':0,'Ngc':0}
			#	Fgc_Ngc_abnorm[i]['Ngc'] += gc_gem_abnormal.count(i)
			for i in xrange(len(gc_rate)):
				if gem_val_abnormal[i]>0:
					if gc_rate[i] not in Fgc_Ngc_abnorm.keys():
						Fgc_Ngc_abnorm[gc_rate[i]] = {'Fgc':0,'Ngc':0}
					Fgc_Ngc_abnorm[gc_rate[i]]['Ngc'] += 1/gem_val_abnormal[i]
			
			#logger.debug(str(Fgc_Ngc_norm.keys()))
			norm_ind = 0
			ind_left = 0
			#logger.debug('len chrom_line = ' + str(len(chrom_line)))
			logger.debug('lenght gc_gem = ' + str(len(gc_gem_normal)))
			logger.debug('lenght gem len = ' + str(len(gem_val_normal)))
			logger.debug('lenght gc_gem_normal = '+str(len(gc_gem_normal)))
			for i in normal_fragments_normal:
				if i > len(gem_line) - median:
					break
				ind_left += 1	
				if ind_left % 50000 == 0:
					logger.debug('Already processed ' + str(ind_left) + ' fragments, left ' + str(len(normal_fragments_normal) - ind_left))
				if i % step == 0:
					if gc_gem_normal[i / step] != 'n':
						Fgc_Ngc_norm[gc_gem_normal[i / step]]['Fgc'] += 1
				ind_left = 0
			for i in normal_fragments_unique:
				if i > len(chrom_line) - max_boundaries:
					break
				ind_left += 1	
				#if ind_left % 50000 == 0:
				#	logger.debug('Already processed ' + str(ind_left) + ' fragments, left ' + str(len(normal_fragments_unique) - ind_left))
				if i % step == 0:
					if gem_val_abnormal[i/step]!=-1:
						Fgc_Ngc_abnorm[gc_rate[i/step]]['Fgc'] += 1
			if count==2:
				break
		lambda_norm = []
		lambda_abnorm = []
		#logger.info(str(Fgc_Ngc_norm))
		#logger.info(str(Fgc_Ngc_abnorm))
		del Fgc_Ngc_norm['n']
		#del Fgc_Ngc_abnorm['n']
		Fgc_Ngc_arr = [Fgc_Ngc_norm,Fgc_Ngc_abnorm]
		for dictionary in Fgc_Ngc_arr:
			if 0.7 not in dictionary.keys():
				dictionary[0.7] = {'Fgc':0,'Ngc':0}
			if 0.3 not in dictionary.keys():
				dictionary[0.3] = {'Fgc':0,'Ngc':0}
			for key in dictionary.keys():
				if key>0.7:
					dictionary[0.7]['Fgc']+=dictionary[key]['Fgc']
					dictionary[0.7]['Ngc']+=dictionary[key]['Ngc']
					del dictionary[key]
				if key<0.3:
					dictionary[0.3]['Fgc']+=dictionary[key]['Fgc']
					dictionary[0.3]['Ngc']+=dictionary[key]['Ngc']
					del dictionary[key]
			
		for i in Fgc_Ngc_norm.keys():
			if Fgc_Ngc_norm[i]['Ngc']!=0:
				logger.debug(str([i, float(Fgc_Ngc_norm[i]['Fgc']) / Fgc_Ngc_norm[i]['Ngc']]))
				lambda_norm.append([i, float(Fgc_Ngc_norm[i]['Fgc']) / Fgc_Ngc_norm[i]['Ngc']])
			else:
				logger.debug(str(i)+str(float(0.0)))
				lambda_norm.append([i, float(0.0)])
		logger.debug('Abnormal lambda')
		for i in Fgc_Ngc_abnorm.keys():
			if Fgc_Ngc_abnorm[i]['Ngc']!=0:
				logger.debug(str([i, float(Fgc_Ngc_abnorm[i]['Fgc']) / Fgc_Ngc_abnorm[i]['Ngc']]))
				lambda_abnorm.append([i, float(Fgc_Ngc_abnorm[i]['Fgc']) / Fgc_Ngc_abnorm[i]['Ngc']])
			else:
				logger.debug(str(i)+str(float(0.0)))
				lambda_abnorm.append([i, float(0.0)])
			
		lambda_norm = filter(partial(is_not, None), lambda_norm)
		lambda_norm.sort(key = SortLambda)
		lambda_abnorm = filter(partial(is_not, None), lambda_abnorm)
		lambda_abnorm.sort(key = SortLambda)

		
		logger.debug('lambda normal, len = '+str(len(lambda_norm)))
		for i in Fgc_Ngc_norm.keys():
			logger.debug(str(i) + ' Fgc = ' + str(Fgc_Ngc_norm[i]['Fgc']) + ' Ngc = ' + str(Fgc_Ngc_norm[i]['Ngc']))
		logger.debug('lambda abnormal, len = '+str(len(lambda_abnorm)))
		for i in Fgc_Ngc_abnorm.keys():
			logger.debug(str(i)+ ' Fgc = '+str(Fgc_Ngc_abnorm[i]['Fgc'])+ ' Ngc = '+str(Fgc_Ngc_abnorm[i]['Ngc']))
		
		#UnloadChrom(config, chrom)
		f = open(config['working_dir']+'lambda_wt_FR.txt','w')
		f.write('Normal lambda '+'\n')
		for i in lambda_norm:
			f.write(str(i[0])+' '+str(i[1]))
			f.write('\n')
		f.write('Abnormal lambda '+'\n')
		for i in lambda_abnorm:
			f.write(str(i[0])+' '+str(i[1]))
			f.write('\n')
		f.close()
		logger.info('Finished calculating lambda')

		return (lambda_norm,lambda_abnorm)
