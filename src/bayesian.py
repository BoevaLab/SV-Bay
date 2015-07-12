import bisect
import logging
import scipy.stats
import scipy.misc
import math
import sys

logger = logging.getLogger('main_logger')

# Get chrom and gem parts for abnormal region
def AbnRegion(link, read_length, median, biggest_normal, chrom1_line, chrom2_line, gem_line_chr1, gem_line_chr2, exp_direction):
	shift = biggest_normal
	if 	link.rightmost_begin + read_length - link.begin > shift or \
		link.end - link.leftmost_end + read_length > shift:
		shift = 2 * biggest_normal

	if link.direction_type[0] == exp_direction[0]:
		bp1_begin = link.rightmost_begin + read_length
		bp1_end = link.begin + shift
	else:
		bp1_begin = link.rightmost_begin + read_length - shift
		bp1_end = link.begin

	if link.direction_type[1] == exp_direction[1]:
		bp2_begin = link.end - shift
		bp2_end = link.leftmost_end - read_length
	else:
		bp2_begin = link.end
		bp2_end = link.leftmost_end - read_length + shift

	# Calculate actual breakpoint simply as the middle of [bp_begin, bp_end]
	bp1_mid = (bp1_begin + bp1_end) / 2
	bp2_mid = (bp2_begin + bp2_end) / 2
	bp1_size = bp1_end - bp1_begin
	bp2_size = bp2_end - bp2_begin

	abn_reg1 = [bp1_mid - shift, bp1_mid] if link.direction_type[0] == exp_direction[0] \
			else [bp1_mid, bp1_mid + shift]
	abn_reg2 = [bp2_mid, bp2_mid + shift] if link.direction_type[1] == exp_direction[1] \
			else [bp2_mid - shift, bp2_mid]

	abn_reg_chrom1 = chrom1_line[abn_reg1[0] : abn_reg1[1]]
	abn_reg_chrom2 = chrom1_line[abn_reg2[0] : abn_reg2[1]]
	abn_reg_gem1 = gem_line_chr1[abn_reg1[0] : abn_reg1[1]]
	abn_reg_gem2 = gem_line_chr1[abn_reg2[0] : abn_reg2[1]]

	return (abn_reg_chrom1, abn_reg_chrom2, abn_reg_gem1, abn_reg_gem2)

# correction is +1 or -1
# allele_set will be modified outside of the function scope,
# because sets are mutable
def AddCorrectedModels(allele_set, obs_num, exp_num, int_sol, correction):
	allele = int_sol + correction
	P = scipy.stats.distributions.poisson.logpmf(obs_num, exp_num * int(allele))
	logger.debug('Initial solution for ' + str(correction) + ' is ' + str(allele))
	logger.debug('Initial P for ' + str(correction) + ' is ' + str(P))

	count = 0
	P_prev = P
	while allele > 0 and count < 3:
		count += 1
		allele_set.append(allele)
		allele = allele + correction
		if allele < 0:
			break
		P_prev = P
		P = scipy.stats.distributions.poisson.logpmf(obs_num, exp_num * int(allele))
		if P_prev < P and count >= 2:
			break

# Create set of models knowing 
# observed and expected numbers and initial solution
def CreateSetModels(obs_num, exp_num, int_sol):
	logger.debug('Initial solution is ' + str(int_sol))
	P = scipy.stats.distributions.poisson.logpmf(obs_num, exp_num * int(int_sol))
	logger.debug('Initial P is ' + str(P))
	allele_set = [int_sol]
	AddCorrectedModels(allele_set, obs_num, exp_num, int_sol, -1)
	AddCorrectedModels(allele_set, obs_num, exp_num, int_sol, +1)
	allele_set.sort()
	logger.debug(allele_set)
	return allele_set

# Calculate number of normal fragments in region
# (observed number)
def CalcNumNormFrag(flank_reg, chrom, normal_fragments):
	normal_fragments_for_chrom = normal_fragments[chrom]
	# normal_fragments_for_chrom is sorted, so we may use binary search to determine first and last nf in flank_reg
	first_nf_in_flanking_region = bisect.bisect_left(normal_fragments_for_chrom, flank_reg[0])
	last_nf_in_flanking_region = bisect.bisect_left(normal_fragments_for_chrom, flank_reg[1])

	return max(last_nf_in_flanking_region - first_nf_in_flanking_region,0)

# Calculate exected and observed fragment numbers for region
def ExpBestObs(chrom_line,reg, gem_line, read_length, median, cummulative_length_probabilities, lambda_gen_index, flag_norm, chrom, normal_fragments, config,stats):
	c1 = reg[0] # related chromosome part begin
	c2 = min(reg[1] + median+1, len(chrom_line)-1) # related chromosome part end
	logger.debug('reg[1] + median+1 = '+str(reg[1] + median+1))
	logger.debug(' len(chrom_line) = '+str(len(chrom_line)))
	chrom_line_part = chrom_line[c1:c2]
	gem_line_part = gem_line[c1:c2]
	exp_numb = ExpNumFrag(chrom_line_part, gem_line_part, read_length, median, cummulative_length_probabilities, lambda_gen_index, flag_norm, config['ploidy'],[],[], stats)
	obs_numb = CalcNumNormFrag(reg, chrom, normal_fragments)
	if exp_numb == 0.0 or config['ploidy']*exp_numb<=5:
		best = -1
	else:
		best = int(round(obs_numb/exp_numb))
	return (exp_numb, obs_numb, best)

# Returns array of 101 elements - for each possible percent (0-100) of g/c it contains 
# number of mapable window positions
def Correction(chrom_line_part, gem_line_part, flank_reg_size, read_length, median,abn_flag,stats, abn_reg_right_gem,abn_reg_right_chrom):
	correction_mapp = [0 for i in xrange(0, 101)]
	curr_lambda = 0
	for k in xrange(median + 1):
		if chrom_line_part[k] == 'c' or chrom_line_part[k] == 'g':
			curr_lambda += 1

	test_mapp_left=[]
	if abn_flag:
		logger.debug('Hey I calculate the expnumfrag for the gamma')
		logger.debug('flank_reg_size = '+str(flank_reg_size))
		logger.debug(len(abn_reg_right_gem))
		common_gem = gem_line_part+abn_reg_right_gem
		logger.debug('Number of mapable positions in left region = '+str(gem_line_part.count('!')))
		logger.debug('Number of mapable positions in right region = '+str(abn_reg_right_gem.count('!')))
		logger.debug('length of left chrom_line = '+str(len(chrom_line_part)))
		logger.debug('length of left chrom_line = '+str(len(abn_reg_right_chrom)))
		test_mapp=[]
		chrom_line_part_abn=chrom_line_part+abn_reg_right_chrom
	if not abn_flag:
		for i in xrange(1, flank_reg_size  - median):
			if chrom_line_part[i - 1] == 'c' or chrom_line_part[i - 1] == 'g':   
				curr_lambda -= 1
			if chrom_line_part[i + median ] == 'c' or chrom_line_part[i + median ] == 'g' :
				curr_lambda += 1
			if gem_line_part[i] == '!':
				correction_mapp[(curr_lambda * 100) / median] += 1		
	else:
		logger.debug('I am in abn_reg!')
		for i in xrange(1, len(chrom_line_part)+len(abn_reg_right_chrom)  - median):
			if chrom_line_part_abn[i - 1] == 'c' or chrom_line_part_abn[i - 1] == 'g':
				curr_lambda -= 1
			if chrom_line_part_abn[i + median ] == 'c' or chrom_line_part_abn[i + median ] == 'g' :
				curr_lambda += 1
			if common_gem[i] == '!':
				mapp_prob = (common_gem[int(stats['a'])+i:int(stats['b'])+i].count('!'))/float(int(stats['b']) - int(stats['a']))
				correction_mapp[(curr_lambda * 100) / median]+=mapp_prob
	if abn_flag:
		logger.debug('curr_lambda = '+str(curr_lambda))	
		logger.debug('Here the mapp_prob would be with the right side from the break_point')
		logger.debug('b = '+str(stats['b']))
		logger.debug('a = '+str(stats['a']))
		logger.debug('mapp_prob = '+str(sum(correction_mapp)))
	if correction_mapp<0:
		logger.debug('What the fuck? correction_mapp = '+str(correction_mapp))
	return (correction_mapp)

# Calculate coefficient to multiply exp_num if abn_flag is set
# Based on length probabilities
def AbnFlagCoef(flank_reg, cummulative_length_probabilities):
	prob_sum = 0.0
	first_ge_len_idx = 0
	for i in xrange(flank_reg):
		if cummulative_length_probabilities[first_ge_len_idx][0] >= i:
			prob_sum += (1.0 - cummulative_length_probabilities[first_ge_len_idx][1])
		else:
			for k in xrange(first_ge_len_idx + 1, len(cummulative_length_probabilities)):
				if cummulative_length_probabilities[k][0] >= i:
					first_ge_len_idx = k
					prob_sum += (1.0 - cummulative_length_probabilities[k][1])
					break
	return prob_sum / (flank_reg)

# Calculate expected number for region
def ExpNumFrag(chrom_line_part, gem_line_part, read_length, median, cummulative_length_probabilities, lambda_index, abn_flag, ploidy,abn_reg_right_gem,abn_reg_right_chrom,stats):
	#logger.debug('In ExpNumFrag, flanking region: ' + str(flank_reg[0]) + ' .. ' + str(flank_reg[1]))
	logger.debug('Flanking region len: ' + str(len(chrom_line_part)))
	#logger.debug('Chromosome len: ' + str(len(chrom_line)))

	correction = []
	flank_reg_size = len(chrom_line_part)
	correction_mapp = Correction(chrom_line_part, gem_line_part, flank_reg_size, read_length, median,abn_flag,stats,abn_reg_right_gem,abn_reg_right_chrom)
	exp_num = 0
	if correction_mapp<0:
		logger.debug('Hell no! correction_mapp = '+str(correction_mapp))
	for i in xrange(101):
		exp_num += correction_mapp[i] * lambda_index[i]# Number of mapable window positions with this g/c percent * corresponding probability
	if exp_num<0:
		logger.debug(lambda_index)
		logger.debug(correction_mapp)
	if abn_flag:
		logger.debug('befor lenght change exp_num_g = '+str(exp_num))
		logger.debug('abn_flag = ' + str(abn_flag))
		exp_num *= AbnFlagCoef(flank_reg_size, cummulative_length_probabilities)
		logger.debug('exp_num_g = '+str(exp_num))
	return exp_num / ploidy

def sortMIN(mod_prob):
	return float(mod_prob[0])

# Build model set for gamma (number of lost alleles)
def BuildGammaSet(gamma_best, biggest_normal, num_abnormal, numb_all_abn, exp_num_G, chrom1_line,prob_zero):
	gamma_set = []
	if gamma_best == 0:
		p = ((biggest_normal)/float(len(chrom1_line)))
		P_G = math.log(int(scipy.misc.comb(numb_all_abn,num_abnormal, exact=1)))+math.log(p)*num_abnormal+math.log(1-p)*(numb_all_abn-num_abnormal)
	else:
		P_G = scipy.stats.distributions.poisson.logpmf(num_abnormal,exp_num_G*int(gamma_best))
	gamma_set.append(gamma_best)
	gamma = gamma_best+1
	P_G_plus = scipy.stats.distributions.poisson.logpmf(num_abnormal,exp_num_G*int(gamma))
	count = 0
	while  P_G_plus > -80 and count<2:
	#while  P_G_plus > 10**(-3) and count<4:
		count+=1
		logger.debug('P_G_plus ='+str(P_G_plus))
		logger.debug('gamma = '+str(gamma))
		gamma_set.append(gamma)
		gamma = gamma +1
		P_G_plus = scipy.stats.distributions.poisson.logpmf(num_abnormal,exp_num_G*int(gamma))
	gamma = gamma_best - 1
	P_G_minus = scipy.stats.distributions.poisson.logpmf(num_abnormal,exp_num_G*int(gamma))
	count = 0
	while  P_G_minus > -80 and count<2:
	#while  P_G_minus > 10**(-3) and count<4:
		gamma_set.append(gamma)
		gamma = gamma - 1
		logger.debug('P_G_minus ='+str(P_G_minus))
		logger.debug('gamma = '+str(gamma))
		if gamma != 0:
			P_G_minus = scipy.stats.distributions.poisson.logpmf(num_abnormal,exp_num_G*int(gamma))
		elif gamma==0:

			p = biggest_normal / float(len(chrom1_line))
			P_G_minus = math.log(int(scipy.misc.comb(numb_all_abn,num_abnormal, exact=1)))*math.log(p)*num_abnormal*math.log(1-p)*numb_all_abn
			#P_G_minus = lg(int(scipy.misc.comb(numb_all_abn,num_abnormal, exact=1))*p**num_abnormal*(1-p)**numb_all_abn
	if 0 not in gamma_set:
		logger.debug('I added 0 to gamma illegaly')
		gamma_set.append(0)
	logger.debug('gamma_best = '+str(gamma_best)+' gamma_set = '+str(gamma_set))
	return gamma_set

# Create set of models to check
def BuildModelsTestSet(gamma_set, alpha_1_set, alpha_2_set, betta_1_set, betta_2_set, exp_dir, direction_type):
	logger.debug('Time to check the BuildModelsTestSet ')
	logger.debug('Expected and actual direction: ' + exp_dir + ' ' + direction_type)

	models_test_set = []
	for gamma in gamma_set:
		pairs_betta2 = []
		pairs_alpha1 = []
		if alpha_1_set != [-1]:
			for alpha_1 in alpha_1_set:
				if direction_type[0] == exp_dir[0] and alpha_1 >= gamma:
					pairs_alpha1.append([alpha_1, alpha_1 - gamma])
				elif direction_type[0] != exp_dir[0]:
					pairs_alpha1.append([alpha_1, alpha_1 + gamma])
		elif alpha_2_set != [-1]:             
			for alpha_2 in alpha_2_set:
				if direction_type[0] == exp_dir[0] :
					pairs_alpha1.append([alpha_2 + gamma, alpha_2])
				elif direction_type[0] != exp_dir[0] and alpha_2 >= gamma:
					pairs_alpha1.append([alpha_2 - gamma, alpha_2])
		else:
			pairs_alpha1.append([-1, -1])

		if betta_2_set != [-1]:
			for betta_2 in betta_2_set:
				if direction_type[1] == exp_dir[1] and betta_2 >= gamma:
					pairs_betta2.append([betta_2 - gamma, betta_2])
				elif direction_type[1] != exp_dir[1]:
					pairs_betta2.append([betta_2 + gamma, betta_2])
		elif betta_1_set != [-1]:
			for betta_1 in betta_1_set:
				if direction_type[1] != exp_dir[1] and betta_1 >= gamma:
					pairs_betta2.append([betta_1, betta_1 - gamma])
				elif direction_type[1] == exp_dir[1] :
					pairs_betta2.append([betta_1, betta_1 + gamma])
		else:
			pairs_betta2.append([-1,-1])

		for alpha_1_pair in pairs_alpha1:
				for betta_2_pair in pairs_betta2:
					models_test_set.append([alpha_1_pair[0], alpha_1_pair[1], betta_2_pair[0], betta_2_pair[1], gamma])

	return models_test_set

# Get probability for this model (alleles number)
# for this flanking region
def ProcessOneFR(fr_name, flank_reg, model, exp_numb, obs_numb, numb_all_abn, chrom_len):
	P = 0.0
	if not flank_reg or flank_reg[1] - flank_reg[0] <= 1000:
		logger.debug('Not processing ' + fr_name)
		return P
	logger.debug('Flanking region ' + fr_name + ': ' + str(flank_reg))
	logger.debug('Flanking region len: ' + str(flank_reg[1] - flank_reg[0]))
	
	if model > 0:
		P = scipy.stats.distributions.poisson.logpmf(obs_numb, exp_numb * int(model))
		logger.debug('Observed and expected numbers: ' + str(obs_numb) + ' ' + str(exp_numb * int(model)))
		logger.debug('P: ' + str(P))
		if math.isnan(P):
			P = -sys.maxint
	elif model == 0:
		logger.debug('The model is equal to zero so we will calculate the probability of missmapings at this region')
		logger.debug('Observed number: ' + str(obs_numb))
		logger.debug('numb_all_abn/chrom_len = ' + str(float(numb_all_abn)/chrom_len))
		exp_numb_mm = float(numb_all_abn) * (flank_reg[1] - flank_reg[0])/ chrom_len
		logger.debug('exp_numb_mm: ' + str(exp_numb_mm))
		P = scipy.stats.distributions.poisson.logpmf(obs_numb, exp_numb_mm)
		logger.debug('P: ' + str(P))
		if math.isnan(P):
			P = -sys.maxint
	return P

# Get probability for this model (alleles number)
# For this cluster
def ProcessModel(model, l, flag_equal, exp_num_G,  biggest_normal, chrom1_line, numb_all_abn, \
				exp_numb_A1, exp_numb_A2, exp_numb_B1, exp_numb_B2, \
				obs_numb_A1, obs_numb_A2, obs_numb_B1, obs_numb_B2, \
				prob_zero, chrom_len_A, chrom_len_B):
	logger.debug('**** ProcessModel ' + str(model) + ' for link ' + str(l.name) + '*****\n')
	P_A1 = ProcessOneFR('A1', l.flanking_regions.A1, model[0], exp_numb_A1, obs_numb_A1, numb_all_abn, chrom_len_A)
	P_B2 = ProcessOneFR('B2', l.flanking_regions.B2, model[3], exp_numb_B2, obs_numb_B2, numb_all_abn, chrom_len_B)
	P_A2 = ProcessOneFR('A2', l.flanking_regions.A2, model[1], exp_numb_A2, obs_numb_A2, numb_all_abn, chrom_len_A)
	P_B1 = 0.0
	if not flag_equal:
		P_B1 = ProcessOneFR('B1', l.flanking_regions.B1, model[2], exp_numb_B1, obs_numb_B1, numb_all_abn, chrom_len_B)

	if model[4] != 0:
		logger.debug('Observed number of fragments in cluster of abnornal fragments = ' + str(l.num_elements))
		logger.debug('Expected number of abnormal  fragments = ' + str(exp_num_G))
		logger.debug('prob_zero = ' + str(prob_zero))
		logger.debug('log(prob_zero) = ' + str(math.log(prob_zero)))
		P_G=scipy.stats.distributions.poisson.logpmf(l.num_elements,exp_num_G*int(model[4])) + math.log(prob_zero)
	else:
		p = (biggest_normal/float(len(chrom1_line)))
		logger.debug('p = '+str(p))
		logger.debug('1 - prob_zero: ' + str(1 - prob_zero) + ', log(1 - prob_zero): ' + str(math.log(1 - prob_zero)))
		P_G = math.log(scipy.misc.comb(numb_all_abn, l.num_elements, exact=1)) + math.log(p)*l.num_elements+math.log(1-p)*(numb_all_abn-l.num_elements)+math.log(1-prob_zero)
		logger.debug('Here we will calculate the probability to observe ' + str(l.num_elements) + ' fragments randomly')

	logger.debug('For G num elements = ' + str(l.num_elements) + ' expected = ' + str(exp_num_G*int(model[4])) + ' region' )
	logger.debug('Here is the prob for P_G=' + str(P_G))
	logger.debug('P_A1 P_A2 P_B1 P_B1 P_G ' + str((P_A1, P_A2, P_B1, P_B2, P_G)))
	logger.debug('P_A1*P_A2*P_B1*P_B2*P_G = ' + str(P_A1+P_A2+P_B1+P_B2+P_G))
	return P_A1+P_A2+P_B1+P_B2+P_G

# Helper to initilaise flanking region
def FRDefaults():
	return (0, 0, -1, [-1]) 

# Initialise expected and observed number 
# and models set for flanking region
def InitFR( fr_name, flank_reg, chrom_line, gem_line, chrom, read_length, median, \
			cummulative_length_probabilities, lambda_norm_index, normal_fragments, config, stats):
	(exp_numb, obs_numb, best, result_set) = FRDefaults()
	if flank_reg and flank_reg[1] - flank_reg[0] > 1000:
		(exp_numb, obs_numb, best) = ExpBestObs(chrom_line, flank_reg, gem_line, read_length, median, \
												cummulative_length_probabilities, lambda_norm_index, 0, \
												chrom, normal_fragments, config, stats)
		result_set = CreateSetModels(obs_numb, exp_numb, best)
		logger.debug('For ' + fr_name + ' exp_numb, obs_numb, best: ' + str(obs_numb) + ' ' + str(exp_numb) + ' ' + str(best))
		# if exp_numb * config['ploidy']<10:
		#     logger.debug(fr_name + ': ' + str(flank_reg) + ' len = ' + str(flank_reg[1] - flank_reg[0]))
	else:
		logger.debug(fr_name + ' is empty')
	logger.debug(fr_name + ' set: ' + str(result_set))
	return (exp_numb, obs_numb, best, result_set)

# Main processing
# Build set of models, calculate every model probability
# and choose best for each cluster
def BayesianModels(input_data, links_to_process, chrom1, chrom2, config, stats):
	# Get data we need from input_data
	lambda_norm_index = input_data.lambda_norm_index
	lambda_abnorm_index = input_data.lambda_abnorm_index
	normal_fragments = input_data.chr_normal_fragments_dict
	cummulative_length_probabilities = input_data.cummulative_length_probabilities
	chrom_line_A = input_data.chr_line_dict[chrom1]
	chrom_line_B = input_data.chr_line_dict[chrom2]
	gem_line_A = input_data.chr_gem_dict[chrom1]
	gem_line_B = input_data.chr_gem_dict[chrom2]

	read_length = config['read_length']

	# TODO think which stats to use when processing translocations
	median = stats['per_chr_stats'][chrom1]['median']
	numb_all_abn = stats['per_chr_stats'][chrom1]['num_all_abn']
	biggest_normal = stats['per_chr_stats'][chrom1]['biggest_normal']
	expected_direction = stats['per_chr_stats'][chrom1]['flag_direction']
	D = stats['per_chr_stats'][chrom1]['R']
	out_put_prob = open(config['working_dir'] + config['links_probabilities_file'], 'a')
	output_links = open(config['working_dir'] + config['valid_links_dir'] + chrom1 + '_' + chrom2 + '_valid_links_4.txt', 'a')
 
	logger.debug('Links to process in total: ' + str(len(links_to_process)))

	for l in links_to_process:
		logger.debug('********** Creating set of models for link: ' + l.name + '**********')
		flag_equal = l.flanking_regions.B1 and l.flanking_regions.A2 and l.flanking_regions.B1 == l.flanking_regions.A2

		(exp_numb_A1, obs_numb_A1, alpha_1_best, alpha_1_set) = \
			InitFR('A1', l.flanking_regions.A1, chrom_line_A, gem_line_A, l.chr1, read_length, median, cummulative_length_probabilities, lambda_norm_index, normal_fragments, config, stats)
		(exp_numb_A2, obs_numb_A2, alpha_2_best, alpha_2_set) = \
			InitFR('A2', l.flanking_regions.A2, chrom_line_A, gem_line_A, l.chr1, read_length, median, cummulative_length_probabilities, lambda_norm_index, normal_fragments, config, stats)
		(exp_numb_B2, obs_numb_B2, betta_2_best, betta_2_set) = \
			InitFR('B2', l.flanking_regions.B2, chrom_line_B, gem_line_B, l.chr2, read_length, median, cummulative_length_probabilities, lambda_norm_index, normal_fragments, config, stats)
		if flag_equal:
			(exp_numb_B1, obs_numb_B1, betta_1_best, betta_1_set) = FRDefaults()
		else:
			(exp_numb_B1, obs_numb_B1, betta_1_best, betta_1_set) = \
				InitFR('B1', l.flanking_regions.B1, chrom_line_B, gem_line_B, l.chr2, read_length, median, cummulative_length_probabilities, lambda_norm_index, normal_fragments, config, stats)

		(abn_reg_chr,abn_reg_right_chrom,abn_reg_gem,abn_reg_right_gem) = AbnRegion(l, read_length, median, biggest_normal, chrom_line_A, chrom_line_B, gem_line_A, gem_line_B, expected_direction)
		exp_num_G = ExpNumFrag(abn_reg_chr, abn_reg_gem, read_length, median, cummulative_length_probabilities, lambda_abnorm_index, 1, config['ploidy'],abn_reg_right_gem,abn_reg_right_chrom, stats)
		if exp_num_G * config['ploidy'] < 1:
			gamma_best = 0
		else:
			gamma_best = int(round(l.num_elements / exp_num_G))
		logger.debug('exp_num_G = '+str(exp_num_G)+ ' obs_numb = '+str(l.num_elements)+' gamma_best = '+ str(gamma_best))
		numb_all_links = len(input_data.links)
		prob_zero = config['exp_num_sv'] / float(numb_all_links)
		gamma_set = BuildGammaSet(gamma_best, biggest_normal, l.num_elements, numb_all_abn, exp_num_G, chrom_line_A,numb_all_links)
		models_test_set = BuildModelsTestSet(gamma_set, alpha_1_set, alpha_2_set, betta_1_set, betta_2_set, expected_direction, l.direction_type)
		logger.debug('All the models are ' + str(models_test_set))
		logger.debug('********** I finished create the set of the models **********')

		prob_diff = []
		all_models_log =[]
		sum_prob = 0
		flag_log = 0
		for model in models_test_set:
			mult = ProcessModel(model, l, flag_equal, exp_num_G, biggest_normal, chrom_line_A, numb_all_abn, \
							exp_numb_A1, exp_numb_A2, exp_numb_B1, exp_numb_B2, \
							obs_numb_A1, obs_numb_A2, obs_numb_B1, obs_numb_B2,prob_zero,len(chrom_line_A),len(chrom_line_B))
			all_models_log.append([mult,model])
		all_models_log.sort(key=sortMIN)
		logger.debug('Before adding all_models_log[-4:]'+str(all_models_log[-4:]))
		logger.debug('all_models_log = '+str(all_models_log))
		max_pr = abs(all_models_log[-1][0])
		logger.debug('all_models_log[-1][0] = '+str(all_models_log[-1][0]))
		all_models_log = [[i[0]+max_pr, i[1]] for i in all_models_log]
		logger.debug('After adding all_models_log[-4:]'+str(all_models_log[-4:]))
		prob_diff = [[math.exp(i[0]),i[1]] for i in all_models_log if i[0]>=-745]
		for i in prob_diff:
			sum_prob+=i[0]
		logger.debug('sum_prob'+str(sum_prob))
		logger.info('Here is the sum_prob: ' + str(sum_prob))
		model_prob = []
		for diff_model in prob_diff:
			model_prob.append([diff_model[0]/sum_prob ,diff_model[-1]])
		model_prob.sort(key=sortMIN)
		logger.info('Elements in models_test_set: ' + str(len(models_test_set)))
		logger.info('Initial solution ' + str([alpha_1_best,alpha_2_best,betta_1_best,betta_2_best,gamma_best]))
		logger.info('Best model for cluster ' + l.name + ' with ' + str(l.num_elements) + ' abnormal fragmetns is ' + str(model_prob[-1]))
		out_put_prob.write('\nInitial solution for cluster'+str(l.name)+ '------ ['+ str(alpha_1_best)+';'+str(alpha_2_best)+';'+str(betta_1_best)+';'+str(betta_2_best)+';'+str(gamma_best)+']')
		out_put_prob.write('\nModels --- '+str(model_prob))
		out_put_prob.write('\nBest model for cluster ' +str(l.name)+' with '+str(l.num_elements)+' abnormal fragmetns is '+ str(model_prob[-1])+ 'with probability '+ str(model_prob[-1][0]))
		l.probability =0.0
		for mod in model_prob:
			if mod[1][4]!=0:
				l.probability+=mod[0]
		l.gamma_alelles=model_prob[-1][1][4]
		l.A1 = str(model_prob[-1][1][0])
		l.B2 = str(model_prob[-1][1][3])
		if model_prob[-1][1][4]>0:
			output_links.write(l.to_string() + '\n')
	out_put_prob.close()
	output_links.close()
