import glob
import bisect
import optparse
import yaml

parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

def LoadChromLine(config, chrom):
	f = open(config['working_dir'] + config['fa_files_dir'] + chrom + '.fa')
	chrom_line = ''
	for line in f:
		chrom_line += line.rstrip()
	f.close()
	return(chrom_line.lower())
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

# Load configuration
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()

# Load serialized stats, calculated by main_clustering.py
serialized_stats_file = open(config['working_dir'] + config['serialized_stats_file'])
stats = yaml.load(serialized_stats_file)
serialized_stats_file.close()
chromosomes = config['chromosomes']
Fgc_Ngc = 	dict()
step = 21
median = stats['chr7']['median']
for chrom in chromosomes:
	print 'Now for this chromosome ', chrom
	gem_val = []
	gc_rate = []
	gc_gem = []
	#Fgc_Ngc[chrom] = {'Fgc':0,'Ngc':0}
	chrom_line = LoadChromLine(config,chrom)
	gem_line = LoadGem(config,chrom)
	normal_fragments = LoadNormalFragments(chrom,config)
	normal_fragments.sort()
	
	print 'Uploading process is finished'
	for i in xrange(0,len(gem_line)-median,step):
		if gem_line[i]!='!' or gem_line[i+median-70]!='!':
			gem_val.append(1)
		else:
			gem_val.append(-1)
	print 'Finished wth gem'
	for i in xrange(0,len(chrom_line)-median,step):
		gc_rate.append(round((chrom_line[i:i+median].count('g')+chrom_line[i:i+median].count('c'))/float(median),2))
	gc_gem = [gc_rate[i]*gem_val[i] for i in xrange(len(gem_val)-1)]
	print 'All windows are counted, now Fgc'
	for i in xrange(len(gc_gem)):
		if gc_gem[i]<0:
			gc_gem[i] = 'n'
	for i in set(gc_gem):
		if i not in Fgc_Ngc.keys():
			Fgc_Ngc[i] = {'Fgc':0,'Ngc':0}
		Fgc_Ngc[i]['Ngc'] += gc_gem.count(i)
	print Fgc_Ngc.keys()
	norm_ind=0
	ind_left=0
	print 'Biggest normal frag = ',normal_fragments[-1]
	print 'len chrom_line =',len(chrom_line)
	print 'lenght gc_gem = ',len(gc_gem)
	print 'lenght gem len = ',len(gem_val)
	for i in normal_fragments:
		if i>len(gem_line)-median:
			break
		ind_left+=1	
		if ind_left%50000==0:
			print 'Already processed ',ind_left,' fragments, left ',len(normal_fragments)-ind_left
		if i%step == 0:
			if gc_gem[i/step]!='n':
				Fgc_Ngc[gc_gem[i/step]]['Fgc'] +=1
			#ind = bisect.bisect_left(normal_fragments[norm_ind:],step*i)+norm_ind
			#if normal_fragments[ind]==step*i:
			#	print i
			#	Fgc_Ngc[gc_gem[i]]['Fgc'] +=1
			#	for frag in normal_fragments[norm_ind+1:]:
			#		if frag!= step*i:
			#			break
			#		else:
			#			Fgc_Ngc[gc_gem[i]]['Fgc'] +=1
			#			norm_ind+=1

lam=[]
for i in Fgc_Ngc.keys():
	print i, 'Fgc = ', Fgc_Ngc[i]['Fgc'], 'Ngc = ', Fgc_Ngc[i]['Ngc']
	lam.append([i, float(Fgc_Ngc[i]['Fgc'])/Fgc_Ngc[i]['Ngc']])













