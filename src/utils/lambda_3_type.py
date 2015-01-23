import glob
import bisect
import optparse
import yaml
from numpy import arange,array,ones,linalg
from pylab import plot,show
from pylab import *
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
Fgc_Ngc_mapp = 	dict()
Fgc_Ngc_unmapp = dict()
Fgc_Ngc_half_unmapp= dict()
step = 101
median = stats['chr3']['median']
for chrom in chromosomes:
	print 'Now for this chromosome ', chrom
	gem_val_mapp = []
	gem_val_unmapp =[]
	gem_val_half_unmapp=[]
	gc_gem_unmapp_half=[]
	gc_rate = []
	gc_gem = []
	#Fgc_Ngc[chrom] = {'Fgc':0,'Ngc':0}
	chrom_line = LoadChromLine(config,chrom)
	gem_line = LoadGem(config,chrom)
	normal_fragments = LoadNormalFragments(chrom,config)
	normal_fragments.sort()
	
	print 'Uploading process is finished'
	for i in xrange(0,len(gem_line)-median,step):
		if gem_line[i]=='!' and gem_line[i+median-70]=='!':
			gem_val_mapp.append(1)
			gem_val_unmapp.append(-1)
			gem_val_half_unmapp.append(1)
		else:
			gem_val_mapp.append(-1)
			if gem_line[i] != '!' and gem_line[i+ median - 70] != '!':
				gem_val_unmapp.append(1)
				gem_val_half_unmapp.append(-1)
			else:
				gem_val_half_unmapp.append(1)
				gem_val_unmapp.append(-1)
	print 'Finished wth gem'
	for i in xrange(0,len(chrom_line)-median,step):
		gc_rate.append(round((chrom_line[i:i+median].count('g')+chrom_line[i:i+median].count('c'))/float(median),2))
	gc_gem_mapp = [gc_rate[i]*gem_val_mapp[i] for i in xrange(len(gem_val_mapp)-1)]
	gc_gem_unmapp = [gc_rate[i]*gem_val_unmapp[i] for i in xrange(len(gem_val_unmapp)-1)]
	gc_gem_unmapp_half = [gc_rate[i]*gem_val_half_unmapp[i] for i in xrange(len(gem_val_half_unmapp)-1)]
	print 'all length = ', len(gc_gem_mapp),
	print 'All windows are counted, now Fgc'
	for i in xrange(len(gc_gem_mapp)):
		if gc_gem_mapp[i]<0:
			gc_gem_mapp[i] = 'n'
	for i in xrange(len(gc_gem_unmapp)):
		if gc_gem_unmapp[i]<0:
			gc_gem_unmapp[i] = 'n'
	for i in xrange(len(gc_gem_unmapp_half)):
		if gc_gem_unmapp_half[i]<0:
			gc_gem_unmapp_half[i] = 'n'

	for i in set(gc_gem_mapp):
		if i not in Fgc_Ngc_mapp.keys():
			Fgc_Ngc_mapp[i] = {'Fgc':0,'Ngc':0}
		Fgc_Ngc_mapp[i]['Ngc'] += gc_gem_mapp.count(i)
	for i in set(gc_gem_unmapp):
		if i not in Fgc_Ngc_unmapp.keys():
			Fgc_Ngc_unmapp[i] = {'Fgc':0,'Ngc':0}
		Fgc_Ngc_unmapp[i]['Ngc'] += gc_gem_unmapp.count(i)
	for i in set(gc_gem_unmapp_half):
		if i not in Fgc_Ngc_half_unmapp.keys():
			Fgc_Ngc_half_unmapp[i] = {'Fgc':0,'Ngc':0}
		Fgc_Ngc_half_unmapp[i]['Ngc'] += gc_gem_unmapp_half.count(i)
	print Fgc_Ngc_unmapp.keys()
	print Fgc_Ngc_mapp.keys()
	print Fgc_Ngc_half_unmapp.keys()
	norm_ind=0
	ind_left=0
	print 'Biggest normal frag = ',normal_fragments[-1]
	print 'len chrom_line =',len(chrom_line)
	print 'lenght gc_gem_mapp = ',len(gc_gem_mapp)
	print 'lenght gem len mapp = ',len(gem_val_mapp)
	print 'lenght gc_gem_unmapp = ',len(gc_gem_unmapp)
	print 'lenght gem len unmapp = ',len(gem_val_unmapp)
	for i in normal_fragments:
		if i>len(gem_line)-median:
			break
		ind_left+=1	
		if ind_left%50000==0:
			print 'Already processed ',ind_left,' fragments, left ',len(normal_fragments)-ind_left
		if i%step == 0:
			if gc_gem_mapp[i/step]!='n':
				Fgc_Ngc_mapp[gc_gem_mapp[i/step]]['Fgc'] +=1
			elif gc_gem_unmapp[i / step] != 'n':
				Fgc_Ngc_unmapp[gc_gem_unmapp[i/step]]['Fgc'] +=1
			else:
				Fgc_Ngc_half_unmapp[gc_gem_unmapp_half[i/step]]['Fgc'] +=1
lam_mapp = []
lam_unmapp = []
lam_unmapp_half =[]
print 'Mappable lambda'
print 'rate 	Fgc 	Ngc'
for i in Fgc_Ngc_mapp.keys():
	print i,  Fgc_Ngc_mapp[i]['Fgc'], Fgc_Ngc_mapp[i]['Ngc']
	lam_mapp.append([i, float(Fgc_Ngc_mapp[i]['Fgc'])/Fgc_Ngc_mapp[i]['Ngc']])
print 'Unmappable lambda'
print 'rate 	Fgc 	Ngc'
for i in Fgc_Ngc_unmapp.keys():
	print i, Fgc_Ngc_unmapp[i]['Fgc'], Fgc_Ngc_unmapp[i]['Ngc']
	lam_unmapp.append([i, float(Fgc_Ngc_unmapp[i]['Fgc'])/Fgc_Ngc_unmapp[i]['Ngc']])
print 'Half Unmappable lambda'
print 'rate 	Fgc 	Ngc'
for i in Fgc_Ngc_unmapp.keys():
	print i, Fgc_Ngc_half_unmapp[i]['Fgc'], Fgc_Ngc_half_unmapp[i]['Ngc']
	lam_unmapp_half.append([i, float(Fgc_Ngc_half_unmapp[i]['Fgc'])/Fgc_Ngc_half_unmapp[i]['Ngc']])
print 'lambda Mappable'
for i in lam_mapp:
	print i

print 'lambda Unmappable'
for i in lam_unmapp:
	print i
print 'lambda Half mapp'
for i in lam_unmapp_half:
	print i
unmapp_coef = polyfit(lam_mapp,lam_unmapp,1)[0]
print 'unmapp_coef =',unmapp_coef









