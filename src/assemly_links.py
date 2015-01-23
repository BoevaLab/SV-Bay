#from cluster import *
import fragment
import yaml
import optparse
import sys
import joblib
import logging
import bisect
import glob
from cluster import * 


parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')

(options, args) = parser.parse_args()

config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()
biggest_normal = 5182
#######################
initial_direction = 'rf'
def sortLinks(link):
	return int(link.begin)

f = initial_direction[0]
r = initial_direction[1]

#f=open('clusters_result.txt','r')
def IsOverlaps (link,links):
	overlaps= []
	for cl in links:
		if cl == link:
			continue
		if link.begin <=cl.begin<=link.end and cl.gamma_alelles-2<=link.gamma_alelles<=cl.gamma_alelles+2:
			overlaps.append(cl)
		elif link.end<cl.begin:
			break
	return overlaps
def IsContain(link,links,biggest_normal):
	contain= []
	for cl in links:
		if cl == link:
			continue
		elif (0<=cl.begin-link.begin<=2*biggest_normal and link.end - cl.end>=0) or (0<=(link.end - cl.end)<= 2*biggest_normal and cl.begin-link.begin>=0) and cl.gamma_alelles-2<=link.gamma_alelles<=cl.gamma_alelles+2:
			contain.append(cl)
		elif link.end<cl.begin:
			break
	return contain
def ReInserionCheck(chr1_pos_1,chr1_pos_2,chr2_pos_1,chr2_pos_2,links_chr1, links_chr2,f,r):
	return_link = []
	for link in links_chr1:
		if link.begin<=chr1_pos_1<=link.end and link.begin<=chr1_pos_2<=link.end and link.direction_type == f+r:
			return_link.append(link)
			break
		elif link.begin>=chr1_pos_1 or link.begin>=chr1_pos_2:
			break
	for link in links_chr2:
		if link.begin<=chr2_pos_1<=link.end and link.begin<=chr2_pos_2<=link.end and link.direction_type == f+r:
			return_link.append(link)
			break
		elif link.begin>=chr2_pos_1 or link.begin>=chr2_pos_2:
			break
	return return_link


def rf_transl(link, links, links_chr1, links_chr2,f,r):
	print 'rf_transl'
	coins = 0
	sv_links=[]
	overlaps = IsOverlaps(link, curr_links)
	overlp_coin =0
	sv_links.append(link)
	if overlp:
		if len(overlp)>1:
		 	print 'We will see len(overlp) = ',len(overlp)
		 	for i in overlp:
		 		print Cluster.to_string(i)
		if len(overlp):
			for overlap in  overlaps: 
				if overlap.direction_type == r+f:
		 			reins = ReInserionCheck(link.begin,overlp.begin,link.end,overlp.end, links_chr1, links_chr2,f,r)
		 			if reins:
		 				if len(reins)>1:
		 					print 'Warning ambiguous re-insertion for cluster ',overlp.Cluster.to_string()
		 					#print overlp.Cluster.to_string()
		 					#type_Sv_ov = 'Unknown'
		 				else:
		 					if not overlp_coin:
		 						type_Sv_ov ='Linking Re-insertion'
		 						#print 'Linking Re-insertion'
		 						#print overlp[0].cluster.Cluster.to_string()
		 						#print reins[0].cluster.Cluster.to_string()
		 						overlp_coin+=1
		 						coins+=1
		 					else:
		 						type_Sv_ov = 'To many overlaps it could be Linking Re-insertion'
		 						coins = 0
		 					sv_links.append(overlp[0])
		 					sv_links.append(reins[0])
					else:
		 				type_Sv_ov = 'Linkin insertion'
		 				coins+=1
		 	contain = IsContain(link, curr_links,biggest_normal) 			
		 	if contain:
		 		if len(contain)==1 and contain[0].direction_type == r+f:
		 			type_Sv_con = 'Balanced translocation'
		 			coins+=1
		 			print link.cluster.Cluster.to_string()
					print contain[0].cluster.Cluster.to_string()
				elif len(contain)>1:
		 			#print 'Unknown SV, contain more then one'
		 			type_Sv_con = 'Unknown'
		 			for i in contain:
		 				print link.cluster.Cluster.to_string()
		 				print Cluster.to_string(i)
	if coins==1:
		if type_Sv_ov:
			type_sv = type_Sv_ov
		else:
			type_sv = type_Sv_con
	elif not coins:
		type_sv = 'Unbalanced translocation'
	else:
		type_sv = 'Unknown'
		print 'Could be both '+type_Sv_con+' '+type_Sv_ov
	return  (type_sv,sv_links)
def rr_transl(link, links,links_chr1, links_chr2,f,r):
	print 'rr_transl'
	sv_links=[]
	type_Sv_con =''
	type_Sv_ov =''
	sv_links.append(link)
	overlp = IsOverlaps(link, links)
	if overlp:
		if len(overlp)>1:
		 	print 'We will see len(overlp) = ',len(overlp)
		 	for i in overlp:
		 		print Cluster.to_string(i)
		else:
		 	if overlp[0].direction_type == f+f:
		 		type_Sv_ov = 'Balanced translocation with inversion'
		 		sv_links.append(IsOverlaps[0])
		 		coins+=1
	contain = IsContain(link, links,biggest_normal) 			
	if contain:
		if len(contain)==1 and contain[0].direction_type == f+f:
		 	reins =  ReInserionCheck(link.begin,overlp.begin,link.end,overlp.end, links_chr1, links_chr2,f,r)
			if reins:
				if len(reins)>1:
		 			print 'Warning ambiguous re-insertion'
		 			print overlp.cluster.Cluster.to_string()	
		 		else:
		 			type_Sv_con =  'Linking Re-insertion with inversion'
		 			sv_links.append(overlp[0])
		 			sv_links.append(reins[0])
		 			print 'Linking Re-insertion with inversion'
		 			print overlp[0].cluster.Cluster.to_string()
					print reins[0].cluster.Cluster.to_string()
					coins+=1
		elif len(contain)>1:
			print 'Unknown SV,contain more then one'
		 	for i in contain:
		 		print link.cluster.Cluster.to_string()
		 		print Cluster.to_string(i)
	if not type_Sv_con and not type_Sv_ov:
		type_sv = type_Sv_con
		print type_sv
		print Cluster.to_string(link)
	elif type_Sv_con:
		print type_Sv_con
		type_sv = type_Sv_con
	elif type_Sv_ov:
		print type_Sv_ov
		type_sv = type_Sv_ov
		print link.cluster.Cluster.to_string()
	else:
		print 'Type of sv is Unknown'
		type_sv = 'Unknown'
		print 'There are two possibilities : '
		print type_Sv_ov +' and '+type_Sv_con
	return (type_sv,sv_links)
def ff_transl(link, links,links_chr1, links_chr2,f,r):
	print 'ff_transl'
	type_Sv_con = ''
	type_Sv_ov = ''
	coins =0
	sv_links.append(link)
	overlp = IsOverlaps(link, links)
	if overlp:
		if len(overlp)>1:
		 	print 'We will see len(overlp) = ',len(overlp)
		 	for i in overlp:
		 		print Cluster.to_string(i)
		else:
		 	if overlp[0].direction_type == f+f:
		 		type_Sv_ov = 'Balanced translocation with inversion'
		 		sv_links.append(overlp[0])
		 		coins+=1
	contain = IsContain(link, links,biggest_normal)
	if contain:
		if len(contain)==1 and contain[0].direction_type == f+f:
		 	reins =  ReInserionCheck(link.begin,contain[0].begin,link.end,contain[0].end, links_chr1, links_chr2,f,r)
			if reins:
				if len(reins)>1:
		 			print 'Warning ambiguous re-insertion'
		 			for i in contain:
		 				print Cluster.to_string(i)
		 			for i in reins:
		 				print Cluster.to_string(i)
		 		else:
		 			print 'Linking Re-insertion with inversion'
		 			print Cluster.to_string(contain[0])
					print Cluster.to_string(reins[0])
					links.append(reins[0])
					links.append(contain[0])
					type_Sv_con = 'Linking Re-insertion with inversion'
					coins+=1

			else:
		 		type_Sv = 'Linkin insertion'
		 		links.append(contain[0])
		 		coins+=1
		elif len(contain)>1:
			print 'Unknown SV,contain more then one'
		 	for i in contain:
		 		print link.cluster.Cluster.to_string()
		 		print Cluster.to_string(i)
	if not type_Sv_con and not type_Sv_ov:
		type_sv = 'Unbalanced translocation with insertion'
		print type_sv
		print Cluster.to_string(link)
	else:
		if coins ==1:
		 	if type_Sv_ov:
		 		type_sv = type_Sv_ov
		 		print type_sv
		 	else:
		 		type_sv = type_Sv_con
		 		print type_sv
		else:
		 	print 'Type of sv is Unknown'
		 	print 'There are two possibilities : '
		 	print type_Sv_ov +' and '+type_Sv_con
		 	type_sv = 'Unknown'
	return (type_sv,sv_links)	

def fr_transl(link, links,links_chr1, links_chr2,f,r):
	print 'fr_transl'
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	overlp = IsOverlaps(link, links)
	sv_links= []
	sv_links.append(link)
	if overlp:
		if len(overlp)>1:
		 	print 'We will see len(overlp) = ',len(overlp)
		 	for i in overlp:
		 		print Cluster.to_string(i)
		else:
			if IsOverlaps[0].direction_type == r+f:
		 		reins =  ReInserionCheck(link.begin,overlp.begin,link.end,overlp.end, links_chr1, links_chr2,f,r)
		 		if reins:
		 			if len(reins)>1:
		 				print 'Warning ambiguous re-insertion'
		 				print overlp.cluster.Cluster.to_string()
		 				type_Sv_ov = 'Unknown'
		 			else:
		 				type_Sv_ov = 'Linking Re-insertion with inversion'
		 				print 'Linking Re-insertion with inversion'
		 				print overlp[0].cluster.Cluster.to_string()
						print reins[0].cluster.Cluster.to_string()
						coins+=1
						sv_links.append(overlp[0])
						sv_links.append(reins[0])
				else:
					sv_links.append(overlp[0])
		 			type_Sv_ov = 'Linkin insertion'
		 			coins+=1
	contain = IsContain(link, links,biggest_normal)
	if contain:
		if len(contain)==1 and contain[0].direction_type == r+f:
		 		type_Sv_con = 'Balanced translocation'
		 		coins+=1
		 		sv_links.append(contain[0])
		elif len(contain)>1:
			print 'Unknown SV,contain more then one'
		 	type_Sv_con = 'Unknown'
		 	for i in contain:
		 		print link.cluster.Cluster.to_string()
		 		print Cluster.to_string(i)
	if coins == 2:
		print 'Type of sv is Unknown'
		print 'There are two possibilities : '
		print type_Sv_ov +' and '+type_Sv_con
	 	type_Sv = 'Unknown'
	elif coins==1:
		if type_Sv_con == 'Unknown':
	 		type_sv = type_Sv_ov
	 	else:
	 		type_sv = type_Sv_con
	else:
		type_sv = 'Unbalanced translocation'
	
	return (type_sv,sv_links)

def fr_intra(link, links,f,r,biggest_normal):
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	overlp = IsOverlaps(link,links)
	sv_links =[]
	sv_links.append(link)
	coins = 0
	if overlp:
		if len(overlp)>2:
		 	print 'We will see len(overlp) = ',len(overlp)
		 	for i in overlp:
		 		print Cluster.to_string(i)
		elif len(overlp) == 2:
			if overlp[0].direction_type == r+f and overlp[1].direction_type == r+f and IsOverlaps(overlp[0],overlp[1]):
				type_Sv_ov = 'Re-insertion'
				sv_links.append(overlp[0])
				sv_links.append(overlp[1])
				coins+=1
			else:
				print 'Warning ambiguous re-insertion'
		 		for i in overlp:
		 			print Cluster.to_string(i)
		 		type_Sv_ov = 'Unknown'
		else:
		 	if overlp[0].direction_type == r+f:
		 		type_Sv_ov = 'Linkin insertion'
		 		sv_links.append(overlp[0])
		 		coins+=1
	if not type_Sv_ov or type_Sv_ov=='Unknown':
		if link.leftmost_end - link.rightmost_begin>0:
			if link.size_type == 'bigger':
		 		type_sv = 'Deletion'
		 	elif link.size_type == 'smaller':
		 		type_sv = 'Small insertion'
		 	else:
 				type_sv = 'Unknown'
		else:
		 	type_sv = 'Small duplication'

		#if type_Sv_ov == 'Unknown':
		#	print type_Sv_simple
		#	type_sv = type_Sv_simple
	else:
		print type_Sv_ov
		type_sv = type_Sv_ov
	return (type_sv,sv_links)

def rf_intra(link, links,f,r,biggest_normal):
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	sv_links =[]
	sv_links.append(link)
	overlp = IsOverlaps(link,links)
	print 'len(IsOverlaps) = ',len(overlp)
	if len(overlp)>=1:
		 	print 'We will see len(overlp) = ',len(overlp)
		 	for i in overlp:
		 		print Cluster.to_string(i)
		 	type_sv = 'Unknown'
	else:
		if link.size_type=='bigger':
			type_sv = 'Large duplication'
		else:
			type_sv = 'Unknown'
	
	return (type_sv,sv_links)

def rr_intra(link,links,f,r,maxlen):
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	type_sv = ''
	contain = IsContain(link,links,biggest_normal)
	sv_links =[]
	sv_links.append(link)
	if len(contain)>1:
		print 'The current link contain to many other links ',len(overlp)
		for i in overlp:
			print Cluster.to_string(i)
			type_sv	 = 'Unknown'
	elif len(contain)==1:
		if contain[0].direction_type == f + f:
			sv_links.append(contain[0])
			if link.size_type == 'bigger':
				if contain[0].size_type =='bigger':
					type_sv	 = 'Linking insertion with inversion'
				else:
					type_sv	 ='Large duplication with inversion'
			else:
				for neibhor in links:
					sv_links.append(contain[0])
					if link.begin - neibhor.end<= 2*maxlen and neibhor.direction_type == f+f:
						type_sv	 = 'Balanced duplication with inversion on one chromosome'
						sv_links.append(neibhor)
						break
				if not type_sv	 :
					type_sv	 = 'Mirror duplication'
		else:
			for neibhor in links:
				if link.begin - neibhor.end<= 2*maxlen and neibhor.direction_type == f+f:
					type_sv	 = 'Balanced duplication with inversion on one chromosome'
					sv_links.append(neibhor)
					break
	if not type_sv:
		type_sv	 = 'Mirror duplication'
	return (type_sv,sv_links)

def ff_intra(link,links,f,r,maxlen):
	coins = 0
	sv_links =[]
	sv_links.append(link)
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	type_Sv_BD = ''
	contain = IsContain(link,links,biggest_normal)
	if contain:
		print 'this link contain something'
	else:
		print 'this link not contain anything'
	if len(contain)>1:
		for i in contain:
			print Cluster.to_string(i)
		#type_sv	 = 'Unknown'
	elif len(contain)==1:
		if contain[0].direction_type == r+r:
			if link.size_type == 'bigger':
				sv_links.append(contain[0])
				if contain[0].size_type =='bigger':
					type_Sv_con = 'Linking insertion with inversion'
					coins+=1
				elif contain[0].size_type =='smaller':
					type_Sv_con ='Large duplication with inversion'
					coins+=1
	
	#if not type_Sv_con:
	#	type_Sv_simple = 'Mirror duplication'
	overlp = IsOverlaps(link,links)
	print '**********',overlp
	if len(overlp)>=2:
		print 'We will see len(overlp) = ',len(overlp)
		for i in overlp:
		 	print Cluster.to_string(i)
	elif len(overlp)==1:
		print 'IsOverlaps'
		if overlp[0].direction_type ==r+r and overlp[0].size_type =='bigger':
			sv_links.append(overlp[0])
			type_Sv_ov = 'Basic inversion'
			coins+=1
			#print type_Sv_ov 
	for neibhor in links:
		if neibhor==link:
			continue
		elif 0<-link.end + neibhor.begin<= 2*maxlen and neibhor.direction_type == r+r:
			sv_links.append(neibhor)
			type_Sv_BD = 'Balanced duplication with inversion on one chromosome'
			coins+=1
			break
	if not coins:
		type_sv = 'Mirror duplication'
	elif coins>=2:
		print 'The type is Unknown it could be both : ',type_Sv_con,' and ',type_Sv_ov, ' and ',type_Sv_BD
		type_sv	 = 'Unknown'
	elif coins==1:
		if type_Sv_BD:
			type_sv	 =type_Sv_BD
		elif type_Sv_con:
			type_sv = type_Sv_con
		else:
			type_sv =type_Sv_ov
	return (type_sv,sv_links)




chromosomes=[]
#chromosomes = ['chr15']
links=[]
#for i in range(1,23,1):
#	chromosomes.append('chr'+str(i))
file_names = glob.glob(config['working_dir'] + config['valid_links_dir'] + '*.txt')
#chromosomes = ['chr15']
print chromosomes
#file_names = [config['working_dir'] +'chr15_chr15_valid_links.txt']
for fname in file_names:
	t = open(fname)
	for line in t:
		links.append(Cluster.from_string(line))
	t.close()
for link in links:
	if link.chr1 not in chromosomes:
		chromosomes.append(link.chr1)
	if link.chr2 not in chromosomes:
		chromosomes.append(link.chr2)
print chromosomes
print len(links)
coamplf_lagre_dupl=[]
#coamplf_lagre_dupl = [link for link in links if link.gamma_alelles > config['numb_allel']]
for link in links: 
	if link.gamma_alelles > int(config['numb_allel']):
		coamplf_lagre_dupl.append(link)
		print Cluster.to_string(link)
links = list(set(links).difference(coamplf_lagre_dupl))
chr_links = dict()
translocations =[]
for chrom in chromosomes:
	chr_links[chrom]=[]
	for link in links:
		if link.chr1 ==link.chr2 == chrom:
			chr_links[chrom].append(link)
		elif link.chr1 !=link.chr2:
			translocations.append(link)
	links = list(set(links).difference(chr_links[chrom]))
	links = list(set(links).difference(translocations))
chrom_combinations = []
#chromosomes = []
#links = sort(key = lambda link: link.begin)
links.sort(key=sortLinks)
curr_links=[]
link_file = open('Assembled_links','w')
for chr1 in chromosomes:
	for chr2 in chromosomes:
		if chr1==chr2:
			continue
		else:
			for link in translocations:
		 		if link.chr1==chr1 and link.chr2 == chr2:
		 			curr_links.append(link)
		 	if curr_links:
		 		print '***********************************************************************************'
				print 'Chromosome = ',chr1 ,' and ',chr2
				print '***********************************************************************************'
				link_file.write('***********************************************************************************'+'\n')
				link_file.write('Chromosome = '+chr1 +' and '+chr2+'\n')
				link_file.write('***********************************************************************************'+'\n')
		 	while curr_links:
		 		link = curr_links[0]
		 		coins = 0
		 		contain = []
		 		overlp =[] 
		 		if link.direction_type ==f+r:
		 			(type_sv,sv_links) = rf_transl(link, curr_links, chr_links[chr1], chr_links[chr2],f,r)
		 		elif link.direction_type ==r+r:
		 			(type_sv,sv_links) =rr_transl(link, curr_links, chr_links[chr1], chr_links[chr2],f,r)
		 		elif link.direction_type == r+f:
		 			(type_sv,sv_links) =rf_transl(link, curr_links, chr_links[chr1], chr_links[chr2],f,r)
		 		elif link.direction_type == f+f:
		 			(type_sv,sv_links) =ff_transl(link, curr_links, chr_links[chr1], chr_links[chr2],f,r)
		 		if type_sv!='Unknown':
					link_file.write('======='+type_sv+'======='+'\n')
					for i in sv_links:
						link_file.write(Cluster.to_string(i)+'\n')
		 		for sv_link in sv_links:
		 			if sv_link.chr1==sv_link.chr2:
		 				links = list(set(sv_links).difference(links))
					else:
		 				curr_links = list(set(sv_links).difference(curr_links))
print chromosomes

for chrom in chromosomes:
	print '***********************************************************************************'
	print 'Chromosome = ',chrom
	print '***********************************************************************************'
	links = chr_links[chrom]
	link_file.write('***********************************************************************************'+'\n')
	link_file.write('Chromosome = '+chrom+'\n')
	link_file.write('***********************************************************************************'+'\n')
	
	while links:
		links.sort(key = sortLinks)
		link = links[0]
		if link.direction_type == r+r:
			(type_sv, sv_links) = rr_intra(link,links,f,r,biggest_normal)
			print '-----------------------------------------------------------------------------------'
			print 'r+r_intra'
			print '!!!', type_sv
			for i in sv_links:
				print Cluster.to_string(i)

		elif link.direction_type == r+f:
			(type_sv, sv_links) = rf_intra(link,links,f,r,biggest_normal)
			print '-----------------------------------------------------------------------------------'
			print 'r+f_intra'
			print '!!!', type_sv

			for i in sv_links:
				print Cluster.to_string(i)
			print '-----------------------------------------------------------------------------------'
		elif link.direction_type == f+f:
			print '-----------------------------------------------------------------------------------'
			print 'f+f_intra'
			(type_sv, sv_links) = ff_intra(link,links,f,r,biggest_normal)
			print '!!!', type_sv
			for i in sv_links:
				print Cluster.to_string(i)
			print '-----------------------------------------------------------------------------------'
		elif link.direction_type == f+r:
			print '-----------------------------------------------------------------------------------'
			print 'f+r_intra'
			(type_sv, sv_links) = fr_intra(link,links,f,r,biggest_normal)
			print '!!!', type_sv
			for i in sv_links:
				print Cluster.to_string(i)
			print '-----------------------------------------------------------------------------------'
		links = list(set(links).difference(sv_links))
		if type_sv!='Unknown':
			link_file.write('======='+type_sv+'======='+'\n')
			for i in sv_links:
				link_file.write(Cluster.to_string(i)+'\n')

print 'Coamplifications : '
for i in coamplf_lagre_dupl:
	Cluster.to_string(i)
	link_file.write(Cluster.to_string(i)+'\n')
link_file.close()















