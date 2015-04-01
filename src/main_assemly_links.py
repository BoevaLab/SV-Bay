#from cluster import *
import fragment
import yaml
import optparse
import sys
import joblib
import logging
import bisect
import glob
import itertools
from cluster import * 

file_names_normal =''
parser = optparse.OptionParser()
parser.add_option('-c', '--config', dest='config_file_name', help='Name of configuration file to use')
parser.add_option('-n', '--normal-dir', dest='file_names_normal', help='Name of the directory with clusters generated out of germ-line data')
#
(options, args) = parser.parse_args()
#
config_file = open(options.config_file_name)
config = yaml.load(config_file)
config_file.close()
serialized_stats_file = open(config['working_dir'] + config['serialized_stats_file'])
stats = yaml.load(serialized_stats_file)
serialized_stats_file.close()
biggest_normal = stats['per_chr_stats'][config['chromosomes'][0]]['biggest_normal']
median = stats['per_chr_stats'][config['chromosomes'][0]]['median']
initial_direction = stats['per_chr_stats'][config['chromosomes'][0]]['flag_direction']
read_len = config['read_length']
#######################
initial_direction = 'rf'
def sortLinks(link):
	return int(link.begin)

f = initial_direction[0]
r = initial_direction[1]

#f=open('clusters_result.txt','r')
def IsOverlaps(link,links,tr_flag,biggest_normal):

	overlaps= []
	if not tr_flag:
		for cl in links:
			if cl == link:
				continue
			if link.begin <=cl.begin<=link.end and cl.end >= link.end  and cl.gamma_alelles-3<=link.gamma_alelles<=cl.gamma_alelles+3:
				overlaps.append(cl)
			elif link.end<cl.begin:
				break
	else:
		for cl in links:
			if cl == link:
				continue
			else:
				if cl.chr1 == link.chr1:
					beg =  abs(cl.begin-link.begin)
					end = abs(cl.end-link.end)
					if (beg<=2*biggest_normal or end<=2*biggest_normal): #or (0<link.begin-cl.begin<=2*biggest_normal or link.end-cl.end<= 2*biggest_normal):
						overlaps.append(cl)
				else:
					beg =  abs(cl.end-link.begin)
					end = abs(link.end-cl.begin)
					if (beg<=2*biggest_normal or end<=2*biggest_normal):
					#if beg and end and (cl.end-link.begin<=2*biggest_normal or -link.end+cl.begin<= 2*biggest_normal): #or (0<-cl.end+link.begin<=2*biggest_normal or -link.end+cl.begin<= 2*biggest_normal):
						overlaps.append(cl)
	return overlaps
def IsContain(link,links,biggest_normal,tr_flag):
	contain= []
	if not tr_flag:
		for cl in links:
			if cl == link:
				continue
			elif (0<=cl.begin-link.begin<=2*biggest_normal and link.end - cl.end>=0) or (0<=(link.end - cl.end)<= 2*biggest_normal and cl.begin-link.begin>=0) and cl.gamma_alelles-3<=link.gamma_alelles<=cl.gamma_alelles+3:
				contain.append(cl)
			elif link.end<cl.begin:
				break
	else:
		for cl in links:
			if cl == link:
				continue
			else:
				if cl.chr1 == link.chr1:
					beg =  (link.begin-cl.begin)
					end = (link.end-cl.end)
					if (0<beg<=2*biggest_normal and 0<end<=2*biggest_normal): #or (0<link.begin-cl.begin<=2*biggest_normal or link.end-cl.end<= 2*biggest_normal):
						contain.append(cl)
				else:
					beg =  (link.begin- cl.end)
					end = (link.end-cl.begin)
					if (0<beg<=2*biggest_normal and 0<end<=2*biggest_normal):
					#if beg and end and (cl.end-link.begin<=2*biggest_normal or -link.end+cl.begin<= 2*biggest_normal): #or (0<-cl.end+link.begin<=2*biggest_normal or -link.end+cl.begin<= 2*biggest_normal):
						contain.append(cl)
				#if (0<=cl.begin-link.begin and link.end - cl.end>=0) and (cl.begin-link.begin<= 2*biggest_normal or link.end - cl.end<= 2*biggest_normal):
				#		contain.append(cl)	
				#else:
				#	if (0<cl.end-link.begin<=2*biggest_normal and 0<link.end-cl.begin<= 2*biggest_normal) or (0<-cl.end+link.begin<=2*biggest_normal and 0<-link.end+cl.begin<= 2*biggest_normal):
				#		contain.append(cl)
	return contain

def ReInserionCheck(chr1_pos_1,chr1_pos_2,chr2_pos_1,chr2_pos_2,links_chr1, links_chr2,f,r):
	return_link = []
	for link in links_chr1:
		if (link.begin<=chr1_pos_1<=link.end or link.begin<=chr1_pos_2<=link.end) and link.direction_type == f+r:
			return_link.append(link)
			break
		elif link.begin>=chr1_pos_2 :
			break
	for link in links_chr2:
		if (link.begin<=chr2_pos_1<=link.end or link.begin<=chr2_pos_2<=link.end) and link.direction_type == f+r:
			return_link.append(link)
			break
		elif link.begin>=chr2_pos_2:
			break
	return return_link

def rf_transl(link, links, links_chr1, links_chr2,f,r):
	coins = 0
	sv_links=[]
	overlaps = IsOverlaps(link, curr_links,1,biggest_normal)
	overlp_coin =0
	sv_links.append(link)
	type_Sv_ov=''
	type_Sv_con=''
	if overlaps:
		#if len(overlaps)>1:

		if len(overlaps):
			for overlap in overlaps: 
				if (overlap.direction_type[0] !=  link.direction_type[0] and overlap.direction_type[1] !=  link.direction_type[1]):
		 			reins = ReInserionCheck(link.begin,overlap.begin,link.end,overlap.end, links_chr1, links_chr2,f,r)
		 			if reins:
		 				if len(reins)>1:
		 					print 'Warning ambiguous re-insertion for cluster ',overlap.Cluster.to_string()
		 					type_Sv_ov = 'Unknown'
		 				else:
		 					if not overlp_coin:
		 						type_Sv_ov ='Linking Re-insertion'
		 						overlp_coin+=1
		 						coins+=1
		 					else:
		 						print 'To many overlaps it could be Linking Re-insertion'
		 						type_Sv_ov = 'To many overlaps it could be Linking Re-insertion'
		 						#print type_Sv_ov
		 						coins = 0
		 					sv_links.append(overlaps[0])
		 					sv_links.append(reins[0])
					else:
						if not overlp_coin:
		 					type_Sv_ov = 'Linkin insertion'
		 					coins+=1
		 					overlp_coin+=1
		 					sv_links.append(overlaps[0])
		 				else:
		 					print 'To many overlaps it could be Linking insertion'
		 					type_Sv_ov = 'To many overlaps it could be Linking insertion'
		 					#print type_Sv_ov
		 					coins = 0
	contain = IsContain(link, curr_links,biggest_normal,1) 	
	coins_cont=0		
	if contain:
		if contain:
			for cont_link in contain:
				#for k in contain:
				if (cont_link.direction_type[0] !=  link.direction_type[0] and cont_link.direction_type[1] !=  link.direction_type[1]):
					if not coins_cont:
						type_Sv_con = 'Balanced translocation'
						coins+=1
						#print Cluster.to_string(link)
						coins_cont+=1
						break
					else:
						type_Sv_con = 'To many contains links it could be Balanced translocation'
						coins = 0
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
	#print 'rr_transl'
	sv_links=[]
	coins =0
	type_Sv_con =''
	type_Sv_ov =''
	sv_links.append(link)
	overlp = IsOverlaps(link, links,1,biggest_normal)
	coins_ov=0
	if overlp:
		if len(overlp):
			for i in overlp:
		 		if i.direction_type == f+f :
		 			if not coins_ov:
		 				print 'Balanced translocation with inversion'
		 				type_Sv_ov = 'Balanced translocation with inversion'
		 				sv_links.append(i)
		 				coins+=1
		 				coins_ov+=1
		 			else:

		 				type_Sv_ov = 'To many overlaps could be Balanced translocation with inversion'
		 				sv_links = sv_links[:-1]
		 				coins=0
	contain = IsContain(link, links,biggest_normal,1) 
	coins_cont=0			
	if contain:

		for cont_link in contain:
			if (cont_link.direction_type[0] !=  link.direction_type[0] and cont_link.direction_type[1] !=  link.direction_type[1]):
			 	reins =  ReInserionCheck(link.begin,cont_link.begin,link.end,cont_link.end, links_chr1, links_chr2,f,r)
				if reins:
					if len(reins)>1:
			 			print 'Warning ambiguous re-insertion'
			 			#print contain.cluster.Cluster.to_string()	
			 		else:
			 			if not coins_cont:
			 				type_Sv_con =  'Linking Re-insertion with inversion'
			 				sv_links.append(cont_link)
			 				sv_links.append(reins[0])
			 				#print 'Linking Re-insertion with inversion'
							coins+=1
							coins_cont+=1
						else:
							type_Sv_con = 'To many contains links it could be Linking Re-insertion with inversion'
							coins = 0
							sv_links = sv_links[:-2]
				else:
					type_Sv_con = 'Linking insertion with inversion'
					sv_links.append(cont_link)
					coins+=1
	if coins==1:
		if type_Sv_ov:
			type_sv = type_Sv_ov
		else:
			type_sv = type_Sv_con
	elif not coins:
		type_sv = 'Unbalanced translocations with inversion'
	else:
		type_sv = 'Unknown'
		print 'Could be both '+type_Sv_con+' '+type_Sv_ov
	return  (type_sv,sv_links)

def fr_intra(link, links,f,r,biggest_normal):
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	overlp = IsOverlaps(link,links,0,biggest_normal)
	sv_links =[]
	sv_links.append(link)
	coins = 0
	if overlp:
		#if len(overlp)>2:
		 	#print 'We will see len(overlp) = ',len(overlp)
		 	#for i in overlp:
		 		#print Cluster.to_string(i)
		if len(overlp) == 2:
			if overlp[0].direction_type == r+f and overlp[1].direction_type == r+f and IsOverlaps(overlp[0],overlp[1],0,biggest_normal):
				type_Sv_ov = 'Re-insertion'
				sv_links.append(overlp[0])
				sv_links.append(overlp[1])
				coins+=1
			else:
				print 'Warning ambiguous re-insertion'
		 		#for i in overlp:
		 		#	print Cluster.to_string(i)
		 		type_Sv_ov = 'Unknown'
		else:
		 	if overlp[0].direction_type == r+f:
		 		type_Sv_ov = 'Linkin insertion'
		 		sv_links.append(overlp[0])
		 		coins+=1
	if not type_Sv_ov or type_Sv_ov=='Unknown':
		#print Cluster.to_string(link)
		if link.leftmost_end - link.rightmost_begin>0:
			if int(link.mean_length_frag) > int(biggest_normal):
		 		type_sv = 'Deletion'
		 	#else link.size_type == 'smaller':
		 	else:
		 		type_sv = 'Insertion'
		 	#else:
 			#	type_sv = 'Unknown'
		else:
		 	type_sv = 'Small duplication'
	else:
		#print type_Sv_ov
		type_sv = type_Sv_ov
	return (type_sv,sv_links)

def rf_intra(link, links,f,r,biggest_normal):
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	sv_links =[]
	sv_links.append(link)
	overlp = IsOverlaps(link,links,0,biggest_normal)
	#print '!!!!rf_intra!!!!'
	##print Cluster.to_string(link)
	#print 'link.size_type = ',link.size_type
	#print 'len(IsOverlaps) = ',len(overlp)
	if len(overlp)>=1:
		 	#print 'We will see len(overlp) = ',len(overlp)
		 	#for i in overlp:
		 		#print Cluster.to_string(i)
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
	contain = IsContain(link,links,biggest_normal,0)
	sv_links =[]
	sv_links.append(link)
	#if len(contain)>1:
	#	continue
		#print 'The current link contain to many other links ',len(overlp)
		#for i in overlp:
		#	print Cluster.to_string(i)
		#	type_sv	 = 'Unknown'
	if len(contain)==1:
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
				if not type_sv	 and link.A1!=link.B2 and link.A1!='-1' and link.B2!='-1':
					type_sv	 = 'Mirror duplication'
		else:
			for neibhor in links:
				if link.begin - neibhor.end<= 2*maxlen and neibhor.direction_type == f+f:
					type_sv	 = 'Balanced duplication with inversion on one chromosome'
					sv_links.append(neibhor)
					break
	if not type_sv and link.A1!=link.B2 and link.A1!='-1' and link.B2!='-1':
		type_sv	 = 'Mirror duplication'
	else:
		type_sv = 'Unknown'
	return (type_sv,sv_links)

def ff_intra(link,links,f,r,maxlen):
	coins = 0
	sv_links =[]
	sv_links.append(link)
	type_Sv_con=''
	type_Sv_ov =''
	type_Sv_simple =''
	type_Sv_BD = ''
	#print Cluster.to_string(link)
	contain = IsContain(link,links,biggest_normal,0)

	#if contain:
	#	print 'this link contain something'
	#else:
	#	print 'this link not contain anything'
	#if len(contain)>1:
	#	for i in contain:
	#		print Cluster.to_string(i)
		#type_sv	 = 'Unknown'
	if len(contain)==1:
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
	overlp = IsOverlaps(link,links,0,biggest_normal)
	#print '**********',overlp
	#if len(overlp)>=2:
		#print 'We will see len(overlp) = ',len(overlp)
		#for i in overlp:
		# 	print Cluster.to_string(i)
	if len(overlp)==1:
		#print 'IsOverlaps'
		#print Cluster.to_string(overlp[0])
		if overlp[0].direction_type ==r+r and overlp[0].size_type !='smaller':
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
	if not coins and link.A1!=link.B2 and link.A1!='-1' and link.B2!='-1':
		type_sv = 'Mirror duplication'
	elif coins==1:
		if type_Sv_BD:
			type_sv	 =type_Sv_BD
		elif type_Sv_con:
			type_sv = type_Sv_con
		else:
			type_sv =type_Sv_ov
	else:
		if coins>=2:
			print 'The type is Unknown it could be both : ',type_Sv_con,' and ',type_Sv_ov, ' and ',type_Sv_BD
		type_sv	 = 'Unknown'
	return (type_sv,sv_links)

chromosomes=[]
links=[]


file_names_tumor = glob.glob(config['working_dir'] + config['valid_links_dir'] + '*.txt')
print file_names_tumor	
#Create links out of strings
for fname in file_names_tumor:
	t = open(fname)
	for line in t:
		if line:
			links.append(Cluster.from_string(line))
	t.close()
#file_names_normal =['/Users/dasha/PhD/results_3_tools/all_ga_germ.txt']

if file_names_normal:
	file_names_normal_dir = open(config['working_dir']+file_names_normal)
else:
	file_names_normal =[]
links_normal=[]
for fname in file_names_normal:
	t = open(fname)
	for line in t:
		#print line
		link = Cluster.from_string(line)
		if link.num_elements>=3:
			links_normal.append(link)
	t.close()


#Create set of chromosomes which were used
for link in links:
	if link.chr1 not in chromosomes:
		chromosomes.append(link.chr1)
	if link.chr2 not in chromosomes:
		chromosomes.append(link.chr2)

links.sort(key = sortLinks)
links_normal.sort(key = sortLinks)
#print 'Number of valid links before comparison with normal clusters = ', len(links)
#print 'Number of normal sv = ',len(links_normal)
#for chrom in chromosomes:
#	curr_links_tumor = [link for link in links if link.chr1 == chrom and link.chr2 == chrom]
#	curr_links_normal = [link for link in links_normal if link.chr1== chrom and link.chr2 == chrom]
#	curr_links_normal.sort(key = sortLinks)
#	curr_links_tumor.sort(key = sortLinks)
#	for link_tumor in curr_links_tumor:
#		for link_normal in curr_links_normal:
#			if link_normal.begin - biggest_normal/2> link_tumor.begin:
#				break
#			else:
#				if abs(link_normal.begin-link_tumor.begin)<=biggest_normal/2 and abs(link_normal.length - link_tumor.length)<=biggest_normal/2 and link_normal.direction_type == link_tumor.direction_type and link_normal.name!=link_tumor.name:
#					#print 'Coincidence tumor and normal'
#					#print 'Normal ==>'
#					#print Cluster.to_string(link_normal)
#					#print 'Tumor ==>'
#					#print Cluster.to_string(link_tumor)
#					Cluster.to_string(link_tumor)
#					#links.remove(link_tumor)
#					Cluster.to_string(link_normal)
#					#links_normal.remove(link_normal)
#					curr_links_normal.remove(link_normal)
#
#					break

#print 'Number of valid links after comparison with normal clusters = ', len(links)



print chromosomes
print len(links)
coamplf_lagre_dupl=[]
#Separate coamplifiacation from other SVs
coamplf_lagre_dupl = [link for link in links if link.gamma_alelles > config['numb_allel']]

links = list(set(links).difference(coamplf_lagre_dupl))
#Create set of tranclocations and separate it from other SVs
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

#Process inter-chromosomal SVs

curr_links=[]
#link_file = open('/Users/dasha/PhD/results_3_tools/assembly_ga_nodup1','w')
print 'chrom1','\t', 'bp_pos1','\t','chrom2','\t','bp_pos2','\t','num elements','\t', 'probability','\t', 'type_sv', '\t','direction'

for (chr1,chr2) in itertools.combinations(chromosomes,2):
	links.sort(key = sortLinks)
	for link in translocations:
		if link.chr1 == chr1 and link.chr2 == chr2 or link.chr1 == chr2 and link.chr2 == chr1:
			curr_links.append(link)
	if curr_links:
		#print '***********************************************************************************'
		#print 'Chromosome = ',chr1 ,' and ',chr2, 'Number of translocations with these chromosomes = ',len(curr_links)
		#print '***********************************************************************************'
		curr_links.sort(key = sortLinks)
		#for link in curr_links:
		#	print Cluster.to_string(link)
		#link_file.write('***********************************************************************************'+'\n')
		#link_file.write('Chromosome = '+chr1 +' and '+chr2+'\n')
		#link_file.write('***********************************************************************************'+'\n')
	while curr_links:
		print '-----------------------------------------------------------------------------------'
		link = curr_links[0]
		coins = 0
		contain = []
		overlp =[] 
		if link.direction_type ==f+r or link.direction_type ==r+f:
			(type_sv,sv_links) = rf_transl(link, curr_links, chr_links[chr1], chr_links[chr2],f,r)
		else: 
			(type_sv,sv_links) =rr_transl(link, curr_links, chr_links[chr1], chr_links[chr2],f,r)
		if type_sv!='Unknown':
			#link_file.write('======='+type_sv+'======='+'\n')
			for i in sv_links:
				direction = ''
				if i.direction_type[0]==initial_direction[0]:
					bp_pos1 = ((int(i.begin)+median)+int(i.rightmost_begin))/2
					bp_pos1  = int(i.rightmost_begin)+read_len
					direction+='-'
				else:
					bp_pos1 = ((int(i.rightmost_begin)-median)+int(i.begin))/2
					direction+='+'
					#bp_pos1 = begin
				if i.direction_type[1]==initial_direction[0]:
					bp_pos2 = ((int(i.end)+median)+int(i.leftmost_end))/2
					direction+='-'
					#bp_pos2 = i.leftmost_end
				else:
					bp_pos2 = ((int(i.leftmost_end)-median)+int(i.end))/2
					direction+='+'
					#bp_pos2 = i.end
				print i.chr1,'\t' ,bp_pos1,'\t' ,i.chr2,'\t' , bp_pos2,'\t', i.num_elements,'\t', i.probability,'\t', type_sv, '\t',direction[0], '\t',direction[1]
				#link_file.write(Cluster.to_string(i)+'\n')
		for sv_link in sv_links:
			if sv_link.chr1==sv_link.chr2:
				links = list(set(links).difference(sv_links))
			else:
				curr_links = list(set(curr_links).difference(sv_links))
		#print '!!!', type_sv
links = [i for i in links if i.num_elements!=3]
#Process intra-chromosomal SVs
for chrom in chromosomes:
	#print '***********************************************************************************'
	#print 'Chromosome = ',chrom
	#print '***********************************************************************************'
	links = chr_links[chrom]
	#link_file.write('***********************************************************************************'+'\n')
	#link_file.write('Chromosome = '+chrom+'\n')
	#link_file.write('***********************************************************************************'+'\n')
	while links:
		links.sort(key = sortLinks)
		link = links[0]
		
		if link.direction_type == r+r:
			(type_sv, sv_links) = rr_intra(link,links,f,r,biggest_normal)
			#print 'r+r_intra'
		elif link.direction_type == r+f:
			(type_sv, sv_links) = rf_intra(link,links,f,r,biggest_normal)
			#print 'r+f_intra'
		elif link.direction_type == f+f:
			#print 'f+f_intra'
			(type_sv, sv_links) = ff_intra(link,links,f,r,biggest_normal)
		elif link.direction_type == f+r:
			#print 'f+r_intra'
			(type_sv, sv_links) = fr_intra(link,links,f,r,biggest_normal)
		
		links = list(set(links).difference(sv_links))

		if type_sv!='Unknown':
			#link_file.write('======='+type_sv+'======='+'\n')
			#if type_sv =='Small insertion' and abs(int(i.mean_length_frag)- biggest_normal)>1200:
			if 'Small' in type_sv or (type_sv == 'Deletion' and abs(int(link.mean_length_frag)- median)<1200):
				continue
			else:
				print '-----------------------------------------------------------------------------------'
				#print '!!!', type_sv

				for i in sv_links:

					direction=''
					if i.direction_type[0]==initial_direction[0]:
						if type_sv == 'Insertion':
							bp_pos1  = int(i.rightmost_begin)+read_len
						else:
							bp_pos1 = ((int(i.begin)+median)+int(i.rightmost_begin))/2
						#bp_pos1  = int(i.rightmost_begin)+read_len
						direction+='-'
					else:
						bp_pos1 = ((int(i.rightmost_begin)-median)+int(i.begin))/2
						direction+='+'
						#bp_pos1 = begin
					if i.direction_type[1]!=initial_direction[1]:
						bp_pos2 = ((int(i.end)+median)+int(i.leftmost_end))/2
						direction+='-'
						#bp_pos2 = i.leftmost_end
					else:
						if type_sv == 'Insertion':
							bp_pos2 = (int(i.leftmost_end))
						else:
							bp_pos2 = ((int(i.leftmost_end)-median)+int(i.end))/2
						direction+='+'
						#bp_pos2 = i.end
					print i.chr1,'\t' ,bp_pos1,'\t' ,i.chr2,'\t' , bp_pos2,'\t',int(bp_pos2)-int(bp_pos1),'\t', i.num_elements,'\t', i.probability,'\t',type_sv, '\t',direction[0], '\t',direction[1]
					#	print Cluster.to_string(i)'\t'
				#for i in sv_links:
					#link_file.write(Cluster.to_string(i)+'\n')

#Process coamplifications

coamplifications =[]

coamplf_lagre_dupl = [link for link in coamplf_lagre_dupl if link.num_elements>30]
coamplf_lagre_dupl.sort(key = sortLinks)
#coamp = [coamplf_lagre_dupl[0]]
coamp_link = coamplf_lagre_dupl[0]
#for k in coamplf_lagre_dupl:
#		print Cluster.to_string(k)
coamp =[]
while coamplf_lagre_dupl:
	#coamp_link = coamplf_lagre_dupl[0]
	coamplf_lagre_dupl.remove(coamp_link)
	coamp.append(coamp_link)
	#coamplf_lagre_dupl.remove(coamp_link)
	flag = 0
	curr_links = [link for link in coamplf_lagre_dupl if link.chr1==coamp_link.chr1 or link.chr2 == coamp_link.chr1 or link.chr1==coamp_link.chr2 or link.chr2 == coamp_link.chr2]
	if curr_links:
			for neibhor_link in curr_links:
				#if coamp_link.direction_type[0]==initial_direction[0]:
				if neibhor_link.chr1 == coamp_link.chr1:
					if neibhor_link.direction_type[1]!=coamp_link.direction_type[0]: #and neibhor_link.begin>coamp_link.begin:
						coamp_link = neibhor_link
						#print 'First---',Cluster.to_string(neibhor_link)
						flag = 1
						break
				if neibhor_link.chr2 == coamp_link.chr1:
					if neibhor_link.direction_type[0]!=coamp_link.direction_type[0]:#and neibhor_link.end>coamp_link.begin:
						coamp_link = neibhor_link
						#print 'Second---',Cluster.to_string(neibhor_link)
						flag = 1
						break
				if neibhor_link.chr2 == coamp_link.chr2:
					if neibhor_link.direction_type[1]!=coamp_link.direction_type[1]:#and neibhor_link.end>coamp_link.begin:
						coamp_link = neibhor_link
						#print 'Second2---',Cluster.to_string(neibhor_link)
						flag = 1
						break
				if neibhor_link.chr1 == coamp_link.chr2:
					if neibhor_link.direction_type[0]!=coamp_link.direction_type[1]:#and neibhor_link.end>coamp_link.begin:
						coamp_link = neibhor_link
						#print 'First2---',Cluster.to_string(neibhor_link)
						flag = 1
						break
	if not flag:
		print '---------------------------------------------------------------------------------'
		#print 'Coamplification'
		for i in coamp:
			direction=''
			if i.direction_type[0]==initial_direction[0]:
				bp_pos1 = ((int(i.begin)+median)+int(i.rightmost_begin))/2
				bp_pos1  = int(i.rightmost_begin)+read_len
				direction+='-'
			else:
				bp_pos1 = ((int(i.rightmost_begin)-median)+int(i.begin))/2
				direction+='+'
				#bp_pos1 = begin
			if i.direction_type[1]==initial_direction[0]:
				bp_pos2 = ((int(i.end)+median)+int(i.leftmost_end))/2
				direction+='-'
				#bp_pos2 = i.leftmost_end
			else:
				bp_pos2 = ((int(i.leftmost_end)-median)+int(i.end))/2
				direction+='+'
				#bp_pos2 = i.end
			print i.chr1,'\t' ,bp_pos1,'\t' ,i.chr2,'\t' , bp_pos2,'\t', i.num_elements,'\t', i.probability,'\t','coamplification', '\t',direction[0], '\t',direction[1]
		print '-----------------------------------------------------------------------------------'	
		if coamplf_lagre_dupl:
			coamp_link = coamplf_lagre_dupl[0]
		coamp = []
		#print 'New coamp link = ',
		#if coamplf_lagre_dupl:
		#	for i in coamp:
		#		if i in coamplf_lagre_dupl:
		#			coamplf_lagre_dupl.remove(i)
		#if coamplf_lagre_dupl:	
		#	coamp = [coamplf_lagre_dupl[0]]
			#	print Cluster.to_string(i)'\t'

			

#link_file.close()














