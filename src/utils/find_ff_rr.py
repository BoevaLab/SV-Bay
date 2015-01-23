
#name_frag= ['read_4283_33181274', 'read_4283_22855075', 'read_4283_33219134', 'read_4283_22855077', 'read_4283_22855080', 'read_4283_22855082', 'read_4283_22855085', 'read_400_1270534', 'read_4283_23002241', 'read_4283_22855089', 'read_4283_22855094', 'read_4283_22694657', 'read_4283_22855099', 'read_4283_22809254', 'read_4283_22247178', 'read_4283_22855105', 'read_4283_22855108', 'read_4283_12056067', 'read_4283_11320250', 'read_4283_22775445']
#This should be from_sam_file
name_frag = ['read_4283_12024048', 'read_4283_22248258', 'read_4283_22651335', 'read_4283_22973131', 'read_4283_22264439', 'read_4283_11355601', 'read_4283_12046565', 'read_4283_12016544', 'read_4283_11355604', 'read_4283_22961949', 'read_4283_11442633', 'read_4283_12183694', 'read_4283_11442636', 'read_4283_11355606', 'read_4283_22302186', 'read_4283_23019381', 'read_4283_22893522', 'read_4283_11442660', 'read_4283_11442662', 'read_4283_11442664']
#name_frag = ['read_4283_12024048', 'read_4283_22248258', 'read_4283_22651335', 'read_4283_22973131', 'read_4283_22264439', 'read_4283_11355601', 'read_4283_12046565', 'read_4283_12016544', 'read_4283_11355604', 'read_4283_22961949', 'read_4283_11442633', 'read_4283_12183694', 'read_4283_11442636', 'read_4283_11355606', 'read_4283_22302186', 'read_4283_23019381', 'read_4283_22893522', 'read_4283_11442660', 'read_4283_11442662', 'read_4283_11442664']
def from_pirs_ind(name_frag,f):
	f = open('/Users/daria/donor/donor_RF_70_4283.read.info','r')
	flag=0
	for i in f:
		if name_frag==[]:
			break
		if i[0]=='@':
			if i.split('/')[0][1:] in name_frag:
				print i
				flag+=1
				if flag==2:
					name_frag.remove(i.split('/')[0][1:])
					flag=0
	f.close()

def from_sam_file(name_frag,f_name):
	flag_print_any_way=0
	f = open(f_name,'r')
	flag =0
	for line in f:
		if flag_print_any_way:
			print 'flag_print_any_way'
			print line
		if name_frag==[]:
			break
		if line.split()[0] in name_frag:
			print line
			flag+=1
			flag_print_any_way=1
			if flag==2:
				name_frag.remove(line.split()[0])
				flag=0
				flag_print_any_way=0
	f.close()

f_name = '/Users/daria/Dropbox/data/sam_files/simulate/by_chr/chr1_filt.sam'
from_sam_file(name_frag,f_name)