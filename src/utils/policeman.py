##########################################
## Policeman model related code (not used by now)
##########################################

class PoliceMan:
	def __init__(self, link,chrom, median,initial_direction):
		def init_trans(begin, end, initial_direction, current_direction):
			if initial_direction == 'rf':
				ind = 1
			else:
				ind = 0
			if current_direction != initial_direction:
				self.begin=int(end  - median)
				self.end=int(end) + ind*read_length
			else:
				self.end=int(begin)
				self.begin=int(begin)-median
			return((self.begin, self.end))

		self.name=link.name
		self.num_elements =link.num_elements
		self.chr1 = link.chr1
		self.chr2 = link.chr2
		self.A1 = []
		self.A2 = []
		self.B1 = []
		self.B2 = []
		#self.flank_reg_chr2 = []
		if link.chr1==chrom:
			(self.chr1_begin,self.chr1_end) = init_trans(link.begin, link.rightmost_begin, initial_direction[0], link.direction_type[0])
			(self.chr2_begin,self.chr2_end) = init_trans(link.leftmost_end, link.end,initial_direction[1], link.direction_type[1])
			print 'self.chr1_begin, self.chr1_end = ', self.chr1_begin, self.chr1_end
			print 'self.chr2_begin, self.chr2_end = ', self.chr2_begin, self.chr2_end
		else:
			(self.chr2_begin,self.chr2_end) = init_trans(link.begin, link.rightmost_begin, initial_direction[0], link.direction_type[0])
			(self.chr1_begin,self.chr1_end) = init_trans(link.leftmost_end, link.end,initial_direction[1], link.direction_type[1])
			print 'self.chr1_begin, self.chr1_end = ', self.chr1_begin, self.chr1_end
			print 'self.chr2_begin, self.chr2_end = ', self.chr2_begin, self.chr2_end

	def __repr__(self):
		return repr((self.name, self.num_elements, self.chr1, self.chr2,self.chr1_begin,self.chr1_end,self.chr2_begin,self.chr2_end, self.A1, self.A2, self.B1, self.B2))

def sortPoliceMenBegin(policeman,chrom):
	# chrom = 'chr20'
	# print '!!!!!!!!', policeman
	if chrom == policeman.chr1:
		begin = policeman.chr1_begin
	else:
		begin = policeman.chr2_begin
	return(begin)

def sortPoliceMenEnd(policeman):
	chrom = 'chr20'
	if chrom == policeman.chr1:
		end = policeman.chr1_end
	else:
		end = policeman.chr2_end
	return(end)

#This part is for translocation, we have to take into account only one part of the cluster
def PolicemanCheck(policemen,flank_reg,curr_num,chrom):
	#flag_beg_end='beg' if for flank region A1 or B1
	#flag_beg_end='end'  if for flank region B2 A2
	#if flag_beg_end=='beg':
	#policemen.sort(key=sortPoliceMenEnd)
	for policeman in policemen:
		if policeman.num_elements>=curr_num:
			print 'Ok number is fine'
			print policeman
			if policeman.chr1 == chrom:
				print 'And chromosome chr1 is fine!'
				print flank_reg
				if flank_reg[0]<=policeman.chr1_begin and policeman.chr1_end<=flank_reg[1]:
					print 'I am here I should write to the policeman!!!!!!!!!!!!!!!'
					flank_reg=[policeman.chr1_end,flank_reg[1]]
					policeman.A1 = [flank_reg[0],policeman.chr1_begin]
					policeman.A2 = [policeman.chr1_end, flank_reg[1]]
				elif policeman.chr1_end>flank_reg[1]:
					break
			else:
				print 'And chromosome chr1 is fine!'
				print flank_reg
				if flank_reg[0]<=policeman.chr2_begin and policeman.chr2_end<=flank_reg[1]:
					flank_reg=[policeman.chr2_end,flank_reg[1]]
					policeman.B1 = [flank_reg[0],policeman.chr2_begin]
					policeman.B2 = [policeman.chr2_end, flank_reg[1]]
				elif policeman.chr2_end>flank_reg[1]:
					break	
	#else:
	#	policemen = sorted(policemen,key=lambda PoliceMan: PoliceMan.chr2_begin )
	#	for policeman in policemen:
	#		if policeman.num_elements>=curr_num:
	#			if policeman.chr1 == chrom:
	#				if flank_reg[0]<=policeman.chr1_begin and policeman.end<=flank_reg[1]:
	#					flank_reg=[flank_reg[0],int(policeman.chr1_begin)]
	#					policeman.flank_reg_chr1 = [flank_reg[0],int(policeman.chr1_begin)]
	#				elif policeman.chr1_begin>flank_reg[1]:
	#					break
	#			else:
	#				if flank_reg[0]<=policeman.chr2_begin and policeman.end<=flank_reg[1]:
	#					flank_reg=[flank_reg[0],int(policeman.chr2_begin)]
	#					policeman.flank_reg_chr2 = [flank_reg[0],int(policeman.chr2_begin)]
	#				elif policeman.chr2_begin>flank_reg[1]:
	#					break
	return(flank_reg)

def PolicemanCheck_A2_B1(policemen,flank_reg):
	flank_B1='Empty'
	flank_A2=flank_reg
	policemen.sort(key=sortByBegins)
	for policeman in policemen:
		if policeman.begin>=flank_reg[0] and policeman.end<=flank_reg[1]:
			flank_A2=[flank_reg[0],int(policeman.begin)]
			flank_B1=[int(policeman.end), flank_reg[1]]
		elif policeman.end>flank_reg[1]:
			break
	return(flank_A2, flank_B1)	
