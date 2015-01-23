import sys

####################################################################
# SubLink class
####################################################################

class SubLink:
	def __init__(self, link, flag_beg_end, biggest_normal,initial_direction):
		self.link = link
		self.left_neighbor_end = 0
		self.right_neighbor_begin = sys.maxint
		self.left_neighbor_name = ''
		self.right_neighbor_name = ''
		self.num_elements = link.num_elements
		self.chrom = ''
		if flag_beg_end == 'beg':
			self.name = link.name + '_1'
			self.chrom = link.chr1
		else:	
			self.name = link.name + '_2'
			self.chrom = link.chr2
		if flag_beg_end == 'beg':
			if link.direction_type[0] == initial_direction[0]:
				self.direction = initial_direction[0]
				self.safe_start = link.rightmost_begin - biggest_normal
				self.safe_end = link.begin + biggest_normal	
			else:
				self.direction = initial_direction[1]
				self.safe_start = link.rightmost_begin - biggest_normal
				self.safe_end = link.begin + biggest_normal
		else:
			
			if link.direction_type[1] == initial_direction[1]:
				self.direction = initial_direction[1]
				self.safe_start = link.end - biggest_normal
				self.safe_end = link.leftmost_end + biggest_normal
			else:
				self.direction = initial_direction[0]
				self.safe_start = link.end - biggest_normal
				self.safe_end = link.leftmost_end + biggest_normal

def sortSubLinksEnd(sub_link):
	return sub_link.safe_end

def sortSubLinksBeg(sub_link):
	return sub_link.safe_start

def sortSubLinksName(sub_link):
	return sub_link.name