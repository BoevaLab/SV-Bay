import logging
import pysam

logger = logging.getLogger('main_logger')

# Helpers to check fragment flags

# Check mapping quality
# AM is basically always present,
# but it's safer to check for errors
# and return 0 if it's not
def has_good_mq(read):
	try:
		return 1 if int(read.opt('AM')) >= 20 else 0
	except KeyError:
		return 0

# We use XT to determine that read is uniquely mapped
# instead of mapping quality, because mapq may take into account read lenghts
# which results in giving low mapq to structure variants and though skipping them
# Unfortunately, XT is not always available (for example, 
# it's available in the output of bwa, but not in the output of novoalign)
# If it's not available we just set unique flag for all fragments - 
# it's better to process several percent of badly mapped reads, 
# than to loose some of the structure variants
def is_unique(read):            
	try:
		return 1 if read.opt('XT') == 'U' else 0
	except KeyError:
		return 1

####################################################################
# Fragment class
####################################################################
class Fragment(object):
	__slots__ = ('first_read_chr', 'second_read_chr', 'begin', 'end', 'length', 'middle', 'direction', 'name','unique_flag','mapp_qul_flag')
		
	def from_reads(self, read_a, read_b, sam_in):
		# r1 and r2 are ordered (r1 - left, r2 - right)
		# Different ordering logic for
		# intra-chr fragment and translocation
		if ((read_a.tid == read_a.rnext and read_a.pos < read_a.pnext) or \
			(read_a.tid != read_a.rnext and sam_in.getrname(read_a.tid) < sam_in.getrname(read_b.tid))):
			(r1, r2) = (read_a, read_b)
		else:
			(r1, r2) = (read_b, read_a)

		self.first_read_chr = sam_in.getrname(r1.tid)
		self.second_read_chr = sam_in.getrname(r2.tid)
		self.begin = r1.pos
		self.end = r2.pos + r2.qlen
		self.length = self.end - self.begin
		self.middle = (self.begin + self.end) / 2
		self.direction = 'r' if r1.is_reverse else 'f'
		self.direction += 'r' if r2.is_reverse else 'f'
		self.name = r1.qname
		self.unique_flag = is_unique(r1) * is_unique(r2)
		self.mapp_qul_flag = has_good_mq(r1) * has_good_mq(r2)

	# Constructor from all fields values
	def from_all_args(self, name, first_read_chr, second_read_chr, \
			begin, end, length, middle, direction, unique_flag, mapp_qul_flag):
		self.name = name
		self.first_read_chr = first_read_chr # First read chromosome
		self.second_read_chr = second_read_chr # Second read chromosome
		self.begin = begin # Fragment begin
		self.end = end # Fragment end
		self.length = length # Length
		self.middle = middle # Middle
		# Supported directions:
		# rf (reverse - forward)
		# fr (forward - reverse)
		# ff (forward - forward)
		# rr (reverse - reverse)
		self.direction = direction
		self.unique_flag = unique_flag
		self.mapp_qul_flag = mapp_qul_flag

	def from_string(self, str):
		(name, first_read_chr, second_read_chr, begin, end, length, middle, direction, uf, mqf) = str.rstrip('\r\n').split(';')
		self.from_all_args(name, first_read_chr, second_read_chr, \
			int(begin), int(end), int(length), int(middle), direction, int(uf), int(mqf))

	def to_string(self):
		return self.name + ';' + self.first_read_chr + ';' + self.second_read_chr + ';' + \
			str(self.begin) + ';' + str(self.end) + ';' + \
			str(self.length) + ';' + str(self.middle) + ';' + \
			self.direction + ';' + str(self.unique_flag) + ';' + str(self.mapp_qul_flag)

	# For human-readable printing
	def __str__(self):
		return 'Fragment ' + self.name + ': [' + str(self.begin) + ';' + str(self.end) + ']' + \
			   ', direction ' + self.direction + ', chromosomes (' + \
			   self.first_read_chr + ',' + self.second_read_chr +' unique_flag ='+str(self.unique_flag)+ ','\
			   + 'mapp_qul_flag = '+str(self.mapp_qul_flag)+ ')'


	def is_abnormal(self, smallest_normal, biggest_normal):
		return self.length < smallest_normal or self.length > biggest_normal
