f = open('/Users/daria/Dropbox/data/sam_files/simulate/by_chr/chr1_filt.sam')
count=2000000
arr=[]
for i in xrange(count):
	arr.append(f.readline())
print 'finished reading!'
for ind in xrange(len(arr)):
	#print arr[ind].split()[0]
	if arr[ind].split()[0]=='read_4283_10539317':
		print 'Found!'
		for i in xrange(ind-10,ind+10):
			print arr[i]
		break