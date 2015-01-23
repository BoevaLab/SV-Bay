f = open('/Users/dasha/Dropbox/GCAss_DB/contigs23.gem.mappability','r')
path = '/Users/dasha/Dropbox/GCAss_DB/gem_files/'
flag=0
for line in f:
	if 'NODE' in line:
		print 'Now for this chromosome ',line[1:-1]
		if flag:
			curr_file.close()
		curr_file = open(path+line[1:-1]+'_gem.txt','w')
		flag = 1
	elif flag:
		curr_file.write(line)
curr_file.close()