f=open('/Users/daria/Dropbox/data/sam_files/chr20_filt.sam','r')
count = 0
for line in f:
	#print line.split()[6], line.split()[2]
	#count+=1
	if int(line.split()[1])&4 or int(line.split()[1])&8:
		count+=1
		print line
	if count == 20:
		break
f.close()
