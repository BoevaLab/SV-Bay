f=open('/Users/daria/Dropbox/data/sam_files/simulate/by_chr/chr2_filt.sam','r')
new_f= open('/Users/daria/Dropbox/data/sam_files/simulate/by_chr/chr2_filt_corr.sam','w')
reads=[]
print 'write to this file chr1_corr.sam'
flag=0
count =0
count_flag=0
for line in f:
	
	count+=1
	reads.append(line)
	if flag:
		#print 'New couple'
		#print reads[0]
		#print reads[1]
		count_flag+=1
		if count_flag<=20:
			print line
		else:
			break
	if not count%1000000:
		print 'Already did ',count
	if len(reads)== 2 and not flag:
		if reads[0].split()[0]==reads[1].split()[0]:
			new_f.write(reads[0])
			new_f.write(reads[1])
			reads=[]
		else:
			print 'That is bad!!!!'
			print reads[0]
			print reads[1]
			print count

			reads.remove(reads[0])
			flag=1

f.close()
new_f.close()