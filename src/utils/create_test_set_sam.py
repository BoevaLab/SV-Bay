f=open('chr20.sam','r')
#new_f=('chr20_test.sam','w')
#for line in f:
#	if line[0]=='@':
#		new_f.write(line)
#	else:
#		line1=line.split()
#		if 13000000<=line1[1]<=15010000:
#			new_f.write(line)
#new_f.close()
#f.close()
count=0
all_lines= 0
qual=0
unn_big_qual=0
unn_law_qual=0
for line in f:
	if line[0]=='@':
		continue
	else:
		all_lines+=1
		line1 = line.split()
		if  int(line1[4])>=28:
			qual+=1
		if 'XT:A:U' in line and int(line1[4])>=28:
			count+=1
		elif 'XT:A:U' not in line and int(line1[4])>=28:
			unn_big_qual+=1
		elif 'XT:A:U' not in line and int(line1[4])<28:
			unn_law_qual+=1
		#line1=line.split()
		#if 13000000<=line1[1]<=15010000:
		#	new_f.write(line)
print 'all_lines = ',all_lines
print 'super mapped reads!!! = ',count
print 'qual = ',qual
print 'unn_big_qual = ',unn_big_qual
print 'unn_law_qual = ',unn_law_qual
f.close()