f=open('/Users/dasha/PhD/data/simulate_sam/chrX_sim.sam','r')
gem = open('/Users/dasha/PhD/data/gem_files_hg38/chrX_gem.txt')
#f_AM = ('AM_chr20','w')

gem_line = ''
for line in gem:
	if '~' not in line:
		if line[-1] == '\n':
			gem_line += line[:-1]
		else:
			gem_line += line
gem.close()
count=0
all_lines= 0
first_mapp = 0
second_mapp = 0
two_side_mapp =0
count_am =0
while  True:
	count+=1
	line1=f.readline()
	if not count%100000:
		print count ,' lines'
	if not line1:
		break
	if line1[0]=='@':
		continue
	else:
		line2 = f.readline().split()
		line1 = line1.split()
		all_lines+=1
		if len(line1)<14:
			continue
		AM = int(line1[14].split(':')[2])
		if AM>=20:
			count_am+=1
			read_1_beg = int(line1[3])
			read_2_beg = int(line1[7])
			if gem_line[read_1_beg]=='!' and gem_line[read_2_beg]=='!':
				two_side_mapp+=1
			elif gem_line[read_1_beg]!='!' and gem_line[read_2_beg]=='!':
				second_mapp+=1
			elif gem_line[read_1_beg]=='!' and gem_line[read_2_beg]!='!':
				first_mapp+=1

		#if AM==37:
		#	count_am+=1
		#if 'XT:A:U' in line1 and 'XT:A:U' in line2 :
			#f_AM.write(str(AM)+'\n')
		#	arr_AM.append(AM)
		#	if AM == 37:
		#		count_am_xa+=1
		#	elif 27<AM<37:
		#		count_what+=1
		#if line1[2]==line2[2]=='chr20':
			#if 'XT:A:U' in line1:
		#	if int(line1[4])>=29:
		#		one_side_mapp+=1
		#		if int(line2[4])<29 and line2[5]!='*':
		#			one_side_unmapp+=1
		#		if line2[5]=='*':
		#			one_side_non_mapp+=1
				#if 'XT:A:U' in line2: 
				#	two_side_mapp+=1
		#line1=line.split()
		#if 13000000<=line1[1]<=15010000:
		#	new_f.write(line)
#bins = set(arr_AM)
#len_AM = len(arr_AM)
#for i in bins:
#	print str(i)+' ' +str(arr_AM.count(i)/float(len_AM))
print 'All fragments with AM>20 count_am = ',count_am
print 'All fragments with both mappable reads = ',two_side_mapp/float(two_side_mapp+first_mapp+second_mapp)
print 'Only first read on the mappable position = ',first_mapp/float(two_side_mapp+first_mapp+second_mapp)
print 'Only second read on the mappable position = ',second_mapp/float(two_side_mapp+first_mapp+second_mapp)

#print 'all_lines = ',all_lines
#print 'one_side_unmapp = ',one_side_unmapp
#print 'one_side_mapp = ',one_side_mapp
#print 'one_side_non_mapp = ',one_side_non_mapp
#print 'one_side_unmapp/one_side_mapp = ',one_side_unmapp/float(one_side_mapp)
#print 'one_side_non_mapp/one_side_mapp = ',float(one_side_non_mapp)/one_side_mapp
#print 'two_side_mapp',two_side_mapp

f.close()