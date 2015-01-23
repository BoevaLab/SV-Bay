f = open('clusters_7_21.txt','r')
new_file =open('clusters_7_21_new_names.txt','w')
count = 0
for line in f:
	line1 = line.split(';')
	new_name = 'cl_'+line[8]+'_'+line[9]+'_'+count
	if line[8]!=line[9]:
		new_name+='_tr_'+line[7]
	else:
		new_name+='_'+line[7]
	new_line = new_name+';'
	for i in line1[1:]:
		new_line+='_'+i
	count+=1
	new_file.write(new_line)
new_file.close()
