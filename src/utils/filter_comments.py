f =open('/Users/dasha/PhD/data/donor/comments_for_donor.fa.txt','r')
arr =[]
new_for_valid = open('hg_38_hg19.txt','w')
for line in f:
	line_sp = line.split('	')
	print line_sp
	if line in arr:
		continue
	else:
		if line_sp[2]!='Duplication' and int(line_sp[4])-int(line_sp[3])>1:
			new_for_valid.write('\n'+line_sp[0].split('_')[1]+':'+line_sp[3]+'-'+line_sp[4])
			arr.append(line)
f.close()
new = open('sv_sim_filt.txt','w')
for i in arr:
	new.write(i)
new.close()
new_for_valid.close()
