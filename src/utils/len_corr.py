a=open('/Users/dasha/PhD/data/length.txt','r')
c=[3240,5327]
new = open('/Users/dasha/PhD/data/len_corr.txt','w')
length = a.readline()[1:-1]
length = length.split(',')
for i in length:
	if 3240<int(i)<	5327:
		new.write(i+'\n')
new.close()
a.close()
