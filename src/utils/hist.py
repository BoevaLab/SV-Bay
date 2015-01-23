f = open('/Users/daria/Dropbox/data/length.txt','r')
normal_lengths = f.readline()[1:-1]
normal_lengths= normal_lengths.split(',')
normal_lengths = map(int, normal_lengths)
normal_lengths.sort()
print normal_lengths[1:1000	]
k=0
bins = set(normal_lengths)
bins_sort=[]
for i in bins:
	bins_sort.append(i)
bins_sort.sort()
historam_file = open('historam.txt','w')
for bin in bins_sort:
	historam_file.write(str(bin) + ' '+str(float(normal_lengths.count(bin))/len(normal_lengths))+'\n')
	print str(bin) + ' '+str(float(normal_lengths.count(bin))/len(normal_lengths))
historam_file.close()
f.close()