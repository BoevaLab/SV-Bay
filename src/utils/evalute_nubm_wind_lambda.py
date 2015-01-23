gen_file = open('/Users/dasha/PhD/data/chrom_files_hg38/chr2.fa','r')
perc_at_arr = dict()
perc_gc_arr = dict()
for i in xrange(0,101,1):
	perc_at_arr[i] = 0
	perc_gc_arr[i] = 0
median = 2520
genome = ''
for i in gen_file:
	genome +=i[:-1]
print 'Finished reading fa file = ', len(genome)
genome = genome.lower()
for i in xrange(len(genome)-(median)-1):
	if not int(i)%1000000:
		print len(genome) - i,' last'
	win = genome[int(i):int(i+median)]
	if win.count('n')<=float(len(win))/5:
		perc_at =(win.count('a')+win.count('t'))/float(len(win)-win.count('n'))
		perc_at = perc_at*100
		perc_gc = 100 - perc_at
		perc_gc_arr[int(perc_gc)]+=1
		perc_at_arr[int(perc_at)]+=1
print 'Lenght of the genome = ',len(genome)
print 'Number of the windows with the given at rate:'
for i in perc_at_arr.keys():
	print 'Percentage = ',i,'Number of windows = ',perc_at_arr[i]
print 'Number of the windows with the given gc rate:'
for i in perc_gc_arr.keys():
	print 'Percentage = ',i,'Number of windows = ',perc_gc_arr[i]






