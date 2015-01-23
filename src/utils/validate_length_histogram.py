##########################################################
# Small utility to validate length histogram
##########################################################

f = open('/opt/bio/data/length_histogram.txt', 'r')
hist_len = []

f.readline()
for line in f:
	(l, p) = line.split()
	hist_len.append((int(l),float(p)))
f.close()

prob_sum = 0.0
for (l, p) in hist_len:
	prob_sum += p

print 'Sum of probabilities: ' + str(prob_sum)
