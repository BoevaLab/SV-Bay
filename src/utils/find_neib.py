import random
import bisect
def sortByBeg(link):
	return int(link[0])

def sortByEnds(link):
	return int(link[1])
max_len_link = 80
chrom_len =5000
links=[]	
for l in xrange(50):
	i=random.randint(0,chrom_len- max_len_link)
	links.append([i,i+random.randint(10,max_len_link)])
links.sort()
links_sort_beg = [i[0] for i in links]
links_end =  sorted(links, key = sortByEnds)
links_sort_end = [i[1] for i in links_end]
print 'sortByBeg'
print links
print 'sortByEnds'
print links_end
for current_link in links:
	print 'Now for this link ',current_link
	flag = 0 
	if links[0]== current_link:
		print 'Flank region is [0,',current_link[0],']'
		flag =1
	else:
		print 'Start look for the Interception'
		for other_link in links[bisect.bisect_left(links_sort_beg,current_link[0]-max_len_link):bisect.bisect_left(links_sort_beg,current_link[0])+1]:
			if current_link == other_link:
				continue
			if current_link[0]<=other_link[1]:
				flag = 1
				print 'Interception!!!'
				print 'Initial lin is ', current_link,' neighbor is ', other_link
				break
		print 'Interception flag = ', flag
		if not flag:
			print 'then neighbor in the left is ', links_end[bisect.bisect_left(links_sort_end, current_link[0])-1]
	flag = 0
	print 'Look for right neighbor'
	if links_end[-1]==current_link:
		print 'Last link'
	else:
		left_boundaries = min(len(links),bisect.bisect_left(links_sort_beg,current_link[1]))
		print 'left_boundaries = ',left_boundaries
		for other_link in links[bisect.bisect_left(links_sort_beg,current_link[1]-max_len_link):left_boundaries]:
			print other_link	
			if other_link[1]>current_link[1]:
				flag = 1
				print 'Interception!!!'
				print 'Initial link is ', current_link,' neighbor is ', other_link
				break
		if not flag:
			left_boundaries = min(len(links)-1,bisect.bisect_left(links_sort_beg, current_link[1]))
			print 'then neighbor in the right is ', links_end[left_boundaries]






