import random
import matplotlib.pyplot as plt
 
f = open('/Users/dasha/PhD/data/donor/donor_RF_70_4283.insertsize.distr')
f.readline()
f.readline()
x=[]
y=[]
for line in f:
	y.append(int(line.split()[1]))
	x.append(int(line.split()[0]))
	
f.close()
plt.bar(x,y)
plt.show()
