from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np




# read data from a text file. One number per line
arch = "/Users/dasha/PhD/data/len_corr.txt"
datos = []
for item in open(arch,'r'):
    item = item.strip()
    if item != '':
        try:
            datos.append(float(item))
        except ValueError:
            pass

# best fit of data
(mu, sigma) = norm.fit(datos)
print mu, sigma

# the histogram of the data
print len(set(datos))
n, bins, patches = plt.hist(datos, len(set(datos)), normed=1, facecolor='green', alpha=0.75)

# add a 'best fit' line
a=mu-(12**0.5)*0.5*sigma
b=mu+(12**0.5)*0.5*sigma
p = 1/(b-a)
print a,b,p

y = mlab.normpdf( bins, mu, sigma)
print 'y = ', y[1000:1200]
print 'n = ',n[1000:1200]
#a = map(int, y)
#print a
print bins[0:100]
for i in bins[0:120]:
	print round(i)

l = plt.plot(bins, y, 'r--')
s = plt.plot([a,b],[p,p], 'b--',[a,a],[0,p],'b--',[b,b],[0,p],'b--')
print mlab.normpdf(p,mu,sigma)
#plot
plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ the \ lengths \ of \ the \ fragments :}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
plt.grid(True)

#plt.show()