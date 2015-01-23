new_f=open('chr21_gem.txt','w')
f = open('out35_hg18.gem','r')
flag=0
for line in f:
 if flag==1 and 'chr' not in line:
  new_f.write(line)
 elif flag==1 and 'chr' in line:
  break
 if 'chr21' in line:
  flag==1
 if 'chr' in line:
  print line
  
new_f.close()