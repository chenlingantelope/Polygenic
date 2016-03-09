from __future__ import print_function
import os
import pandas
import numpy as np
import sys

#write this as an argument
seqlength=int(sys.argv[1])

os.system("egrep -x '[0-1]+' out2.msms > haplo.txt")
os.system("sed -i -e 's/0/A/g' haplo.txt")
os.system("sed -i -e 's/1/T/g' haplo.txt")
nind=int(str.split(os.popen('wc -l haplo.txt').read())[0])

segsites=os.popen('grep \'positions\' out2.msms').read()
pos=map(float,segsites.split(': ')[1].split(' \n')[0].split(' '))
pos=np.multiply(pos,seqlength).astype(int)

f=open('snp','w')
for i in range(len(pos)):
    print('  - rs'+str(i)+': '+str(pos[i]),file=f)

f.close()


f = open('haplo.txt', 'r')
haplo=f.read()
haplotype=haplo.split("\n")
haplotype=haplotype[0:len(haplotype)-1]

f=open('haplotype',"w")
print('phased_haplotypes:', file=f)
for ind in range(len(haplotype)/2):
    print("  - NA1111",ind,"_c1: ",haplotype[ind*2],sep='',file=f)
    print("  - NA1111",ind,"_c2: ",haplotype[ind*2+1],sep='',file=f)

f.close()


nsites=len(haplotype[0])


f=open("Ans_Der"+'.txt',"w")
g=open('snp','w')
for site in range(nsites):
    print("rs"+str(site),'A','T',sep='\t',file=f)
    print('  - rs'+str(site)+': '+ str(pos[site]),file=g)

f.close()
g.close()

f=open('sample'+'.txt',"w")
for ind in range(len(haplotype)/2):
    print("NA1111",ind,sep='',file=f)

f.close()

f=open('header','w')
print('---',file=f)
print('pop: '+'msms',file=f)
print('build: msms',file=f)
#change name to include simulation script and 
print('hapmap_release: bla',file=f)
print('start: 0',file=f)
print('end: '+str(seqlength),file=f)
print('snps:',file=f)
f.close()


os.system('cat header snp haplotype > Haplotype'+'.txt')
os.system('rm header snp haplotype haplo.txt haplo.txt-e')
os.system('../nSL_03032016/nSL -samfile sample.txt -hapfile Haplotype.txt -adfile Ans_Der.txt -maxLen '+str(seqlength)+' > out2.nSL')
os.system('rm sample.txt Haplotype.txt Ans_Der.txt')
