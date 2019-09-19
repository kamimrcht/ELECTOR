#!/usr/bin/python2.7
import sys
import os

if len(sys.argv)>=2:
    in_filename = sys.argv[1]

else:
    print("Remove quality information and separate the reads from their names for performance reasons")
    print("usage: ./ConverToPacBio_q2a.py input_filename")
    print("or python ConverToPacBio_q2a.py input_filename")
    sys.exit(1)
    
################################################################

f = open(in_filename,'r')
l = f.readline()

o = open('LR.fasta','w')
i = 0
while l:
    if l[0] !='>':
        print ("Err: invalid LR fastq format")
        exit(1)
    temp = l[1:]
    #length = int(temp.split('=')[1])
    name = "pacbio_LR_"+str(i)
    index = str(i)
    
    l=f.readline()
    ll = len(l)-1
    
    #if(length != ll):
    #    print "ERROR",i
    
    s = '>'+name+'/'+index+'/1_'+str(ll)+'\n'
    o.write(s)
    o.write(l)

    i = i + 1
    
    #l=f.readline()
    #l=f.readline()
    
    l=f.readline()
f.close()
o.close()

print (i)
