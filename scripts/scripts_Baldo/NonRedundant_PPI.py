import math
import string
import sys
import os
import re

from SeqIO import *


def printfasta (out,name,seq):
   out.write(">%s\n"%(name))
   n=0;
   w=""
   for s in seq:
    w += s
    n+=1
    if n == 80:
     out.write("%s\n"%(w))
     n=0
     w=""
   if n>0: out.write("%s\n"%(w))

def main():

 if len(sys.argv)>3:
  fa=sys.argv[1]
  ed=open(sys.argv[2],"r")
  out=open(sys.argv[3],'w')
 elif len(sys.argv)==3:
  fa=sys.argv[1]
  ed=open(sys.argv[2],"r")
  out=sys.stdout
 else:
  print("Execution requires two input files: 'NonRedundant_PPI.py  FASTA PPI [output]'")
  sys.exit(0)


 ppi={}
 for line in ed:
  p=re.compile(r"^(\w+)\s*(\w+)",re.IGNORECASE)
  m=p.match(line)
  if m is None:
    pass
  else:
   word=set(m.group(0).strip().split())
   seq=[s for s in FASTA_iterator(fa) if s.identifier in word]
   if len(seq) == 2:
    ppi.setdefault(seq[0],set()).add(seq[1])
    ppi.setdefault(seq[1],set()).add(seq[0])
   elif len(seq) == 1:
    ppi.setdefault(seq[0],set()).add(seq[0])
   else:
    print("Error number of Sequences on the interaction %s cannot be handled: %d \n"%("".join(word),len(seq)))

 nr_seq=set()
 pairs=set()
 for s in ppi.keys():
  if s not in nr_seq:
   nr_seq.add(s)
   uE1=s
  else:
   es=[x for x in nr_seq if x==s]
   uE1=es[0]
  for p in ppi[s]:
   if p not in nr_seq:
    nr_seq.add(p)
    uE2=p
   else:
    es=[x for x in nr_seq if x==p]
    uE2=es[0]
   if (uE1,uE2) not in pairs and (uE2,uE1) not in pairs:
    pairs.add((uE1,uE2))


 out.write("# SEQUENCES\n")
 out.write("#--------------------------------------------\n")
 for ns in nr_seq:
  printfasta(out,ns.identifier,ns.sequence)

 out.write("#--------------------------------------------\n")
 out.write("# INTERACTIONS\n")
 out.write("#--------------------------------------------\n")

 for p in pairs:
  out.write("%s::%s\n"%(p[0].get_identifier(),p[1].get_identifier()))


if  __name__ == "__main__":
 main()
