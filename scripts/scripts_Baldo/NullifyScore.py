import sys

def main():

 if len(sys.argv)>3:
  input=sys.argv[1]
  out=open(sys.argv[2],'w')
  formated=sys.argv[3]
 elif len(sys.argv)>2:
  input=sys.argv[1]
  out=open(sys.argv[2],'w')
  formated='guild'
 else:
  input=sys.argv[1]
  out=sys.stdout
  formated='guild'

 fd=open(input,'r')
 for line in fd:
  word=line.strip().split()
  p=word[0]
  if formated == 'netscore' :
   out.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(p,1.,1.,0.))
  else:
   out.write("{0} {1:10.5f}\n".format(p,0.01))

 fd.close()
 out.close()

if  __name__ == "__main__":
    main()


