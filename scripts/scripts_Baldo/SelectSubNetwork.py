import sys
import argparse
import os


def main():

    options = parse_user_arguments()
    select(options)


def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Select a subnetwork with a particular set of nodes.",
        epilog      = "@oliva's lab 2014")
    parser.add_argument('-e','--network_input_file',dest='input_edge',action = 'store',default='edges.txt',
                        help = 'Input file of edges (default is edges.txt)')
    parser.add_argument('-n','--node_input_file',dest='input_node',action = 'store',default='nodes.txt',
                        help = 'Input file of nodes (default is nodes.txt)')
    parser.add_argument('-iformat','--input_format',dest='iformat',action = 'store',default='guild',
                        help = 'Format of input files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-oformat','--output_format',dest='oformat',action = 'store',default='guild',
                        help = 'Format of output files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-oe','--network_output_file',dest='network_output',action = 'store', default=sys.stdout,
                        help = 'Output file with edges (default is standard output)')
    parser.add_argument('-on','--nodes_output_file',dest='nodes_output',action = 'store', default=sys.stdout,
                        help = 'Output file with nodes (default is standard output)')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Flag to use verbose mode')
    options=parser.parse_args()

    return options

def fileExist (file):               #Checks if a file exists AND is a file
    return os.path.exists(file) and os.path.isfile(file)

def prune(ppi):
 pairs=set()
 for p,q in ppi:
   if (p,q) not in pairs and (q,p) not in pairs: pairs.add((p,q))
 return  pairs


def select(options):
    if not fileExist(options.input_edge):
     print("File with input network is missing\n")
     sys.exit(10)
    if not fileExist(options.input_node):
     print("File with selected nodes is missing\n")
     sys.exit(10)

    ppi=set()
    infonet={}
    fd=open(options.input_edge,"r")
    for line in fd:
      if options.iformat == 'netscore' :
       word=line.split()
       p=word[0]
       q=word[1]
       score=float(word[2])
      else:
       word=line.split()
       p=word[0]
       score=float(word[1])
       q=word[2]
      info=(score,"".join([str(word[i]) for i in range(3,len(word))]))
      infonet.setdefault((p,q),info)
      ppi.add((p,q))
      ppi.add((q,p))
    fd.close()

    ppi=prune(ppi)

    print("Num Original Interactions: %d\n"%(len(ppi)))

    nodes={}
    if fileExist(options.input_node):
     fd=open(options.input_node,"r")
     for line in fd:
      if options.iformat == 'netscore' :
       word=line.split()
       p=word[0]
       score=float(word[3])
      else:
       word=line.split()
       p=word[0]
       score=float(word[1])
      if nodes.has_key(p): nodes[p]+=score
      else: nodes[p]=score

#Check all nodes of the network
    realnode=set()
    for (p,q) in ppi:
     realnode.add(p)
     realnode.add(q)

    print("Total number of nodes in the network: %d\n"%(len(realnode)))
    print("Total number of nodes with information: %d\n"%(len(nodes)))


#Check all the nodes are in the network, otherwise prune the file
    n=0
    for p in nodes.keys():
        if p not in realnode:
         nodes[p]=None
         n+=1
    print("Total number of pruned nodes: %d\n"%(n))

    out=options.network_output
    if not options.network_output == sys.stdout: out=open(options.network_output,"w")

    for a,b in ppi:
     if a in nodes and b in nodes:
      if nodes[a] is not None and nodes[b] is not None:
        if infonet.has_key((a,b)):   score,info=infonet[(a,b)]
        elif infonet.has_key((b,a)): score,info=infonet[(b,a)]
        else:                        score,info=1.0,""
        if options.oformat == 'netscore' :
            out.write("\t{0}\t{1}\t{2:f}\t\t{3:s}\n".format(a,b,score,info))
        else:
            out.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
    out.close()

    out=options.nodes_output
    if not options.nodes_output == sys.stdout: out=open(options.nodes_output,"w")
    for p in nodes:
      if nodes[p] is not None:
        if options.oformat == 'netscore' :
         out.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(p,1.,1.,nodes[p]))
        else:
         out.write("{0} {1:10.5f}\n".format(p,nodes[p]))
    out.close()


if  __name__ == "__main__":
    main()


