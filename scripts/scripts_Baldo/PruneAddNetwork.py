import sys
import argparse
import os


def main():

    options = parse_user_arguments()
    pruneadd(options)


def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Add a network or modify the scores by adding the scores of two files. The final network is prunned to skip redundant edges",
        epilog      = "@oliva's lab 2014")
    parser.add_argument('-e','--network_input_file',dest='input_edge',action = 'store',default='edges.txt',
                        help = 'Input file of edges (default is edges.txt)')
    parser.add_argument('-n','--node_input_file',dest='input_node',action = 'store',default='nodes.txt',
                        help = 'Input file of nodes (default is nodes.txt)')
    parser.add_argument('-ae','--add_edges_file',dest='add_edge',action = 'store',default='add_edges.txt',
                        help = 'File with the new edges or adding scores')
    parser.add_argument('-an','--add_nodes_file',dest='add_node',action = 'store',default='add_nodes.txt',
                        help = 'File with the new nodes or adding scores')
    parser.add_argument('-iformat','--input_format',dest='iformat',action = 'store',default='guild',
                        help = 'Format of input files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-oformat','--output_format',dest='oformat',action = 'store',default='guild',
                        help = 'Format of output files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-oe','--network_output_file',dest='network_output',action = 'store', default=sys.stdout,
                        help = 'Output file with edges (default is standard output)')
    parser.add_argument('-on','--nodes_output_file',dest='nodes_output',action = 'store', default=sys.stdout,
                        help = 'Output file with nodes (default is standard output)')
    parser.add_argument('-score','--non_seed_score',dest='score',action = 'store', default=0.01,type=float,
                        help = 'Score of non seed nodes (default is 0.01)')
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


def pruneadd(options):

    if not fileExist(options.input_edge):
     print("File with input network is missing\n")
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

    ppi2=set()
    infonet2={}
    if fileExist(options.add_edge):
     fd=open(options.add_edge,"r")
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
      infonet2.setdefault((p,q),info)
      ppi2.add((p,q))
      ppi2.add((q,p))
     fd.close()
     ppi2=prune(ppi2)
     print("Num Interactions in the adding set: %d\n"%(len(ppi2)))

    ppi3=set()
    for (p,q) in ppi2:
      if infonet2.has_key((p,q)):
       score,info2=infonet2[(p,q)]
      elif infonet2.has_key((q,p)):
       score,info2=infonet2[(q,p)]
      else:
       sys.stderr.write("Inconsistent added pair {0} - {1} without score\n".format(p,q))
       pass
      if (p,q) in ppi or (q,p) in ppi:
       if options.verbose: print("Found existing pair %s %s "%(p,q))
       if options.verbose: print("   Add Score %f "%(score))
       if infonet.has_key((p,q)):
         scr_old,info_old=infonet[(p,q)]
         score_new = scr_old + score
         if info2 != "": info_new  = info_old + "+" + info2
         else:info_new  = info_old
         infonet[(p,q)]=(score_new,info_new)
       elif infonet.has_key((q,p)):
         scr_old,info_old=infonet[(q,p)]
         score_new = scr_old + score
         if info2 != "": info_new  = info_old + "+" + info2
         else:info_new  = info_old
         infonet[(q,p)]=(score_new,info_new)
       else:
         sys.stderr.write("Inconsistent pair {0} - {1} without score\n".format(p,q))
         pass
       if options.verbose: print("   Final Score %f "%(score_new))
      else:
       info=(score,info2)
       infonet.setdefault((p,q),info)
       ppi3.add((p,q))
       ppi3.add((q,p))

    print("Num Added Interactions: %d\n"%(len(ppi3)))
    ppi.update(ppi3)
    print("Num Final Interactions(redundant): %d\n"%(len(ppi)))
    ppi=prune(ppi)
    print("Num Final Interactions(non-redundant): %d\n"%(len(ppi)))

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
    else:
     for p,q in ppi:
       nodes.setdefault(p,options.score)
       nodes.setdefault(q,options.score)

    if fileExist(options.add_node):
     fd=open(options.add_node,"r")
     for line in fd:
      if options.iformat == 'netscore' :
       word=line.split()
       p=word[0]
       score=float(word[3])
      else:
       word=line.split()
       p=word[0]
       score=float(word[1])
      if nodes.has_key(p):
        nodes[p]+=score
      else:
        nodes[p]=score

#Check all nodes of the network exist in nodes file otherwise add them
    realnode=set()
    for (p,q) in ppi:
     realnode.add(p)
     realnode.add(q)
    for p in realnode:
     if not nodes.has_key(p):
        nodes[p]=options.score

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



