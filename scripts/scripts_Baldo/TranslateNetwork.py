import sys
import argparse
import os


def main():

    options = parse_user_arguments()
    translate(options)


def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Translate a bPPI network into your selected codes",
        epilog      = "@oliva's lab 2014")
    parser.add_argument('-i','--network_input_file',dest='input_edge',action = 'store',default='edges.txt',
                        help = 'Input file of network (default is edges.txt)')
    parser.add_argument('-n','--nodes_input_file',dest='input_node',action = 'store',default='nodes.txt',
                        help = 'Input file of nodes (default is nodes.txt)')
    parser.add_argument('-trans','--translation_of_nodes_file',dest='translation_file',action = 'store',default='translation_nodes.txt',
                        help = 'File with the translation of codes from BIANA to the selected type for all nodes')
    parser.add_argument('-iformat','--input_format',dest='iformat',action = 'store',default='guild',
                        help = 'Format of input files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-oformat','--output_format',dest='oformat',action = 'store',default='guild',
                        help = 'Format of output files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-oe','--output_file_edges',dest='output_edge',action = 'store', default=sys.stdout,
                        help = 'Output file with edges(default is standard output)')
    parser.add_argument('-on','--output_file_nodes',dest='output_node',action = 'store', default=sys.stdout,
                        help = 'Output file with nodes(default is standard output)')
    options=parser.parse_args()

    return options

def fileExist (file):               #Checks if a file exists AND is a file
    return os.path.exists(file) and os.path.isfile(file)

def translate(options):
    if not fileExist(options.translation_file):
     print("File with translation is missing\n")
     sys.exit(10)
    if not fileExist(options.input_edge):
     print("File with input network is missing\n")
     sys.exit(10)

    ft=open(options.translation_file,"r")
    new={}
    for line in ft:
     word=line.strip().split("\t")
     node=word[0]
     new.setdefault(node,set())
     for trans in word[1].split("'"):
       if trans != "" and trans != ",":
         name="_".join([str(x) for x in trans.split()])
         new[node].add(name)
    ft.close()

    fd=open(options.input_edge,"r")
    out_network=options.output_edge
    if not options.output_edge == sys.stdout: out_network=open(options.output_edge,"w")

    nodes=set()
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
      info="".join([str(word[i]) for i in range(3,len(word))])
      for a in new[p]:
       nodes.add(a)
       for b in new[q]:
          nodes.add(b)
          if options.oformat == 'netscore' :
            out_network.write("\t{0}\t{1}\t{2:f}\t\t{3:s}\n".format(a,b,score,info))
          else:
            out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
    fd.close()
    if not options.output_edge == sys.stdout: out_network.close()

    out_network=options.output_node
    if not options.output_node == sys.stdout: out_network=open(options.output_node,"w")

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
      for a in new[p]:
          if options.oformat == 'netscore' :
           out_network.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(a,1.,1.,score))
          else:
           out_network.write("{0} {1:10.5f}\n".format(a,score))
     fd.close()
    else:
     for a in nodes:
          if options.oformat == 'netscore' :
           out_network.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(a,1.,1.,0.))
          else:
           out_network.write("{0} {1:10.5f}\n".format(a,0.01))

    out_network.close()

if  __name__ == "__main__":
    main()



