"""
***
GML
***

Read graphs in GML format.
See http://www.infosun.fim.uni-passau.de/Graphlet/GML/gml-tr.html
for format specification.

Example graphs in GML format:
http://www-personal.umich.edu/~mejn/netdata/

"""
__author__ = """Aric Hagberg (hagberg@lanl.gov)"""
#    Copyright (C) 2008 by 
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    Distributed under the terms of the GNU Lesser General Public License
#    http://www.gnu.org/copyleft/lesser.html

__all__ = ['read_gml', 'parse_gml', 'write_gml']

import sys
import time
import networkx
from networkx.exception import NetworkXException, NetworkXError
from networkx.utils import _get_fh, is_string_like

	
def read_gml(path):
    """Read graph in GML format from path.
    Returns an Graph or DiGraph.

    This doesn't implement the complete GML specification for
    nested attributes for graphs, edges, and nodes. 

    """
    fh=_get_fh(path,mode='r')        
    G=parse_gml(fh)
    return G


def parse_gml(lines):
    """Parse GML format from string or iterable.
    Returns an Graph or DiGraph.

    This doesn't implement the complete GML specification for
    nested attributes for graphs, edges, and nodes. 
    """
    try:
        from pyparsing import ParseException
    except ImportError:
        raise ImportError, \
          "Import Error: not able to import pyparsing: http://pyparsing.wikispaces.com/"

    try:
        data = "".join(lines)
        gml = pyparse_gml()
        tokens =gml.parseString(data)
    except ParseException, err:
        print err.line
        print " "*(err.column-1) + "^"
        print err

    graph_attr=tokens.asDict()
    # determine directed or undirected and init corresponding NX class
    directed=graph_attr.get('directed',0)
    if directed==1:
        G=networkx.DiGraph()
    else:
        G=networkx.Graph()
    G.node_attr={} # store node attributes here
    G.graph_attr=graph_attr

    # first pass, nodes and labels
    label={}
    for item in tokens:
        if item[0]=='node':
            d=item.asDict()
            id=d['id']
            if 'label'in d:
                label[id]=d['label']
                del d['label']
            else:
                label[id]=id
                del d['id']
            G.add_node(label[id])
            if len(d)==0:
                d=None
            G.node_attr[label[id]]=d

    # second pass, edges            
    for item in tokens:
        if item[0]=='edge':
            d=item.asDict()
            source=d['source']
            target=d['target']
            del d['source']
            del d['target']
            if len(d)==0:
                d=None
            G.add_edge(label[source],label[target],d)
    return G
            
graph = None
def pyparse_gml():
    """pyparser tokenizer for GML graph format

    This doesn't implement the complete GML specification for
    nested attributes for graphs, edges, and nodes. 

    """  
    global graph
    
    try:
        from pyparsing import \
             Literal, CaselessLiteral,Word,\
             ZeroOrMore, Group, Dict, Optional, Combine,\
             ParseException, restOfLine, White, alphanums, nums,\
             OneOrMore,quotedString,removeQuotes,dblQuotedString
    except ImportError:
        raise ImportError, \
          "Import Error: not able to import pyparsing: http://pyparsing.wikispaces.com/"

    if not graph:
        graphkey = Literal("graph").suppress()
        lbrack = Literal("[").suppress()
        rbrack = Literal("]").suppress()
        pound = ("#")
        comment = pound + Optional( restOfLine )
        white = White(" \t\n")
        point = Literal(".")
        e = CaselessLiteral("E")
        integer = Word(nums).setParseAction(lambda s,l,t:[ int(t[0])])
        real = Combine( Word("+-"+nums, nums )+ 
                        Optional(point+Optional(Word(nums)))+
                        Optional(e+Word("+-"+nums, nums))).setParseAction(
                                        lambda s,l,t:[ float(t[0]) ])
        key=Word(alphanums)
        value=integer^real^Word(alphanums)^quotedString.setParseAction(removeQuotes)
        keyvalue = Dict(Group(key+OneOrMore(white).suppress()\
                   +value+OneOrMore(white).suppress()))
        node = Group(Literal("node") + lbrack + OneOrMore(keyvalue) + rbrack)
        edge = Group(Literal("edge") + lbrack + OneOrMore(keyvalue) + rbrack)
        graph = graphkey + lbrack + OneOrMore(edge|node|keyvalue) + rbrack
        graph.ignore(comment)
        
    return graph

def write_gml(G, path):
    """
    Write the graph G in GML format to the file or file handle path.

    Examples
    ---------

    >>> G=nx.path_graph(4)
    >>> nx.write_gml(G,"test.gml")

    path can be a filehandle or a string with the name of the file.

    >>> fh=open("test.gml",'w')
    >>> nx.write_gml(G,fh)

    Filenames ending in .gz or .bz2 will be compressed.

    >>> nx.write_gml(G,"test.gml.gz")

    
    The output file will use the default text encoding on your system.
    It is possible to write files in other encodings by opening
    the file with the codecs module.  See doc/examples/unicode.py
    for hints.

    >>> import codecs
    >>> fh=codecs.open("test.gml",'w',encoding='iso8859-1')# use iso8859-1
    >>> nx.write_gml(G,fh)

    GML specifications indicate that the file should only use
    7bit ASCII text encoding.iso8859-1 (latin-1). 

    Only a single level of attributes for graphs, nodes, and edges,
    is supported.

    """
    fh=_get_fh(path,mode='w')        
#    comments="#"
#    pargs=comments+" "+' '.join(sys.argv)
#    fh.write("%s\n" % (pargs))
#    fh.write(comments+" GMT %s\n" % (time.asctime(time.gmtime())))
#    fh.write(comments+" %s\n" % (G.name))

    # check for attributes or assign empty dict
    if hasattr(G,'graph_attr'):
        graph_attr=G.graph_attr
    else:
        graph_attr={}
    if hasattr(G,'node_attr'):
        node_attr=G.node_attr
    else:
        node_attr={}

    indent=2*' '
    count=iter(range(G.number_of_nodes()))
    node_id={}

    fh.write("graph [\n")
    # write graph attributes 
    for k,v in graph_attr.items():
        if is_string_like(v):
            v='"'+v+'"'
        fh.write(indent+"%s %s\n"%(k,v))
    # write nodes        
    for n in G:
        fh.write(indent+"node [\n")
        # get id or assign number
        if n in node_attr:
            nid=node_attr[n].get('id',count.next())
        else:
            nid=count.next()
        node_id[n]=nid
        fh.write(2*indent+"id %s\n"%nid)
        fh.write(2*indent+"label \"%s\"\n"%n)
        if n in node_attr:
          for k,v in node_attr[n].items():
              if is_string_like(v): v='"'+v+'"'
              if k=='id': continue
              fh.write(2*indent+"%s %s\n"%(k,v))
        fh.write(indent+"]\n")
    # write edges
    for e in G.edges_iter():
        if len(e)==3:
            u,v,d=e
            # try to guess what is on the edge and do something reasonable
            edgedata={}
            if d is None:  # no data 
                pass
            elif hasattr(d,'iteritems'): # try dict-like
                edgedata=dict(d.iteritems())
            elif hasattr(d,'__iter__'): # some kind of container
                edgedata['data']=d.__repr__()
            else: # something else - string, number
                edgedata['data']=d
        else:
            u,v=e
            edgedata={}
        fh.write(indent+"edge [\n")
        fh.write(2*indent+"source %s\n"%node_id[u])
        fh.write(2*indent+"target %s\n"%node_id[v])
        for k,v in edgedata.items():
            if k=='source': continue
            if k=='target': continue
            if is_string_like(v): v='"'+v+'"'
            fh.write(2*indent+"%s %s\n"%(k,v))
        fh.write(indent+"]\n")
    fh.write("]\n")
