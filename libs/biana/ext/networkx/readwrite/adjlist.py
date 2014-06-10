"""
***************
Adjacency Lists
***************

Read and write NetworkX graphs as adjacency lists.

Note that NetworkX graphs can contain any hashable Python object as
node (not just integers and strings).  So writing a NetworkX graph
as a text file may not always be what you want: see write_gpickle
and gread_gpickle for that case.

This module provides the following :

Adjacency list with single line per node:
Useful for connected or unconnected graphs without edge data.

    write_adjlist(G, path)
    G=read_adjlist(path)

Adjacency list with multiple lines per node:
Useful for connected or unconnected graphs with or without edge data.

    write_multiline_adjlist(G, path)
    read_multiline_adjlist(path)

"""
__author__ = """Aric Hagberg (hagberg@lanl.gov)\nDan Schult (dschult@colgate.edu)"""
#    Copyright (C) 2004-2008 by 
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    Distributed under the terms of the GNU Lesser General Public License
#    http://www.gnu.org/copyleft/lesser.html

__all__ = ['read_multiline_adjlist', 'write_multiline_adjlist',
           'read_adjlist', 'write_adjlist']

import cPickle 
import codecs
import locale
import string
import sys
import time

from networkx.utils import is_string_like,_get_fh
import networkx

def write_multiline_adjlist(G, path, delimiter=' ', comments='#'):
    """
    Write the graph G in multiline adjacency list format to the file
    or file handle path.

    See read_multiline_adjlist for file format details.

    Examples
    --------

    >>> G=nx.path_graph(4)
    >>> nx.write_multiline_adjlist(G,"test.adjlist")

    path can be a filehandle or a string with the name of the file.

    >>> fh=open("test.adjlist",'w')
    >>> nx.write_multiline_adjlist(G,fh)

    Filenames ending in .gz or .bz2 will be compressed.

    >>> nx.write_multiline_adjlist(G,"test.adjlist.gz")

    The file will use the default text encoding on your system.
    It is possible to write files in other encodings by opening
    the file with the codecs module.  See doc/examples/unicode.py
    for hints.

    >>> import codecs
    >>> fh=codecs.open("test.adjlist",'w',encoding='utf=8') # utf-8 encoding
    >>> nx.write_multiline_adjlist(G,fh)

    """
    fh=_get_fh(path,mode='w')        
    pargs=comments+" "+string.join(sys.argv,' ')
    fh.write("%s\n" % (pargs))
    fh.write(comments+" GMT %s\n" % (time.asctime(time.gmtime())))
    fh.write(comments+" %s\n" % (G.name))

    def make_str(t):
        if is_string_like(t): return t
        return str(t)

    if G.directed:
        if G.multigraph:
            for s,nbrs in G.adjacency_iter():
                nbr_edges=[ (u,d) for u,dl in nbrs.iteritems() for d in dl]
                deg=len(nbr_edges)
                fh.write(make_str(s)+delimiter+"%i\n"%(deg))
                for u,d in nbr_edges:
                    if d is None:    
                        fh.write(make_str(u)+'\n')
                    else:   
                        fh.write(make_str(u)+delimiter+make_str(d)+"\n")
        else: # directed single edges
            for s,nbrs in G.adjacency_iter():
                deg=len(nbrs)
                fh.write(make_str(s)+delimiter+"%i\n"%(deg))
                for u,d in nbrs.iteritems():
                    if d is None:    
                        fh.write(make_str(u)+'\n')
                    else:   
                        fh.write(make_str(u)+delimiter+make_str(d)+"\n")
    else: #undirected
        if G.multigraph:
            seen={}  # helper dict used to avoid duplicate edges
            for s,nbrs in G.adjacency_iter():
                nbr_edges=[ (u,d) for u,dl in nbrs.iteritems() if u not in seen for d in dl]
                deg=len(nbr_edges)
                fh.write(make_str(s)+delimiter+"%i\n"%(deg))
                for u,d in nbr_edges:
                    if d is None:    
                        fh.write(make_str(u)+'\n')
                    else:   
                        fh.write(make_str(u)+delimiter+make_str(d)+"\n")
                seen[s]=1
        else: # undirected single edges
            seen={}  # helper dict used to avoid duplicate edges
            for s,nbrs in G.adjacency_iter():
                nbr_edges=[ (u,d) for u,d in nbrs.iteritems() if u not in seen]
                deg=len(nbr_edges)
                fh.write(make_str(s)+delimiter+"%i\n"%(deg))
                for u,d in nbr_edges:
                    if d is None:    
                        fh.write(make_str(u)+'\n')
                    else:   
                        fh.write(make_str(u)+delimiter+make_str(d)+"\n")
                seen[s]=1
            
def read_multiline_adjlist(path, comments="#", delimiter=' ',
                           create_using=None,
                           nodetype=None, edgetype=None):
    """Read graph in multi-line adjacency list format from path.

    Examples
    --------

    >>> G=nx.path_graph(4)
    >>> nx.write_multiline_adjlist(G,"test.adjlist")
    >>> G=nx.read_multiline_adjlist("test.adjlist")

    path can be a filehandle or a string with the name of the file.

    >>> fh=open("test.adjlist")
    >>> G=nx.read_multiline_adjlist(fh)

    Filenames ending in .gz or .bz2 will be compressed.

    >>> nx.write_multiline_adjlist(G,"test.adjlist.gz")
    >>> G=nx.read_multiline_adjlist("test.adjlist.gz")

    nodetype is an optional function to convert node strings to nodetype

    For example

    >>> G=nx.read_multiline_adjlist("test.adjlist", nodetype=int)

    will attempt to convert all nodes to integer type

    Since nodes must be hashable, the function nodetype must return hashable
    types (e.g. int, float, str, frozenset - or tuples of those, etc.) 

    edgetype is a function to convert edge data strings to edgetype

    >>> G=nx.read_multiline_adjlist("test.adjlist", edgetype=int)

    create_using is an optional networkx graph type, the default is
    Graph(), a simple undirected graph 

    >>> G=nx.read_multiline_adjlist("test.adjlist", create_using=nx.DiGraph())

    The comments character (default='#') at the beginning of a
    line indicates a comment line.

    The entries are separated by delimiter (default=' ').
    If whitespace is significant in node or edge labels you should use
    some other delimiter such as a tab or other symbol.  
    

    Example multiline adjlist file format

    No edge data::

     # source target for Graph or DiGraph
     a 2
     b
     c
     d 1
     e

     Wiht edge data::

     # source target for XGraph or XDiGraph with edge data
     a 2
     b edge-ab-data
     c edge-ac-data
     d 1
     e edge-de-data

    Reading the file will use the default text encoding on your system.
    It is possible to read files with other encodings by opening
    the file with the codecs module.  See doc/examples/unicode.py
    for hints.

    >>> import codecs
    >>> fh=codecs.open("test.adjlist",'r',encoding='utf=8') # utf-8 encoding
    >>> G=nx.read_multiline_adjlist(fh)
    """
    if create_using is None:
        G=networkx.Graph()
    else:
        try:
            G=create_using
            G.clear()
        except:
            raise TypeError("Input graph is not a networkx graph type")

    inp=_get_fh(path)        

    for line in inp:
#        if line.startswith("#") or line.startswith("\n"):
#            continue
#        line=line.strip() #remove trailing \n 
        line = line[:line.find(comments)].strip()
        if not line: continue
        try:
            (u,deg)=line.split(delimiter)
            deg=int(deg)
        except:
            raise TypeError("Failed to read node and degree on line (%s)"%line)
        if nodetype is not None:
            try:
                u=nodetype(u)
            except:
                raise TypeError("Failed to convert node (%s) to type %s"\
                      %(u,nodetype))
        G.add_node(u)
        for i in range(deg):
            line=inp.next().strip()
            vlist=line.split(delimiter)
            numb=len(vlist)
            if numb>0:
                v=vlist[0]
                if nodetype is not None:
                    try:
                        v=nodetype(v)
                    except:
                        raise TypeError("Failed to convert node (%s) to type %s"\
                                        %(v,nodetype))
            if numb==1:
                G.add_edge(u,v)
            elif numb==2:
                d=vlist[1]
                if edgetype is not None:
                    try:
                        d=edgetype(d)
                    except:
                        raise TypeError\
                              ("Failed to convert edge data (%s) to type %s"\
                                %(d, edgetype))
                G.add_edge(u,v,d)
            else:
                raise TypeError("Failed to read line: %s"%vlist)
    return G


def write_adjlist(G, path, comments="#", delimiter=' '):
    """Write graph G in single-line adjacency-list format to path.

    See read_adjlist for file format details.

    Examples
    --------

    >>> G=nx.path_graph(4)
    >>> nx.write_adjlist(G,"test.adjlist")

    path can be a filehandle or a string with the name of the file.

    >>> fh=open("test.adjlist",'w')
    >>> nx.write_adjlist(G, fh)

    Filenames ending in .gz or .bz2 will be compressed.

    >>> nx.write_adjlist(G, "test.adjlist.gz")

    The file will use the default text encoding on your system.
    It is possible to write files in other encodings by opening
    the file with the codecs module.  See doc/examples/unicode.py
    for hints.

    >>> import codecs

    fh=codecs.open("test.adjlist",encoding='utf=8') # use utf-8 encoding
    nx.write_adjlist(G,fh)

    Does not handle edge data. 
    Use 'write_edgelist' or 'write_multiline_adjlist'
    """
    fh=_get_fh(path,mode='w')        
    pargs=comments+" "+string.join(sys.argv,' ')
    fh.write("%s\n" % (pargs))
    fh.write(comments+" GMT %s\n" % (time.asctime(time.gmtime())))
    fh.write(comments+" %s\n" % (G.name))

    def make_str(t):
        if is_string_like(t): return t
        return str(t)
    directed=G.directed

    seen={}
    for s,nbrs in G.adjacency_iter():
        fh.write(make_str(s)+delimiter)
        for t,data in nbrs.iteritems():
            if not directed and t in seen: continue
            if G.multigraph:
                for d in data:
                    fh.write(make_str(t)+delimiter)
            else:
                fh.write(make_str(t)+delimiter)
        fh.write("\n")            
        if not directed: 
            seen[s]=1


def read_adjlist(path, comments="#", delimiter=' ',
                 create_using=None, nodetype=None):
    """Read graph in single line adjacency list format from path.

    Examples
    --------

    >>> G=nx.path_graph(4)
    >>> nx.write_adjlist(G, "test.adjlist")
    >>> G=nx.read_adjlist("test.adjlist")

    path can be a filehandle or a string with the name of the file.

    >>> fh=open("test.adjlist")
    >>> G=nx.read_adjlist(fh)

    Filenames ending in .gz or .bz2 will be compressed.

    >>> nx.write_adjlist(G, "test.adjlist.gz")
    >>> G=nx.read_adjlist("test.adjlist.gz")

    nodetype is an optional function to convert node strings to nodetype

    For example

    >>> G=nx.read_adjlist("test.adjlist", nodetype=int)

    will attempt to convert all nodes to integer type

    Since nodes must be hashable, the function nodetype must return hashable
    types (e.g. int, float, str, frozenset - or tuples of those, etc.) 

    create_using is an optional networkx graph type, the default is
    Graph(), an undirected graph. 

    >>> G=nx.read_adjlist("test.adjlist", create_using=nx.DiGraph())

    Does not handle edge data: use 'read_edgelist' or 'read_multiline_adjlist'

    The comments character (default='#') at the beginning of a
    line indicates a comment line.

    The entries are separated by delimiter (default=' ').
    If whitespace is significant in node or edge labels you should use
    some other delimiter such as a tab or other symbol.  

    Sample format::

     # source target
     a b c
     d e

    """
    if create_using is None:
        G=networkx.Graph()
    else:
        try:
            G=create_using
            G.clear()
        except:
            raise TypeError("Input graph is not a networkx graph type")

    fh=_get_fh(path)        

    for line in fh.readlines():
        line = line[:line.find(comments)].strip()
        if not len(line): continue
#        if line.startswith("#") or line.startswith("\n"):
#            continue
#        line=line.strip() #remove trailing \n 
        vlist=line.split(delimiter)
        u=vlist.pop(0)
        # convert types
        if nodetype is not None:
            try:
                u=nodetype(u)
            except:
                raise TypeError("Failed to convert node (%s) to type %s"\
                                %(u,nodetype))
        G.add_node(u)
        try:
            vlist=map(nodetype,vlist)
        except:
            raise TypeError("Failed to convert nodes (%s) to type %s"\
                            %(','.join(vlist),nodetype))
        for v in vlist:
            G.add_edge(u,v)
    return G


