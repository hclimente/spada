"""
***************
Graphviz AGraph
***************

Interface to pygraphviz AGraph class.

Usage 

 >>> G=nx.complete_graph(5)
 >>> A=nx.to_agraph(G)
 >>> H=nx.from_agraph(A)

Pygraphviz: http://networkx.lanl.gov/pygraphviz


"""
__author__ = """Aric Hagberg (hagberg@lanl.gov)"""
#    Copyright (C) 2004-2008 by 
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    Distributed under the terms of the GNU Lesser General Public License
#    http://www.gnu.org/copyleft/lesser.html

__all__ = ['from_agraph', 'to_agraph', 
           'write_dot', 'read_dot', 
           'graphviz_layout',
           'pygraphviz_layout']

import os
import sys
from networkx.utils import _get_fh
try:
    import pygraphviz
except ImportError:
    raise

def from_agraph(A,create_using=None):
    """Return a NetworkX Graph or DiGraph from a pygraphviz graph.


    >>> G=nx.complete_graph(5)
    >>> A=nx.to_agraph(G)
    >>> X=nx.from_agraph(A)

    The Graph X will have a dictionary X.graph_attr containing
    the default graphviz attributes for graphs, nodes and edges.

    Default node attributes will be in the dictionary X.node_attr
    which is keyed by node.

    Edge attributes will be returned as edge data in the graph X.

    If you want a Graph with no attributes attached instead of an XGraph
    with attributes use

    >>> G=nx.Graph(X)

    """
    import networkx
    if A.is_strict():
        multiedges=False
        selfloops=False
    else:
        multiedges=True
        selfloops=True
        
    if create_using is None:        
        if A.is_directed():
            if A.is_strict():
                create_using=networkx.DiGraph()
            else:
                create_using=networkx.MultiDiGraph()
        else:
            if A.is_strict():
                create_using=networkx.Graph()
            else:
                create_using=networkx.MultiGraph()

    # assign defaults        
    N=networkx.empty_graph(0,create_using)
    N.name=str(A)
    node_attr={}
    # add nodes, attributes to N.node_attr
    for n in A.nodes():
        N.add_node(n)
        node_attr[n]=n.attr

    # add edges, attributes attached to edge
    for e in A.edges():
        if len(e)==2:
            u,v=e
            N.add_edge(u,v)
        else:
            u,v,k=e
            attr=dict((k,v) for k,v in e.attr.items() if v!='')
            N.add_edge(u,v,attr)
        
    # add default attributes for graph, nodes, and edges       
    # hang them on N.graph_attr
    N.graph_attr={}
    N.graph_attr['graph']=A.graph_attr
    N.graph_attr['node']=A.node_attr
    N.graph_attr['edge']=A.edge_attr
    N.node_attr=node_attr
    return N        

def to_agraph(N, graph_attr=None, node_attr=None, strict=True):
    """Return a pygraphviz graph from a NetworkX graph N.

    If N is a Graph or DiGraph, graphviz attributes can
    be supplied through the arguments

    graph_attr:  dictionary with default attributes for graph, nodes, and edges
                 keyed by 'graph', 'node', and 'edge' to attribute dictionaries

    node_attr: dictionary keyed by node to node attribute dictionary

    If N has an dict N.graph_attr an attempt will be made first
    to copy properties attached to the graph (see from_agraph)
    and then updated with the calling arguments if any.

    """
    directed=N.directed
    strict=not N.multigraph # FIXME and N.number_of_selfloops()==0
    A=pygraphviz.AGraph(name=N.name,strict=strict,directed=directed)

    # default graph attributes            
    try:             
        A.graph_attr.update(N.graph_attr['graph'])
    except:
        pass
    try:
        A.graph_attr.update(graph_attr['graph'])
    except:
        pass
    # default node attributes            
    try:        
        A.node_attr.update(N.graph_attr['node'])
    except:
        pass
    try:
        A.node_attr.update(graph_attr['node'])
    except:
        pass
    # default edge attributes            
    try:        
        A.edge_attr.update(N.graph_attr['edge'])
    except:
        pass
    try:
        A.edge_attr.update(graph_attr['edge'])
    except:
        pass

    # add nodes
    for n in N:
        A.add_node(n)
        node=pygraphviz.Node(A,n)
        # try node attributes attached to graph
        try:
            if n in N.node_attr:
                node.attr.update(N.node_attr[n])
        except:
            pass
        # update with attributes from calling parameters
        try:
            if n in node_attr:
                node.attr.update(node_attr[n])
        except:
            pass

    # loop over edges
    for e in N.edges_iter():
        if len(e)==2:
            (u,v)=e
            d=None
        else:
            (u,v,d)=e

        if d is None: 
            # no data, just add edge
            A.add_edge(u,v)
        else: 
            if hasattr(d,"__iter__"):
                # edge data is dictionary-like, treat it as attributes
                # check for user assigned key
                if 'key' in d:
                    key=d['key']
                    del d['key']
                else:
                    key=str(d)
                data=d
            else:
                # edge data is some other object
                key=str(d)
                data={'data':key}
            A.add_edge(u,v,key=key,**data)
            edge=pygraphviz.Edge(A,u,v,key)
    return A

def write_dot(G,path):
    """Write NetworkX graph G to Graphviz dot format on path.

    Path can be a string or a file handle.
    """
    A=to_agraph(G)
    A.write(path)
    return

def read_dot(path,create_using=None):
    """Return a NetworkX XGraph or XdiGraph from a dot file on path.

    Path can be a string or a file handle.

    """
    A=pygraphviz.AGraph(file=path)
    return from_agraph(A)


def graphviz_layout(G,prog='neato',root=None, args=''):
    """
    Create layout using graphviz.
    Returns a dictionary of positions keyed by node.

    >>> G=nx.petersen_graph()
    >>> pos=nx.graphviz_layout(G)
    >>> pos=nx.graphviz_layout(G,prog='dot')
    
    This is a wrapper for pygraphviz_layout.

    """
    return pygraphviz_layout(G,prog=prog,root=root,args=args)

def pygraphviz_layout(G,prog='neato',root=None, args=''):
    """
    Create layout using pygraphviz and graphviz.
    Returns a dictionary of positions keyed by node.

    >>> G=nx.petersen_graph()
    >>> pos=nx.pygraphviz_layout(G)
    >>> pos=nx.pygraphviz_layout(G,prog='dot')
    
    """
    A=to_agraph(G)
    if root is not None:
        args+="-Groot=%s"%root
    A.layout(prog=prog,args=args)
    node_pos={}
    for n in G:
        node=pygraphviz.Node(A,n)
        try:
            xx,yy=node.attr["pos"].split(',')
            node_pos[n]=(float(xx),float(yy))
        except:
            print "no position for node",n
            node_pos[n]=(0.0,0.0)
    return node_pos

