"""
**********
Matplotlib
**********

Draw networks with matplotlib (pylab).

References:
 - matplotlib:     http://matplotlib.sourceforge.net/
 - pygraphviz:     http://networkx.lanl.gov/pygraphviz/

"""
__author__ = """Aric Hagberg (hagberg@lanl.gov)"""
#    Copyright (C) 2004-2008 by 
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    Distributed under the terms of the GNU Lesser General Public License
#    http://www.gnu.org/copyleft/lesser.html

__all__ = ['draw',
           'draw_networkx',
           'draw_networkx_nodes',
           'draw_networkx_edges',
           'draw_networkx_labels',
           'draw_circular',
           'draw_random',
           'draw_spectral',
           'draw_spring',
           'draw_shell',
           'draw_graphviz']


import networkx

import sys

try:
    import matplotlib
    import matplotlib.cbook as cb
    from matplotlib.colors import colorConverter,normalize,Colormap
    from matplotlib.collections import LineCollection
    from matplotlib.numerix import sin, cos, pi, sqrt, arctan2, asarray
    from matplotlib.numerix.mlab import amin, amax, ravel
    import matplotlib.pylab
except ImportError:
    raise ImportError, "Import Error: not able to import matplotlib."
except RuntimeError:
    pass # unable to open display

def draw(G, pos=None, ax=None, hold=None, **kwds):
    """Draw the graph G with matplotlib (pylab).

    This is a pylab friendly function that will use the
    current pylab figure axes (e.g. subplot).

    pos is a dictionary keyed by vertex with a two-tuple
    of x-y positions as the value.
    See networkx.layout for functions that compute node positions.

    Usage:

    >>> from networkx import *
    >>> G=dodecahedral_graph()
    >>> draw(G)
    >>> pos=graphviz_layout(G)
    >>> draw(G,pos)
    >>> draw(G,pos=spring_layout(G))

    Also see doc/examples/draw_*

    :Parameters:

      - `nodelist`: list of nodes to be drawn (default=G.nodes())
      - `edgelist`: list of edges to be drawn (default=G.edges())
      - `node_size`: scalar or array of the same length as nodelist (default=300)
      - `node_color`: single color string or numeric/numarray array of floats (default='r')
      - `node_shape`: node shape (default='o'), or 'so^>v<dph8' see pylab.scatter
      - `alpha`: transparency (default=1.0) 
      - `cmap`: colormap for mapping intensities (default=None)
      - `vmin,vmax`: min and max for colormap scaling (default=None)
      - `width`: line width of edges (default =1.0)
      - `edge_color`: scalar or array (default='k')
      - `edge_cmap`: colormap for edge intensities (default=None) 
      - `edge_vmin,edge_vmax`: min and max for colormap edge scaling (default=None)
      - `style`: edge linestyle (default='solid') (solid|dashed|dotted,dashdot)
      - `labels`: dictionary keyed by node of text labels (default=None)
      - `font_size`: size for text labels (default=12)
      - `font_color`: (default='k')
      - `font_weight`: (default='normal')
      - `font_family`: (default='sans-serif')
      - `ax`: matplotlib axes instance

    for more see pylab.scatter

    NB: this has the same name as pylab.draw so beware when using

    >>> from networkx import *

    since you will overwrite the pylab.draw function.

    A good alternative is to use

    >>> import pylab as P
    >>> import networkx as NX
    >>> G=NX.dodecahedral_graph()

    and then use

    >>> NX.draw(G)  # networkx draw()

    and
    >>> P.draw()    # pylab draw()

    """
    if pos is None:
        pos=networkx.drawing.spring_layout(G) # default to spring layout

    cf=matplotlib.pylab.gcf()
    cf.set_facecolor('w')
    if ax is None:
        if cf._axstack() is None:
            ax=cf.add_axes((0,0,1,1))
        else:
            ax=cf.gca()

 # allow callers to override the hold state by passing hold=True|False
    b = matplotlib.pylab.ishold()
    h = kwds.pop('hold', None)
    if h is not None:
        matplotlib.pylab.hold(h)
    try:
        draw_networkx(G,pos,ax=ax,**kwds)
        ax.set_axis_off()
        matplotlib.pylab.draw_if_interactive()
    except:
        matplotlib.pylab.hold(b)
        raise
    matplotlib.pylab.hold(b)
    return


def draw_networkx(G, pos, with_labels=True, **kwds):
    """Draw the graph G with given node positions pos

    Usage:

    >>> from networkx import *
    >>> import pylab as P
    >>> ax=P.subplot(111)
    >>> G=dodecahedral_graph()
    >>> pos=spring_layout(G)
    >>> draw_networkx(G,pos,ax=ax)

    This is same as 'draw' but the node positions *must* be
    specified in the variable pos.
    pos is a dictionary keyed by vertex with a two-tuple
    of x-y positions as the value.
    See networkx.layout for functions that compute node positions.

    An optional matplotlib axis can be provided through the
    optional keyword ax.

    with_labels contols text labeling of the nodes

    Also see:

    draw_networkx_nodes()
    draw_networkx_edges()
    draw_networkx_labels()
    """
    from matplotlib.pylab import draw_if_interactive 
    node_collection=draw_networkx_nodes(G, pos, **kwds)
    edge_collection=draw_networkx_edges(G, pos, **kwds) 
    if with_labels:
        draw_networkx_labels(G, pos, **kwds)
    draw_if_interactive()

def draw_networkx_nodes(G, pos,
                        nodelist=None,
                        node_size=300,
                        node_color='r',
                        node_shape='o',
                        alpha=1.0,
                        cmap=None,
                        vmin=None,
                        vmax=None, 
                        ax=None,
                        linewidths=None,
                        **kwds):
    """Draw nodes of graph G

    This draws only the nodes of the graph G.

    pos is a dictionary keyed by vertex with a two-tuple
    of x-y positions as the value.
    See networkx.layout for functions that compute node positions.

    nodelist is an optional list of nodes in G to be drawn.
    If provided only the nodes in nodelist will be drawn.
    
    see draw_networkx for the list of other optional parameters.

    """
    if ax is None:
        ax=matplotlib.pylab.gca()

    if nodelist is None:
        nodelist=G.nodes()

    if not nodelist or len(nodelist)==0:  # empty nodelist, no drawing
        return None 

    try:
        xy=asarray([pos[v] for v in nodelist])
    except KeyError,e:
        raise networkx.NetworkXError('Node %s has no position.'%e)
    except ValueError:
        raise networkx.NetworkXError('Bad value in node positions.')


    node_collection=ax.scatter(xy[:,0], xy[:,1],
                               s=node_size,
                               c=node_color,
                               marker=node_shape,
                               cmap=cmap, 
                               vmin=vmin,
                               vmax=vmax,
                               alpha=alpha,
                               linewidths=linewidths)
                               
    matplotlib.pylab.sci(node_collection)
    node_collection.set_zorder(2)            
    return node_collection


def draw_networkx_edges(G, pos,
                        edgelist=None,
                        width=1.0,
                        edge_color='k',
                        style='solid',
                        alpha=1.0,
                        edge_cmap=None,
                        edge_vmin=None,
                        edge_vmax=None, 
                        ax=None,
                        arrows=True,
                        **kwds):
    """Draw the edges of the graph G

    This draws only the edges of the graph G.

    pos is a dictionary keyed by vertex with a two-tuple
    of x-y positions as the value.
    See networkx.layout for functions that compute node positions.

    edgelist is an optional list of the edges in G to be drawn.
    If provided, only the edges in edgelist will be drawn. 

    edgecolor can be a list of matplotlib color letters such as 'k' or
    'b' that lists the color of each edge; the list must be ordered in
    the same way as the edge list. Alternatively, this list can contain
    numbers and those number are mapped to a color scale using the color
    map edge_cmap.
    
    For directed graphs, "arrows" (actually just thicker stubs) are drawn
    at the head end.  Arrows can be turned off with keyword arrows=False.

    See draw_networkx for the list of other optional parameters.

    """
    if ax is None:
        ax=matplotlib.pylab.gca()

    if edgelist is None:
        edgelist=G.edges()

    if not edgelist or len(edgelist)==0: # no edges!
        return None

    # set edge positions
    edge_pos=asarray([(pos[e[0]],pos[e[1]]) for e in edgelist])

    if not cb.iterable(width):
        lw = (width,)
    else:
        lw = width

    if not cb.is_string_like(edge_color) \
           and cb.iterable(edge_color) \
           and len(edge_color)==len(edge_pos):
        if matplotlib.numerix.alltrue([cb.is_string_like(c) 
                                       for c in edge_color]):
            # (should check ALL elements)
            # list of color letters such as ['k','r','k',...]
            edge_colors = tuple([colorConverter.to_rgba(c,alpha) 
                                 for c in edge_color])
        elif matplotlib.numerix.alltrue([not cb.is_string_like(c) 
                                         for c in edge_color]):
            # numbers (which are going to be mapped with a colormap)
            edge_colors = None
        else:
            raise ValueError('edge_color must consist of either color names or numbers')
    else:
        if len(edge_color)==1:
            edge_colors = ( colorConverter.to_rgba(edge_color, alpha), )
        else:
            raise ValueError('edge_color must be a single color or list of exactly m colors where m is the number or edges')

    edge_collection = LineCollection(edge_pos,
                                colors       = edge_colors,
                                linewidths   = lw,
                                antialiaseds = (1,),
                                linestyle    = style,     
                                transOffset = ax.transData,             
                                )
    edge_collection.set_alpha(alpha)

    # need 0.87.7 or greater for edge colormaps
    mpl_version=matplotlib.__version__
    if mpl_version.endswith('svn'):
        mpl_version=matplotlib.__version__[0:-3]
    if mpl_version.endswith('pre'):
        mpl_version=matplotlib.__version__[0:-3]
    if map(int,mpl_version.split('.'))>=[0,87,7]:
        if edge_colors is None:
            if edge_cmap is not None: assert(isinstance(edge_cmap, Colormap))
            edge_collection.set_array(asarray(edge_color))
            edge_collection.set_cmap(edge_cmap)
            if edge_vmin is not None or edge_vmax is not None:
                edge_collection.set_clim(edge_vmin, edge_vmax)
            else:
                edge_collection.autoscale()
            matplotlib.pylab.sci(edge_collection)

#    else:
#        sys.stderr.write(\
#            """matplotlib version >= 0.87.7 required for colormapped edges.
#        (version %s detected)."""%matplotlib.__version__)
#        raise UserWarning(\
#            """matplotlib version >= 0.87.7 required for colormapped edges.
#        (version %s detected)."""%matplotlib.__version__)

    arrow_collection=None

    if G.directed and arrows:

        # a directed graph hack
        # draw thick line segments at head end of edge
        # waiting for someone else to implement arrows that will work 
        arrow_colors = ( colorConverter.to_rgba('k', alpha), )
        a_pos=[]
        p=1.0-0.25 # make head segment 25 percent of edge length
        for src,dst in edge_pos:
            x1,y1=src
            x2,y2=dst
            dx=x2-x1 # x offset
            dy=y2-y1 # y offset
            d=sqrt(float(dx**2+dy**2)) # length of edge
            if d==0: # source and target at same position
                continue
            if dx==0: # vertical edge
                xa=x2
                ya=dy*p+y1
            if dy==0: # horizontal edge
                ya=y2
                xa=dx*p+x1
            else:
                theta=arctan2(dy,dx)
                xa=p*d*cos(theta)+x1
                ya=p*d*sin(theta)+y1
                
            a_pos.append(((xa,ya),(x2,y2)))

        arrow_collection = LineCollection(a_pos,
                                colors       = arrow_colors,
                                linewidths   = [4*ww for ww in lw],
                                antialiaseds = (1,),
                                transOffset = ax.transData,             
                                )
        
    # update view        
    minx = amin(ravel(edge_pos[:,:,0]))
    maxx = amax(ravel(edge_pos[:,:,0]))
    miny = amin(ravel(edge_pos[:,:,1]))
    maxy = amax(ravel(edge_pos[:,:,1]))



    w = maxx-minx
    h = maxy-miny
    padx, pady = 0.05*w, 0.05*h
    corners = (minx-padx, miny-pady), (maxx+padx, maxy+pady)
    ax.update_datalim( corners)
    ax.autoscale_view()

    edge_collection.set_zorder(1) # edges go behind nodes            
    ax.add_collection(edge_collection)
    if arrow_collection:
        arrow_collection.set_zorder(1) # edges go behind nodes            
        ax.add_collection(arrow_collection)


    return edge_collection


def draw_networkx_labels(G, pos,
                         labels=None,
                         font_size=12,
                         font_color='k',
                         font_family='sans-serif',
                         font_weight='normal',
                         alpha=1.0,
                         ax=None,
                         **kwds):
    """Draw node labels on the graph G

    pos is a dictionary keyed by vertex with a two-tuple
    of x-y positions as the value.
    See networkx.layout for functions that compute node positions.

    labels is an optional dictionary keyed by vertex with node labels
    as the values.  If provided only labels for the keys in the dictionary
    are drawn.
    
    See draw_networkx for the list of other optional parameters.

    """
    if ax is None:
        ax=matplotlib.pylab.gca()

    if labels is None:
        labels=dict(zip(G.nodes(),G.nodes()))

    text_items={}  # there is no text collection so we'll fake one        
    for (n,label) in labels.items():
        (x,y)=pos[n]
        if not cb.is_string_like(label):
            label=str(label) # this will cause "1" and 1 to be labeled the same
        t=ax.text(x, y,
                label,
                size=font_size,
                color=font_color,
                family=font_family,
                weight=font_weight,
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transData,
                )
        text_items[n]=t

    return text_items

def draw_circular(G, **kwargs):
    """Draw the graph G with a circular layout"""
    from networkx.drawing.layout import circular_layout
    draw(G,circular_layout(G),**kwargs)
    
def draw_random(G, **kwargs):
    """Draw the graph G with a random layout."""
    from networkx.drawing.layout import random_layout
    draw(G,random_layout(G),**kwargs)

def draw_spectral(G, **kwargs):
    """Draw the graph G with a spectral layout."""
    from networkx.drawing.layout import spectral_layout
    draw(G,spectral_layout(G),**kwargs)

def draw_spring(G, **kwargs):
    """Draw the graph G with a spring layout"""
    from networkx.drawing.layout import spring_layout
    draw(G,spring_layout(G),**kwargs)

def draw_shell(G, **kwargs):
    """Draw networkx graph with shell layout"""
    from networkx.drawing.layout import shell_layout
    nlist = kwargs.get('nlist', None)
    if nlist != None:        
        del(kwargs['nlist'])
    draw(G,shell_layout(G,nlist=nlist),**kwargs)

def draw_graphviz(G, prog="neato", **kwargs):
    """Draw networkx graph with graphviz layout"""
    pos=networkx.drawing.graphviz_layout(G,prog)
    draw(G,pos,**kwargs)

def draw_nx(G,pos,**kwds):
    """For backward compatibility; use draw or draw_networkx"""
    draw(G,pos,**kwds)
