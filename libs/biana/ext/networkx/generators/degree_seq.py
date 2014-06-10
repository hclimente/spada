"""
Generate graphs with a given degree sequence or expected degree sequence.

"""
#    Copyright (C) 2004-2008 by 
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    Distributed under the terms of the GNU Lesser General Public License
#    http://www.gnu.org/copyleft/lesser.html

__author__ = """Aric Hagberg (hagberg@lanl.gov)\nPieter Swart (swart@lanl.gov)\nDan Schult (dschult@colgate.edu)"""


__all__ = ['configuration_model',
           'expected_degree_graph',
           'havel_hakimi_graph',
           'degree_sequence_tree',
           'is_valid_degree_sequence',
           'create_degree_sequence',
           'double_edge_swap',
           'connected_double_edge_swap',
           'li_smax_graph',
           's_metric']


import random
import networkx
import networkx.utils
from networkx.generators.classic import empty_graph
import heapq
#---------------------------------------------------------------------------
#  Generating Graphs with a given degree sequence
#---------------------------------------------------------------------------

def configuration_model(deg_sequence,seed=None):
    """Return a random pseudograph with the given degree sequence.

      - `deg_sequence`: degree sequence, a list of integers with each entry
                        corresponding to the degree of a node (need not be
                        sorted). A non-graphical degree sequence (i.e. one
                        not realizable by some simple graph) will raise an
                        Exception.
      - `seed`: seed for random number generator (default=None)


    >>> from networkx.utils import powerlaw_sequence
    >>> z=nx.create_degree_sequence(100,powerlaw_sequence)
    >>> G=nx.configuration_model(z)

    The pseudograph G is a networkx.MultiGraph that allows multiple (parallel) 
    edges between nodes and self-loops (edges from a node to itself).

    To remove parallel edges:

    >>> G=nx.Graph(G)

    Steps:

     - Check if deg_sequence is a valid degree sequence.
     - Create N nodes with stubs for attaching edges
     - Randomly select two available stubs and connect them with an edge.

    As described by Newman [newman-2003-structure].
    
    Nodes are labeled 1,.., len(deg_sequence),
    corresponding to their position in deg_sequence.

    This process can lead to duplicate edges and loops, and therefore
    returns a pseudograph type.  You can remove the self-loops and
    parallel edges (see above) with the likely result of
    not getting the exat degree sequence specified.
    This "finite-size effect" decreases as the size of the graph increases.

    References:
    
    [newman-2003-structure]  M.E.J. Newman, "The structure and function
    of complex networks", SIAM REVIEW 45-2, pp 167-256, 2003.
        
    """
    if not sum(deg_sequence)%2 ==0:
        raise networkx.NetworkXError, 'Invalid degree sequence'

    if not seed is None:
        random.seed(seed)

    # start with empty N-node graph
    N=len(deg_sequence)
#    G=networkx.empty_graph(N,create_using=networkx.Graph()) # no multiedges or selfloops

    # allow multiedges and selfloops
    G=networkx.empty_graph(N,create_using=networkx.MultiGraph())

    if N==0 or max(deg_sequence)==0: # done if no edges
        return G 

    # build stublist, a list of available degree-repeated stubs
    # e.g. for deg_sequence=[3,2,1,1,1]
    # initially, stublist=[1,1,1,2,2,3,4,5]
    # i.e., node 1 has degree=3 and is repeated 3 times, etc.
    stublist=[]
    for n in G:
        for i in range(deg_sequence[n-1]):
            stublist.append(n)

    # shuffle stublist and assign pairs by removing 2 elements at a time      
    random.shuffle(stublist)
    while stublist:
        n1 = stublist.pop()
        n2 = stublist.pop()
        G.add_edge(n1,n2)

    G.name="configuration_model %d nodes %d edges"%(G.order(),G.size())
    return G


def expected_degree_graph(w, seed=None):
    """Return a random graph G(w) with expected degrees given by w.

    :Parameters:
       - `w`: a list of expected degrees
       - `seed`: seed for random number generator (default=None)

    >>> z=[10 for i in range(100)]
    >>> G=nx.expected_degree_graph(z)


    Reference::

      @Article{connected-components-2002,
        author =	{Fan Chung and L. Lu},
        title = 	{Connected components in random graphs
        with given expected degree sequences},
        journal = 	{Ann. Combinatorics},
        year = 		{2002},
        volume = 	{6},
        pages = 	{125-145},
        }


	"""

    n = len(w)
    # allow self loops
    G=networkx.empty_graph(n,create_using=networkx.Graph())
    G.name="random_expected_degree_graph"

    if n==0 or max(w)==0: # done if no edges
        return G 

    d = sum(w)
    rho = 1.0 / float(d) # Vol(G)
    for i in xrange(n):
        if (w[i])**2 > d:
            raise networkx.NetworkXError,\
                  "NetworkXError w[i]**2 must be <= sum(w)\
                  for all i, w[%d] = %f, sum(w) = %f" % (i,w[i],d)

    if seed is not None:
        random.seed(seed)
	
    for u in xrange(n):
        for v in xrange(u,n):
            if random.random() < w[u]*w[v]*rho:
                G.add_edge(u,v)
    return G 


def havel_hakimi_graph(deg_sequence):
    """Return a simple graph with given degree sequence, constructed using the
    Havel-Hakimi algorithm.

      - `deg_sequence`: degree sequence, a list of integers with each entry
         corresponding to the degree of a node (need not be sorted).
         A non-graphical degree sequence (not sorted).
         A non-graphical degree sequence (i.e. one
         not realizable by some simple graph) raises an Exception.    

    The Havel-Hakimi algorithm constructs a simple graph by
    successively connecting the node of highest degree to other nodes
    of highest degree, resorting remaining nodes by degree, and
    repeating the process. The resulting graph has a high
    degree-associativity.  Nodes are labeled 1,.., len(deg_sequence),
    corresponding to their position in deg_sequence.

    See Theorem 1.4 in [chartrand-graphs-1996].
    This algorithm is also used in the function is_valid_degree_sequence.

    References:

    [chartrand-graphs-1996] G. Chartrand and L. Lesniak, "Graphs and Digraphs",
                            Chapman and Hall/CRC, 1996.

    """
    if not is_valid_degree_sequence(deg_sequence):
        raise networkx.NetworkXError, 'Invalid degree sequence'


    N=len(deg_sequence)
    G=networkx.empty_graph(N) # always return a simple graph

    if N==0 or max(deg_sequence)==0: # done if no edges
        return G 
 
    # form list of [stubs,name] for each node.
    stublist=[ [deg_sequence[n],n] for n in G]
    #  Now connect the stubs
    while stublist:
        stublist.sort()
        if stublist[0][0]<0: # took too many off some vertex
            return False     # should not happen if deg_seq is valid

        (freestubs,source) = stublist.pop() # the node with the most stubs
        if freestubs==0: break          # the rest must also be 0 --Done!
        if freestubs > len(stublist):  # Trouble--can't make that many edges
            return False               # should not happen if deg_seq is valid

        # attach edges to biggest nodes
        for stubtarget in stublist[-freestubs:]:
            G.add_edge(source, stubtarget[1])  
            stubtarget[0] -= 1  # updating stublist on the fly
    
    G.name="havel_hakimi_graph %d nodes %d edges"%(G.order(),G.size())
    return G

def degree_sequence_tree(deg_sequence):
    """
    Make a tree for the given degree sequence.

    A tree has #nodes-#edges=1 so
    the degree sequence must have
    len(deg_sequence)-sum(deg_sequence)/2=1
    """

    if not len(deg_sequence)-sum(deg_sequence)/2.0 == 1.0:
        raise networkx.NetworkXError,"Degree sequence invalid"

    G=empty_graph(0)
    # single node tree
    if len(deg_sequence)==1:
        return G
    deg=[s for s in deg_sequence if s>1] # all degrees greater than 1
    deg.sort(reverse=True)

    # make path graph as backbone
    n=len(deg)+2
    G=networkx.path_graph(n)
    last=n

    # add the leaves
    for source in range(1,n-1):
        nedges=deg.pop()-2
        for target in range(last,last+nedges):
            G.add_edge(source, target)
        last+=nedges
        
    # in case we added one too many 
    if len(G.degree())>len(deg_sequence): 
        G.delete_node(0)
    return G
        

def is_valid_degree_sequence(deg_sequence):
    """Return True if deg_sequence is a valid sequence of integer degrees
    equal to the degree sequence of some simple graph.
       
      - `deg_sequence`: degree sequence, a list of integers with each entry
         corresponding to the degree of a node (need not be sorted).
         A non-graphical degree sequence (i.e. one not realizable by some
         simple graph) will raise an exception.
                        
    See Theorem 1.4 in [chartrand-graphs-1996]. This algorithm is also used
    in havel_hakimi_graph()

    References:

    [chartrand-graphs-1996] G. Chartrand and L. Lesniak, "Graphs and Digraphs",
                            Chapman and Hall/CRC, 1996.

    """
    # some simple tests 
    if deg_sequence==[]:
        return True # empty sequence = empty graph 
    if not networkx.utils.is_list_of_ints(deg_sequence):
        return False   # list of ints
    if min(deg_sequence)<0:
        return False      # each int not negative
    if sum(deg_sequence)%2:
        return False      # must be even
    
    # successively reduce degree sequence by removing node of maximum degree
    # as in Havel-Hakimi algorithm
        
    s=deg_sequence[:]  # copy to s
    while s:      
        s.sort()    # sort in non-increasing order
        if s[0]<0: 
            return False  # check if removed too many from some node

        d=s.pop()             # pop largest degree 
        if d==0: return True  # done! rest must be zero due to ordering

        # degree must be <= number of available nodes
        if d>len(s):   return False

        # remove edges to nodes of next higher degrees
        s.reverse()  # to make it easy to get at higher degree nodes.
        for i in range(d):
            s[i]-=1

    # should never get here b/c either d==0, d>len(s) or d<0 before s=[]
    return False


def create_degree_sequence(n, sfunction=None, max_tries=50, **kwds):
    """ Attempt to create a valid degree sequence of length n using
    specified function sfunction(n,**kwds).

      - `n`: length of degree sequence = number of nodes
      - `sfunction`: a function, called as "sfunction(n,**kwds)",
         that returns a list of n real or integer values.
      - `max_tries`: max number of attempts at creating valid degree
         sequence.

    Repeatedly create a degree sequence by calling sfunction(n,**kwds)
    until achieving a valid degree sequence. If unsuccessful after
    max_tries attempts, raise an exception.
    
    For examples of sfunctions that return sequences of random numbers,
    see networkx.Utils.

    >>> from networkx.utils import uniform_sequence
    >>> seq=nx.create_degree_sequence(10,uniform_sequence)

    """
    tries=0
    max_deg=n
    while tries < max_tries:
        trialseq=sfunction(n,**kwds)
        # round to integer values in the range [0,max_deg]
        seq=[min(max_deg, max( int(round(s)),0 )) for s in trialseq]
        # if graphical return, else throw away and try again
        if is_valid_degree_sequence(seq):
            return seq
        tries+=1
    raise networkx.NetworkXError, \
          "Exceeded max (%d) attempts at a valid sequence."%max_tries

def double_edge_swap(G, nswap=1):
    """Attempt nswap double-edge swaps on the graph G.

    Return count of successful swaps.
    The graph G is modified in place.
    A double-edge swap removes two randomly choseen edges u-v and x-y
    and creates the new edges u-x and v-y::

     u--v            u  v
            becomes  |  |
     x--y            x  y


    If either the edge u-x or v-y already exist no swap is performed so
    the actual count of swapped edges is always <= nswap
    
    Does not enforce any connectivity constraints.
    """
    # this algorithm and connected_double_edge_swap avoid choosing
    # uniformly at random from a generated edge list by instead
    # choosing nonuniformly from the set nodes (probability weighted by degree)
    n=0
    swapcount=0
    deg=G.degree(with_labels=True)
    dk=deg.keys() # key labels 
    cdf=networkx.utils.cumulative_distribution(deg.values())  # cdf of degree
    if len(cdf)<4:
        raise networkx.NetworkXError, \
          "Graph has less than four nodes."
    while n < nswap:
#        if random.random() < 0.5: continue # trick to avoid periodicities?
        # pick two randon edges without creating edge list
        # choose source node indices from discrete distribution
        (ui,xi)=networkx.utils.discrete_sequence(2,cdistribution=cdf) 
        if ui==xi: continue # same source, skip
        u=dk[ui] # convert index to label
        x=dk[xi] 
        v=random.choice(G.neighbors(u)) # choose target uniformly from nbrs
        y=random.choice(G.neighbors(x)) # Note: dan't use G[u] because choice can't use dict 
        if v==y: continue # same target, skip
        if (x not in G[u]) and (y not in G[v]):
            G.add_edge(u,x)
            G.add_edge(v,y)
            G.delete_edge(u,v)
            G.delete_edge(x,y)
            swapcount+=1
        n+=1
    return swapcount

def connected_double_edge_swap(G, nswap=1):
    """Attempt nswap double-edge swaps on the graph G.

    Returns count of successful swaps.  Enforces connectivity.
    The graph G is modified in place.

    A double-edge swap removes two randomly choseen edges u-v and x-y
    and creates the new edges u-x and v-y::

     u--v            u  v
            becomes  |  |
     x--y            x  y


    If either the edge u-x or v-y already exist no swap is performed so
    the actual count of swapped edges is always <= nswap

    The initial graph G must be connected and the resulting graph is connected.

    Reference::

     @misc{gkantsidis-03-markov,
      author = "C. Gkantsidis and M. Mihail and E. Zegura",
      title = "The Markov chain simulation method for generating connected
               power law random graphs",
      year = "2003",
      url = "http://citeseer.ist.psu.edu/gkantsidis03markov.html"
     }


    """
    import math
    if not networkx.is_connected(G):
       raise networkx.NetworkXException("Graph not connected")

    n=0
    swapcount=0
    deg=G.degree(with_labels=True)
    dk=deg.keys() # key labels 
    ideg=G.degree()
    cdf=networkx.utils.cumulative_distribution(G.degree()) 
    if len(cdf)<4:
        raise networkx.NetworkXError, \
          "Graph has less than four nodes."
    window=1
    while n < nswap:
        wcount=0
        swapped=[]
        while wcount < window and n < nswap:
            # pick two randon edges without creating edge list
            # chose source nodes from discrete distribution
            (ui,xi)=networkx.utils.discrete_sequence(2,cdistribution=cdf) 
            if ui==xi: continue # same source, skip
            u=dk[ui] # convert index to label
            x=dk[xi] 
            v=random.choice(G.neighbors(u)) # choose target uniformly from nbrs
            y=random.choice(G.neighbors(x)) # Note: dan't use G[u] because choice can't use dict 
            if v==y: continue # same target, skip
            if (not G.has_edge(u,x)) and (not G.has_edge(v,y)):
                G.delete_edge(u,v)
                G.delete_edge(x,y)
                G.add_edge(u,x)
                G.add_edge(v,y)
                swapped.append((u,v,x,y))
                swapcount+=1
            n+=1
            wcount+=1
        if networkx.is_connected(G): # increase window
            window+=1
        else: # undo changes from previous window, decrease window
            while swapped:
                (u,v,x,y)=swapped.pop()
                G.add_edge(u,v)
                G.add_edge(x,y)
                G.delete_edge(u,x)
                G.delete_edge(v,y)
                swapcount-=1
            window = int(math.ceil(float(window)/2))
        assert G.degree() == ideg
    return swapcount



def li_smax_graph(degree_seq):
    """Generates a graph based with a given degree sequence and maximizing
    the s-metric.  Experimental implementation.

    Maximum s-metrix  means that high degree nodes are connected to high
    degree nodes. 
        
    - `degree_seq`: degree sequence, a list of integers with each entry
       corresponding to the degree of a node.
       A non-graphical degree sequence raises an Exception.    

    Reference::      
    
      @unpublished{li-2005,
       author = {Lun Li and David Alderson and Reiko Tanaka
                and John C. Doyle and Walter Willinger},
       title = {Towards a Theory of Scale-Free Graphs:
               Definition, Properties, and  Implications (Extended Version)},
       url = {http://arxiv.org/abs/cond-mat/0501169},
       year = {2005}
      }

    The algorithm::

     STEP 0 - Initialization
     A = {0}
     B = {1, 2, 3, ..., n}
     O = {(i; j), ..., (k, l),...} where i < j, i <= k < l and 
             d_i * d_j >= d_k *d_l 
     wA = d_1
     dB = sum(degrees)

     STEP 1 - Link selection
     (a) If |O| = 0 TERMINATE. Return graph A.
     (b) Select element(s) (i, j) in O having the largest d_i * d_j , if for 
             any i or j either w_i = 0 or w_j = 0 delete (i, j) from O
     (c) If there are no elements selected go to (a).
     (d) Select the link (i, j) having the largest value w_i (where for each 
             (i, j) w_i is the smaller of w_i and w_j ), and proceed to STEP 2.

     STEP 2 - Link addition
     Type 1: i in A and j in B. 
             Add j to the graph A and remove it from the set B add a link 
             (i, j) to the graph A. Update variables:
             wA = wA + d_j -2 and dB = dB - d_j
             Decrement w_i and w_j with one. Delete (i, j) from O
     Type 2: i and j in A.
         Check Tree Condition: If dB = 2 * |B| - wA. 
             Delete (i, j) from O, continue to STEP 3
         Check Disconnected Cluster Condition: If wA = 2. 
             Delete (i, j) from O, continue to STEP 3
         Add the link (i, j) to the graph A 
         Decrement w_i and w_j with one, and wA = wA -2
     STEP 3
         Go to STEP 1


    The article states that the algorithm will result in a maximal s-metric. 
    This implementation can not guarantee such maximality. I may have 
    misunderstood the algorithm, but I can not see how it can be anything 
    but a heuristic. Please contact me at sundsdal@gmail.com if you can 
    provide python code that can guarantee maximality.
    Several optimizations are included in this code and it may be hard to read.
    Commented code to come.
    """
    
    if not is_valid_degree_sequence(degree_seq):
        raise networkx.NetworkXError, 'Invalid degree sequence'
    degree_seq.sort() # make sure it's sorted
    degree_seq.reverse()
    degrees_left = degree_seq[:]
    A_graph = networkx.Graph()
    A_graph.add_node(0)
    a_list = [False]*len(degree_seq)
    b_set = set(range(1,len(degree_seq)))
    a_open = set([0])
    O = []
    for j in b_set:
        heapq.heappush(O, (-degree_seq[0]*degree_seq[j], (0,j)))
    wa = degrees_left[0] #stubs in a_graph
    db = sum(degree_seq) - degree_seq[0] #stubs in b-graph
    a_list[0] = True #node 0 is now in a_Graph
    bsize = len(degree_seq) -1 #size of b_graph
    selected = []
    weight = 0
    while O or selected:
        if len(selected) <1 :
            firstrun = True
            while O:
                (newweight, (i,j)) = heapq.heappop(O)
                if degrees_left[i] < 1 or degrees_left[j] < 1 :
                    continue
                if firstrun:
                    firstrun = False
                    weight = newweight
                if not newweight == weight:
                    break
                heapq.heappush(selected, [-degrees_left[i], \
                                    -degrees_left[j], (i,j)])
            if not weight == newweight:
                heapq.heappush(O,(newweight, (i,j)))
            weight *= -1
        if len(selected) < 1:
            break
        
        [w1, w2, (i,j)] = heapq.heappop(selected)
        if degrees_left[i] < 1 or degrees_left[j] < 1 :
            continue
        if a_list[i] and j in b_set:
            #TYPE1
            a_list[j] = True
            b_set.remove(j)
            A_graph.add_node(j)
            A_graph.add_edge(i, j)
            degrees_left[i] -= 1
            degrees_left[j] -= 1
            wa += degree_seq[j] - 2
            db -= degree_seq[j]
            bsize -= 1
            newweight = weight
            if not degrees_left[j] == 0:
                a_open.add(j)
                for k in b_set:
                    if A_graph.has_edge(j, k): continue
                    w = degree_seq[j]*degree_seq[k]
                    if w > newweight:
                        newweight = w
                    if weight == w and not newweight > weight:
                        heapq.heappush(selected, [-degrees_left[j], \
                                            -degrees_left[k], (j,k)])
                    else:
                        heapq.heappush(O, (-w, (j,k)))
                if not weight == newweight:
                    while selected:
                        [w1,w2,(i,j)] = heapq.heappop(selected)
                        if degrees_left[i]*degrees_left[j] > 0:
                            heapq.heappush(O, [-degree_seq[i]*degree_seq[j],(i,j)])
            if degrees_left[i] == 0:
                a_open.discard(i)
                    
        else:
            #TYPE2
            if db == (2*bsize - wa):
                #tree condition
                #print "removing because tree condition    "
                continue
            elif db < 2*bsize -wa:
                raise networkx.NetworkXError, "THIS SHOULD NOT HAPPEN! - not graphable"
                continue
            elif wa == 2 and bsize > 0:
                #print "removing because disconnected  cluster"
                #disconnected cluster condition
                continue
            elif wa == db - (bsize)*(bsize-1):
                #print "MYOWN removing because disconnected  cluster"
                continue
            A_graph.add_edge(i, j)
            degrees_left[i] -= 1
            degrees_left[j] -= 1
            if degrees_left[i] < 1:
                a_open.discard(i)
            if degrees_left[j] < 1:
                a_open.discard(j)
            wa -=  2
            if not degrees_left[i] < 0 and not degrees_left[j] < 0:
                selected2 = (selected)
                selected = []
                while selected2:
                    [w1,w1, (i,j)] = heapq.heappop(selected2)
                    if degrees_left[i]*degrees_left[j] > 0:
                        heapq.heappush(selected, [-degrees_left[i], \
                                        -degrees_left[j], (i,j)])
    return A_graph 

def connected_smax_graph(degree_seq):
    """
    Not implemented.
    """
    # incomplete implementation
    
    if not is_valid_degree_sequence(degree_seq):
        raise networkx.NetworkXError, 'Invalid degree sequence'
   
    # build dictionary of node id and degree, sorted by degree, largest first
    degree_seq.sort() 
    degree_seq.reverse()
    ddict=dict(zip(xrange(len(degree_seq)),degree_seq))

    G=empty_graph(1) # start with single node

    return False

def s_metric(G):
    """
    Return the "s-Metric" of graph G:
    the sum of the product deg(u)*deg(v) for every edge u-v in G
    
    Reference::

     @unpublished{li-2005,
      author = {Lun Li and David Alderson and
               John C. Doyle and Walter Willinger},
      title = {Towards a Theory of Scale-Free Graphs:
               Definition, Properties, and  Implications (Extended Version)},
      url = {http://arxiv.org/abs/cond-mat/0501169},
      year = {2005}
      }

    """
    # this function doesn't belong in this module
    return sum([G.degree(u)*G.degree(v) for (u,v) in G.edges_iter()])
   
