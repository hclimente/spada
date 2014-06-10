"""
    BIANA: Biologic Interactions and Network Analysis
    Copyright (C) 2009  Javier Garcia-Garcia, Emre Guney, Baldo Oliva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#########################################################################
# Graph Utility Methods
# Methods to
# - create a network
# - filter a network based on degree
# - find paths and (connected) components of given network
# - randomize a network
# - read a network from SIF file
# - analyze network 
#
#########################################################################

#import networkx
import biana.ext.biana_networkx as networkx
import random

#MIN_NUMBER_OF_PERTURBATION = 25
MAX_NUMBER_OF_TRIAL = 6

def create_graph_with_same_type(G):
    return create_empty_copy(G)

def create_graph():
    """
        Creates & returns a graph
    """
    return networkx.Graph()


def get_number_of_distinct_edges(G):
    edge_list = G.edges()
    edge_set = set()
    for id1, id2, data in edge_list:
        edge_set.add((id1, id2))
    return len(edge_set)

def get_shortest_path_between(G, source_id, target_id):
    return networkx.shortest_path(G, source_id, target_id)
        
def get_all_paths_from(G, source_id): 
    """
        get all paths from source node to all possible nodes in a dictionary
    """
    return networkx.dijkstra_path(G, source_id)

def get_path_network(G, listNodes, path_length_cutoff=10000):
    """
    Returns a subgraph containing only nodes in listNodes

    Nodes are connected if exist a path between them, and weight of the edge consists in the length of shortest path

    If the shortest between two nodes includes another node in the list, this edge is not added
    """
    # First, check all nodes in listNodes are in network
    new_graph = networkx.MultiGraph()

    for x in xrange(len(listNodes)):
        for y in xrange(x):
            sp = networkx.shortest_path(G, listNodes[x], listNodes[y])
            if sp:
                if (len(sp)-1)<=path_length_cutoff:
                    if len(set(listNodes).intersection(set(sp[1:-1])))==0:
                        new_graph.add_edge(listNodes[x],listNodes[y],len(sp)-1)
            
    return new_graph

def get_connected_components(G, return_as_graph_list=True):
    """
        Finds (strongly in the case of directed network) connected components of graph
        returnAsGraphList: returns list of graph objects corresponding to connected components (from larger to smaller)
        otherwise returns list of node list corresponding nodes in connected components
    """
    result_list = []

    if return_as_graph_list:
        result_list = networkx.connected_component_subgraphs(G)
    else:
        result_list = networkx.connected_components(G)

    return result_list

def randomize_graph(graph, randomization_type, allow_self_edges = True):
    """
    Creates a random network from given network as a networkx graph
    randomization_type: 
        - "random": add same number of edges randomly between nodes of original graph
        - "preserve_topology": keep edges, shuffle nodes of original graph
        - "preserve_topology_and_node_degree": keep edges, shuffle nodes of original graph with the nodes of same degree
        - "preserve_degree_distribution": remove an edge between two random nodes with degrees k, l then add to two nodes with degrees k-1 & l-1, then shuffle nodes
        - "preserve_degree_distribution_and_node_degree": remove 2 random edges between a-b and c-d where degree(a)=degree(c) and degree(b)=degree(d) then add 2 edges between a-d and b-c, then shuffle nodes with the same degree
	- "erdos_renyi": creates a graph where edges are redistributed based on erdos renyi random model. 
	- "barabasi_albert": creates a graph where edges are redistributed based on barabasi albert model (preferential attachment). 
    """

    debug = False

    n_node = graph.number_of_nodes()
    n_edge = graph.number_of_edges()

    if randomization_type == "erdos_renyi":
	#raise Exception("Work in progress")
	p = float(2 * n_edge) / (n_node*n_node - 2*n_node)
	# Chooses each of the possible [n(n-1)]/2 edges with probability p 
	new_graph = networkx.erdos_renyi_graph(n_node, p)
	mapping = dict(zip(new_graph.nodes(), graph.nodes()))
	new_graph = networkx.relabel_nodes(new_graph, mapping)
	available_edges = graph.edges()
	
	# Map graph from random model to new graph
        for edge in new_graph.edges():
	    if len(available_edges) > 0:
		edge_org = available_edges.pop()
		if debug:
		    print "From random:", (edge[0], edge[1])
		new_graph.add_edge(edge[0], edge[1], graph.get_edge(edge_org[0], edge_org[1]))
	    # If the random model added too many edges
	    else:
		if debug:
		    print "Removing:", edge
		new_graph.remove_edge(edge[0], edge[1])

	# If the random model failed to add enough edges
	nodes = new_graph.nodes()
	for edge_org in available_edges:
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
	    if debug:
		print "Adding:", (source_id, target_id)
	    new_graph.add_edge(source_id, target_id, graph.get_edge(edge_org[0], edge_org[1]))
	return new_graph

    if randomization_type == "barabasi_albert":
	#raise Exception("Work in progress")
	if n_edge >= n_node:
	    # A graph of n nodes is grown by attaching new nodes each with m edges that are preferentially attached to existing nodes with high degree
	    new_graph = networkx.barabasi_albert_graph(n_node, n_edge / n_node)
	    mapping = dict(zip(new_graph.nodes(), graph.nodes()))
	    new_graph = networkx.relabel_nodes(new_graph, mapping)
	else:
	    new_graph = networkx.create_empty_copy(graph) 

	available_edges = graph.edges() 
	degree_map = networkx.degree(new_graph, with_labels=True)
	nodes = new_graph.nodes()

	# Map graph from random model to new graph
        for edge in new_graph.edges():
	    if len(available_edges) > 0:
		edge_org = available_edges.pop()
		if debug:
		    print "From random:", (edge[0], edge[1])
		new_graph.add_edge(edge[0], edge[1], graph.get_edge(edge_org[0], edge_org[1]))
	    # If the random model added too many edges
	    else:
		nodes_to_select = [ id for id, d in degree_map.items() for j in xrange(d+1) ]
		source_id = random.choice(nodes())
		target_id = random.choice(nodes_to_select)
		if debug:
		    print "Removing:", (source_id, target_id)
		new_graph.remove_edge(source_id, target_id)
		degree_map[source_id] -= 1 
		degree_map[target_id] -= 1 

	# If the random model failed to add enough edges
	for edge_org in available_edges:
	    nodes_to_select = [ id for id, d in degree_map.items() for j in xrange(d+1) ]
            source_id = random.choice(nodes)
            target_id = random.choice(nodes_to_select)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
		source_id = random.choice(nodes)
		target_id = random.choice(nodes_to_select)
	    if debug:
		print "Adding:", (source_id, target_id)
	    new_graph.add_edge(source_id, target_id, graph.get_edge(edge_org[0], edge_org[1]))
	    degree_map[source_id] += 1 
	    degree_map[target_id] += 1 

	return new_graph

    new_graph = networkx.create_empty_copy(graph) 
    #new_graph.add_nodes_from(graph.nodes())

    if randomization_type == "random":
	nodes = new_graph.nodes()
        for edge in graph.edges():
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            new_graph.add_edge(source_id, target_id, graph.get_edge(edge[0],edge[1]))
        
    elif randomization_type=="preserve_topology": # shuffle_nodes
        nodes = graph.nodes()
        random_nodes = graph.nodes()
        random.shuffle(random_nodes)
        equivalences = dict([(nodes[i],random_nodes[i]) for i in xrange(len(nodes))])
        new_graph.add_edges_from([ (equivalences[current_edge[0]],equivalences[current_edge[1]],graph.get_edge(current_edge[0],current_edge[1])) for current_edge in graph.edges() ])

    elif randomization_type=="preserve_topology_and_node_degree": # shuffle_nodes_within_same_degree
        nodes_by_degree = dict( (degree,[]) for degree in graph.degree() )
        graph_degree = graph.degree(with_labels=True)
        [ nodes_by_degree[graph_degree[node]].append(node) for node in graph_degree ]
        equivalences = {}
        for current_degree in nodes_by_degree.keys():
            nodes = nodes_by_degree[current_degree]
            random_nodes = list(nodes)
            random.shuffle(random_nodes)
            equivalences.update(dict([(nodes[i],random_nodes[i]) for i in xrange(len(nodes))]))
        new_graph.add_edges_from([ (equivalences[current_edge[0]],equivalences[current_edge[1]], graph.get_edge(current_edge[0],current_edge[1])) for current_edge in graph.edges() ])
        
    elif randomization_type=="preserve_degree_distribution":
        ## add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge(current_node1, current_node2))
        max_degree = sorted(graph.degree())[-1]
        #nodes_by_degree = dict( (degree,{}) for degree in graph.degree() )
        nodes_by_degree = dict( (degree,{}) for degree in xrange(max_degree+1) )
        graph_degree = graph.degree(with_labels=True)
        [ nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree ]
        #print new_graph.nodes(), new_graph.edges()
        #print nodes_by_degree
        #if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(n_edge/2, n_edge)
        #else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)
        n_perturbation = random.randint(n_edge/2, n_edge)
        for i in xrange(n_perturbation):
            n_trial = 0
            while True:
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
		    if debug:
			print "Warning: Max number of trials exceeded in perturbation ", i
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                while source_degree < 1:  
                    source_id = random.choice(new_graph.nodes())
                    source_degree = new_graph.degree(source_id)
                target_id = random.choice(new_graph.neighbors(source_id))
                target_degree = new_graph.degree(target_id)
                del nodes_by_degree[source_degree][source_id] 
                nodes_by_degree[source_degree-1].setdefault(source_id)
                del nodes_by_degree[target_degree][target_id] 
                nodes_by_degree[target_degree-1].setdefault(target_id)
                ## not very important to check for cases where new_source = source (v.v. for targets) 
                new_source_id = random.choice(nodes_by_degree[source_degree-1].keys())
                new_target_id = random.choice(nodes_by_degree[target_degree-1].keys())
		if debug:
		    print source_id, target_id, " / ", new_source_id, new_target_id
                ## check if going to add an existing edge or self edge
                if new_graph.has_edge(new_source_id, new_target_id) or new_source_id == new_target_id:
                    del nodes_by_degree[source_degree-1][source_id] 
                    nodes_by_degree[source_degree].setdefault(source_id)
                    del nodes_by_degree[target_degree-1][target_id] 
                    nodes_by_degree[target_degree].setdefault(target_id)
                    continue
		if debug:
		    print "rm %d %d" % (source_id, target_id)
                edge_data = new_graph.get_edge(source_id, target_id)
                new_graph.delete_edge(source_id, target_id)
		if debug:
		    print "add %d %d" % (new_source_id, new_target_id)
                new_graph.add_edge(new_source_id, new_target_id, edge_data)
                del nodes_by_degree[source_degree-1][new_source_id] 
                nodes_by_degree[source_degree].setdefault(new_source_id)
                del nodes_by_degree[target_degree-1][new_target_id] 
                nodes_by_degree[target_degree].setdefault(new_target_id)
                break
        #self.randomize_graph(new_graph, "preserve_topology")

    elif randomization_type=="preserve_degree_distribution_and_node_degree":
        ## add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge(current_node1, current_node2))
        nodes_by_degree = dict( (degree,{}) for degree in graph.degree() )
        graph_degree = graph.degree(with_labels=True)
        [ nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree ]
        
        #if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(1, n_edge)
        #else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)
        n_perturbation = random.randint(n_edge/2, n_edge)
        for i in xrange(n_perturbation):
            source_id = random.choice(new_graph.nodes())
            source_degree = new_graph.degree(source_id)
            ## find a node for which another node with the same degree exists
            #available_neighbors = []
            n_trial = 0
            while True: #(len(nodes_by_degree[source_degree]) < 2 or len(available_neighbors) < 1):
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
		    if debug:
			print "Warning: Max number of trials exceeded in perturbation ", i
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                if len(nodes_by_degree[source_degree]) < 2:
                    continue
                available_neighbors = []
                ## find a neighbor for which another node with the same degree exists
                for neighbor_id in new_graph.neighbors_iter(source_id):
                    if source_degree == new_graph.degree(neighbor_id): 
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 2:
                            available_neighbors.append(neighbor_id)
                    else:
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 1:
                            available_neighbors.append(neighbor_id)
                if len(available_neighbors) < 1:
                    continue
                target_id = random.choice(available_neighbors)
                target_degree = new_graph.degree(target_id)
                ## select a new source node with different id
		n_trial2 = 0
		inner_break = False
                while True:
		    n_trial2 += 1
		    if n_trial2 > MAX_NUMBER_OF_TRIAL:
			if debug:
			    print "Warning: Max number of trials exceeded in perturbation ", i
			inner_break = True
			break
                    new_source_id = random.choice(nodes_by_degree[source_degree].keys())
                    while new_source_id == source_id:
                        new_source_id = random.choice(nodes_by_degree[source_degree].keys())
                    new_available_neighbors = []
                    ## find a neighbor as new target node for which id is different from target and has an id equivalent to target
                    for neighbor_id in new_graph.neighbors_iter(new_source_id):
                        if target_degree == new_graph.degree(neighbor_id): 
                            new_available_neighbors.append(neighbor_id)
                    if len(new_available_neighbors) < 1:
                        continue
                    new_target_id = random.choice(new_available_neighbors)
                    if len(new_available_neighbors) > 1:
                        while new_target_id == target_id:
                            new_target_id = random.choice(new_available_neighbors)
                            #print new_available_neighbors, new_target_id
                    else:
                        new_target_id = new_available_neighbors[0]
                    break
		if inner_break:
		    break
		if debug:
		    print source_id, target_id, " / ", new_source_id, new_target_id
                if source_id == new_target_id or new_source_id == target_id:
                    continue
                if new_graph.has_edge(source_id, new_target_id) or new_graph.has_edge(new_source_id, target_id):
                    continue
		if debug:
		    print "rm %d %d" % (source_id, target_id)
		    print "rm %d %d" % (new_source_id, new_target_id)
                edge_data_1 = new_graph.get_edge(source_id, target_id)
                edge_data_2 = new_graph.get_edge(new_source_id, new_target_id)
                new_graph.delete_edge(source_id, target_id)
                new_graph.delete_edge(new_source_id, new_target_id)
		if debug:
		    print "add %d %d" % (source_id, new_target_id)
		    print "add %d %d" % (new_source_id, target_id)
                new_graph.add_edge(source_id, new_target_id, edge_data_1)
                new_graph.add_edge(new_source_id, target_id, edge_data_2)

    else:
        raise Exception("Unknown randomization type %s" % randomization_type)

    return new_graph

def create_network_from_sif_file(network_file_in_sif, use_edge_data = False, delim = None):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type = use_edge_data, delim = delim)
    g=networkx.Graph()
    if use_edge_data:
	for e,w in dictEdge.iteritems():
	    u,v = e
	    g.add_edge(u,v,w)
    else:
	g.add_edges_from(setEdge)
    return g

def output_network_in_sif(g, output_file_name, delim = " "):
    f = open(output_file_name, 'w')
    for u,v in g.edges_iter():
	f.write("%s%s%s%s%s\n" % (u, delim, g.get_edge(u,v), delim, v) )
    f.close()
    return

def analyze_network(g, calculate_radius =  False):
    # print networkx.info(g) # currently buggy but fixed in the next version of networkx
    print "(V,E): ", g.number_of_nodes(), g.number_of_edges()
    print "V/E: ", g.number_of_nodes() / float(g.number_of_edges())
    degrees = g.degree() 
    degrees.sort()
    print "Average degree: ", sum(degrees)/float(len(degrees))
    print "Most connected 20 nodes: ", degrees[-20:]
    connected_components = networkx.connected_components(g)
    print "Connected component sizes: ", map(len, connected_components)
    # Radius calculation is very time consuming
    if calculate_radius:
	print "Radius (of the largest connected component): ", get_network_radius(g.subgraph(connected_components[0]))
    #print get_network_degree_histogram(g)
    return

def filter_network(g, degree_threshold=None, largest_connected_component=True):
    print "V,E:", g.number_of_nodes(), g.number_of_edges()
    degrees = g.degree(with_labels=True)
    subgraph_nodes = []
    for id, d in degrees.iteritems():
	if degree_threshold is None or d <= degree_threshold:
	    subgraph_nodes.append(id)
    g_filtered = g.subgraph(subgraph_nodes)
    if largest_connected_component:
	component_nodes = networkx.connected_components(g_filtered)[0]
	g_filtered = g_filtered.subgraph(component_nodes )
    print "V,E filtered:", g_filtered.number_of_nodes(), g_filtered.number_of_edges()
    return g_filtered

def get_edge_values_from_sif_attribute_file(file_name, store_edge_type=False, delim=None):
    """
	store_edge_type: if True returns dict in [(u, v, "pp")] = val format, if False returns dict in [(u,v)] = val format
    """
    edge_to_values = {}
    f=open(file_name)
    line = f.readline() # read attribute name
    line = f.readline()
    while line:
	if delim is None:
	    words = line[:-1].split()
	else:
	    words = line[:-1].split(delim)
	if len(words) != 5 or words[3] != "=":
	    print "format error", line
	    continue
	if store_edge_type:
	    edge_to_values.setdefault((words[0], words[2], words[1]), set()).add(words[4])
	else:
	    edge_to_values.setdefault((words[0], words[2]), set()).add(words[4])
	line = f.readline()
    f.close()
    return edge_to_values

def get_nodes_and_edges_from_sif_file(file_name, store_edge_type = False, delim=None):
    """
	Parse sif file into node and edge sets and dictionaries
	returns setNode, setEdge, dictNode, dictEdge
	store_edge_type: if True, dictEdge[(u,v)] = edge_value
	delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    f=open(file_name)
    for line in f:
	if delim is None:
	    words = line[:-1].split()
	else:
	    words = line[:-1].split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            score = float(words[1])
            dictNode[id1] = score
        elif len(words) == 3: 
            id2 = words[2]
            setNode.add(id2)
            if store_edge_type:
                dictEdge[(id1, id2)] = words[1]
		#setEdge.add((id1, id2, words[1]))
	    else:
		setEdge.add((id1, id2))
    f.close()
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    return setNode, setEdge, dictNode, dictEdge

def get_jaccard_index_map(g):
    edge_to_jaccard = {}
    for u,v in g.edges_iter():
	u_neighbors = set(g.neighbors(u))
	v_neighbors = set(g.neighbors(v))
	edge_to_jaccard[(u,v)] = float(len(u_neighbors & v_neighbors)) / len(u_neighbors | v_neighbors)
    return edge_to_jaccard 

def get_clustering_coefficient_map(g):
    return networkx.clustering(g, with_labels=True)

def get_network_radius(g):
    return networkx.radius(g)

def get_network_degree_histogram(g):
    return networkx.degree_histogram(g)

#def create_arff_file_with_network_metrics(file_node, file_network, file_arff):
def create_arff_file_with_network_metrics(g, node_to_score, seeds, arff_file_name):
    delim = ","
    header = "@RELATION aneurysm\n@ATTRIBUTE id STRING\n@ATTRIBUTE score NUMERIC\n" + \
	        "@ATTRIBUTE degree INTEGER\n@ATTRIBUTE linker_degree INTEGER\n" + \
		"@ATTRIBUTE ld_ratio NUMERIC\n@ATTRIBUTE clustering_coefficient NUMERIC\n" + \
		"@ATTRIBUTE betweenness_centrality NUMERIC\n" + \
	        "@ATTRIBUTE degree2 INTEGER\n@ATTRIBUTE linker_degree2 INTEGER\n" + \
		"@ATTRIBUTE ld_ratio2 NUMERIC\n" + \
		"@ATTRIBUTE class {involved,not-involved}\n@DATA\n" 

    #seeds, setDummy, node_to_score, dictDummy = get_nodes_and_edges_from_sif_file(file_node, False)
    #setNode, setEdge, dictDummy, dictDummy = get_nodes_and_edges_from_sif_file(file_network, False)
    #g=networkx.Graph()
    #g.add_edges_from(setEdge)
    seeds = set(seeds)

    print "Calculating betweenness centrality.."
    mapB = networkx.betweenness_centrality(g) 
    #mapB = dict(zip(g.nodes(), range(len(g.nodes())))) 

    print "Calculating clustering coefficients.."
    mapC = networkx.clustering(g, with_labels=True) 

    #print mapB["IL10"], mapC["IL10"]
    
    #print "connected component sizes: ", map(len, networkx.connected_components(g))
    #cliques = networkx.find_cliques(g) # high computational cost
    
    f = open(arff_file_name, 'w')
    f.write(header)
    # for v in setNodes:
    for v in g.nodes_iter():
	s="?"
	d=0
	ld=0
	cc=0.0
	bc=0.0
	d2=0
	ld2=0
	n1 = 0.0
	n2 = 0.0
	c="?"
	if g.has_node(v):
	    setNeighbor = set(g.neighbors(v))
	    setNeighbor2 = set(g.neighbors(v))
	    d = len(setNeighbor)
	    ld = len(setNeighbor & seeds)
	    cc = mapC[v]
	    bc = mapB[v]
	    n1 = float(ld)/d
	    for u in setNeighbor:
		setNeighbor2 |= set(g.neighbors(u))
	    setNeighbor2.remove(v)
	    d2 = len(setNeighbor2)
	    ld2 = len(setNeighbor2 & seeds)
	    if d2 == 0:
		n2 = 0.0
	    else:
		n2 = float(ld2)/d2
	if v in seeds:
	    s=node_to_score[v]
	    c="involved"
	else:
	    s="?"
	    c="not-involved"
	#f.write("%s%s%s%s%i%s%i%s%.3f%s%.3f%s%.3f%s%i%s%i%s%.3f%s%s\n" % (v, delim, str(s), delim, d, delim, ld, delim, n1, delim, cc, delim, flaot(bc), delim, d2, delim, ld2, delim, n2, delim, c))
	# id score degree linker_degree ld_ratio clustering_coeff betweenness_cent d2 ld2 ld_ratio2 class
	# v s d ld n1 cc bc d2 ld2 n2 c
	f.write( ("%s" % delim).join( map(str, [v, s, d, ld, n1, cc, bc, d2, ld2, n2, c]) ) + "\n" )
    f.close()
    return



if __name__ == "__main__":

    test_network = networkx.Graph()
    test_network.add_nodes_from([1,2,3,4,5,6,7,8,9])
    test_network.add_edges_from([(1,2),(1,3),(1,5),(2,6),(3,4),(5,6),(4,7),(7,8)]) # (1,4)
    test_network2.add_edges_from([(1,2),(1,3),(2,3),(2,4),(4,5)])

    print "original network:"
    print test_network.edges()

    print "preserve topology:"
    random_network = randomize_graph(graph=test_network, randomization_type="preserve_topology")
    print random_network.edges()

    print "preserve topology and node degrees:"
    random_network = randomize_graph(graph=test_network, randomization_type="preserve_topology_and_node_degree")
    print random_network.edges()

    print "preserve degree distribution:"
    randomn = randomize_graph(graph=test_network, randomization_type="preserve_degree_distribution")
    print randomn.edges()

    print "preserve degree distribution and node degree:"
    randomn = randomize_graph(graph=test_network, randomization_type="preserve_degree_distribution_and_node_degree")
    print randomn.edges()

    print "creating big graph..."
    test_power_graph = networkx.Graph()
    test_power_graph.add_nodes_from(range(100000))
    test_power_graph.add_edges_from([ (random.randint(0,99999),random.randint(0,99999)) for x in xrange(1000000) ])
    print "randomizing big network by preserve_topology..."
    randomn = randomize_graph(graph=test_power_graph, randomization_type="preserve_topology")
    print "randomizing by preserve_topology_and_node_degree..."
    randomn = randomize_graph(graph=test_power_graph, randomization_type="preserve_topology_and_node_degree")   
    print "randomizing by preserve_degree_distribution..."
    randomn = randomize_graph(graph=test_power_graph, randomization_type="preserve_degree_distribution")


