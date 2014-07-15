#!/soft/devel/python-2.7/bin/python

import networkx as nx

class Network:
	def __init__(self):
		self._net 			= nx.Graph()
		self._rejectedNodes = set()
	
	## Getters ##
	def n(self): 		return self._net
	
	## Pseudo-getters ##
	def nodes(self, **kwds): 	return self._net.nodes(**kwds)
	def edges(self, **kwds):	return self._net.edges(**kwds)

	def _update_node(self, nodeId, key, value):
		if isinstance( self._net.node[nodeId][key], list):
			self._net.node[nodeId][key].append( value )
		elif isinstance( self._net.node[nodeId][key], set):
			self._net.node[nodeId][key].add( value )
		else: 
			self._net.node[nodeId][key] = value

	def _update_edge(self, node1, node2, key, value):
		if isinstance( self._net.node[nodeId][key], list):
			self._net.edge[node1][node2][key].append( value )
		elif isinstance( self._net.node[nodeId][key], set):
			self._net.edge[node1][node2][key].add( value )
		else: 
			self._net.edge[node1][node2][key] = value