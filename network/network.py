#!/soft/devel/python-2.7/bin/python

from libs import options

import networkx
import cPickle
import abc
import logging

class Network:

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):
		
		self._name			= name
		self._net 			= networkx.Graph()
		self._rejectedNodes = set()
		self.createLogger()
	
	## Getters ##
	def n(self): 		return self._net
	
	## Pseudo-getters ##
	def nodes(self, **kwds): 	return self._net.nodes(**kwds)
	def edges(self, **kwds):	return self._net.edges(**kwds)

	def _update_node(self, nodeId, key, value):

		finalValue = value

		if nodeId not in self.nodes():
			logging.warning("Tried to update node {0}, but does not exist.".format(nodeId))
			return False

		if isinstance( self._net.node[nodeId][key], list):
			self._net.node[nodeId][key].append( finalValue )
		elif isinstance( self._net.node[nodeId][key], set):
			self._net.node[nodeId][key].add( finalValue )
		else: 
			if key is "score":
				finalValue = min(1.0, self._net.node[nodeId]["score"] + finalValue)
			if self._net.node[nodeId] is not None and self._net.node[nodeId] != finalValue:
				logging.debug("Node {0}, {1} had a value of {2}. Updated to {3}.".format(
									nodeId, key, self._net.node[nodeId][key], finalValue) )
			self._net.node[nodeId][key] = finalValue

		return True

	def _update_edge(self, node1, node2, key, value):

		finalValue = value

		if not [ (x,y) for (x,y) in self.edges() if (node1,node2) == (x,y) or (node1,node2) == (y,x) ]:
			logging.warning("Tried to update edge {0} - {1}, but it does not exist.".format(node1, node2))
			return False

		if isinstance( self._net.edge[node1][node2][key], list):
			self._net.edge[node1][node2][key].append( finalValue )
		elif isinstance( self._net.edge[node1][node2][key], set):
			self._net.edge[node1][node2][key].add( finalValue )
		else: 
			if key is "score":
				finalValue = min(1.0, self._net.edge[node1][node2]["score"] + finalValue)
			if self._net.edge[node1][node2][key] is not None and self._net.edge[node1][node2][key] != finalValue:
				logging.debug("Edge {0} - {1}, {2} had a value of {3}. Updated to {4} .".format(
								node1, node2, key, self._net.edge[node1][node2][key], finalValue) )

			self._net.edge[node1][node2][key] = finalValue

		return True

	def _add_edge(self, node1, node2, **kwds):

		if [ (x,y) for (x,y) in self.edges() if (node1,node2) == (x,y) or (node1,node2) == (y,x) ]:
			logging.debug("Tried to add edge {0} - {1}, but it exists.".format(node1, node2))
			return False

		self._net.add_edge(node1, node2, **kwds)

		return True

	def saveNetwork(self, filename):
		#Unattach logger to allow saving
		self.removeLogger()
		logging.debug("Saving network at {0}{1}.".format(options.Options().qout,filename))
		with open(options.Options().qout + filename, "wb") as NET_DUMP:
			cPickle.dump(self, NET_DUMP, -1)

		self.createLogger()

	def createLogger(self):
		self.logger	= logging.getLogger(self._name)
	def removeLogger(self):
		self.logger = None