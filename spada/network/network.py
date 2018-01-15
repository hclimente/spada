import networkx as nx
import pickle
import abc
import logging

class Network:

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):

		self._name			= name
		self._net 			= nx.Graph()
		self._rejectedNodes = set()
		self.createLogger()

	## Getters ##
	def n(self): 		return self._net

	## Pseudo-getters ##
	def nodes(self, **kwds): 	return self._net.nodes(**kwds)
	def edges(self, **kwds):	return self._net.edges(**kwds)

	def _update_node(self, nodeId, key, value, secondKey=""):

		if nodeId not in self.nodes():
			self.logger.warning("Tried to update node {}, but does not exist.".format(nodeId))
			return False

		if isinstance( self._net.node[nodeId][key], list):
			self._net.node[nodeId][key].append( value )
		elif isinstance( self._net.node[nodeId][key], set):
			self._net.node[nodeId][key].add( value )
		elif isinstance( self._net.node[nodeId][key], dict):
			if secondKey not in self._net.node[nodeId][key]:
				self._net.node[nodeId][key].setdefault(secondKey,set())
			self._net.node[nodeId][key][secondKey].add(value)
		else:
			if self._net.node[nodeId] is not None and self._net.node[nodeId] != value:
				self.logger.debug("Node {}, {} had a value of {}. Updated to {}.".format(
									nodeId, key, self._net.node[nodeId][key], value) )
			self._net.node[nodeId][key] = value

		return True

	def _update_edge(self, node1, node2, key, value):

		if not [ (x,y) for (x,y) in self.edges() if (node1,node2) == (x,y) or (node1,node2) == (y,x) ]:
			self.logger.warning("Tried to update edge {} - {}, but it does not exist.".format(node1, node2))
			return False

		if isinstance( self._net.edge[node1][node2][key], list):
			self._net.edge[node1][node2][key].append( value )
		elif isinstance( self._net.edge[node1][node2][key], set):
			self._net.edge[node1][node2][key].add( value )
		else:
			if self._net.edge[node1][node2][key] is not None and self._net.edge[node1][node2][key] != value:
				self.logger.debug("Edge {} - {}, {} had a value of {}. Updated to {} .".format(
								node1, node2, key, self._net.edge[node1][node2][key], value) )

			self._net.edge[node1][node2][key] = value

		return True

	def _add_edge(self, node1, node2, **kwds):

		if [ (x,y) for (x,y) in self.edges() if (node1,node2) == (x,y) or (node1,node2) == (y,x) ]:
			self.logger.debug("Tried to add edge {} - {}, but it exists.".format(node1, node2))
			return False

		self._net.add_edge(node1, node2, **kwds)

		return True

	def saveNetwork(self,filename):

		self.logger.debug("Saving network at {}.".format(filename))
		#Unattach logger to save without thread problems
		self.removeLogger()
		with open(filename, "wb") as NET_DUMP:
			pickle.dump(self, NET_DUMP, -1)

		self.createLogger()

	def createLogger(self):
		self.logger	= logging.getLogger(self._name)
	def removeLogger(self):
		self.logger = None