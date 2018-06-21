import networkx as nx
import pickle
import abc
import logging

class Network:

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):

		self._name	= name
		self._net 	= nx.Graph()
		self.createLogger()

	@abc.abstractmethod
	def accept(self, **kwds):
		raise NotImplementedError()

	## Getters ##
	def n(self): 		return self._net

	## Pseudo-getters ##
	def nodes(self, **kwds): 	return self._net.nodes(**kwds)
	def edges(self, **kwds):	return self._net.edges(**kwds)

	def _update_node(self, node, key, value, secondKey = "", override = False):

		if node not in self.nodes():
			self.logger.debug("Tried to update node {}, but does not exist.".format(node))
			return False

		if not override and isinstance(self._net.node[node][key], list):
			self._net.node[node][key].append( value )
		elif not override and isinstance(self._net.node[node][key], set):
			self._net.node[node][key].add( value )
		elif not override and isinstance(self._net.node[node][key], dict):
			if secondKey not in self._net.node[node][key]:
				self._net.node[node][key].setdefault(secondKey,set())
			self._net.node[node][key][secondKey].add(value)
		else:
			if self._net.node[node] is not None and self._net.node[node] != value:
				self.logger.debug("Node {}, {} had a value of {}. Updated to {}.".format(
									node, key, self._net.node[node][key], value) )
			self._net.node[node][key] = value

		return True

	def _add_edge(self, node1, node2, **kwds):

		if [ (x,y) for (x,y) in self.edges() if (node1,node2) == (x,y) or (node1,node2) == (y,x) ]:
			self.logger.debug("Tried to add edge {} - {}, but it exists.".format(node1, node2))
			return False

		self._net.add_edge(node1, node2, **kwds)

		return True

	def createLogger(self):
		self.logger	= logging.getLogger(self._name)
	def removeLogger(self):
		self.logger = None
