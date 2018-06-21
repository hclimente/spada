from spada.biological_entities.switch import IsoformSwitch, LiteSwitch
from spada.io import io
from spada.network.network import Network

import abc
import numpy as np
import operator
import random

class GeneNetwork(Network):
	"""docstring for GeneNetwork
	GeneNetwork contains a network of genes.

	Node information:
		id(str) 						Gene Id
		symbol(str) 					Gene Symbol
		switches(list,[]) 				List of detected isoform switches.

	Edge information:
		id1(str) 						Gene id of interactor 1.
		id2(str) 						Gene id of interactor 2.

	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):
		Network.__init__(self, name)

	@abc.abstractmethod
	def nameFilter(self, **kwds):
		"""Receive a gene identifier and convert it to the consensus identifier for the network."""
		raise NotImplementedError()

	def add_node(self, full_name = "", gene_id = "", gene_symbol = ""):
		"""Adds a node to the network. Return True if succesful; else, return False.
		The value of the attributes are the default, specified in GeneNetwork documentation."""

		self.logger.debug("Importing node: full name {0} id {1} symbol {2}".format(
								full_name, gene_id, gene_symbol) )
		geneID,geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id, gene_symbol=gene_symbol)

		if geneID in self._net.nodes():
			self.logger.debug("Node {0} already exist.".format(geneID))
			return True
		elif geneID is None:
			self.logger.debug("Could not retrieve name from node: full name {0} id {1} symbol {2}".format(
								full_name, gene_id, gene_symbol) )
			return False
		else:
			self.logger.debug("Node {0} imported.".format(geneID))
			self._net.add_node( geneID,
								symbol 				= geneSymbol,
								switches 			= [] )

			return True

	def update_node(self, key, value, full_name = "", gene_id = "", secondKey=""):
		"""Changes the value of a node attribute, specified by the key argument.
		Returns True if succesful; else, returns False."""

		geneID, geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id)

		if geneID is None:
			self.logger.error("Unable to get gene id from {} {}".format(full_name, geneID))
			return False

		return self._update_node(geneID,key,value,secondKey)

	def update_nodes(self, key, values):
		for gene, value in values.items():
			self.update_node(key, value, gene_id = gene)

	def add_edge(self, full_name1 = "", gene_id1 = "", symbol1 = "",
					   full_name2 = "", gene_id2 = "", symbol2 = ""):
		"""Adds an edge to the network. Return True if succesful; else, return False.
		The value of the attributes are the default, specified in GeneNetwork documentation."""

		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1, gene_symbol=symbol1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2, gene_symbol=symbol2)[0]

		if (node_id1 is None or node_id1 is "") or (node_id2 is None or node_id2 is ""):
			self.logger.debug( "Cannot add edge {} - {}.".format(full_name1, full_name2))
			return False
		elif node_id1 not in self.nodes():
			self.logger.debug("Node {} does not exist.".format(node_id1))
			return False
		elif node_id2 not in self.nodes():
			self.logger.debug("Node {} does not exist.".format(node_id2))
			return False

		return self._add_edge(node_id1, node_id2)

	def flushSwitches(self):

		self.logger.debug("Cleaning imported network.")

		# removing isoform switches
		for gene,info in self.genes():
			self._net.node[gene]["switches"] = []

	def readSwitches(self, switchesFile, tx_network):
		"""Import a set of genes with an isoform switch from candidateList.tsv.
		"""
		self.logger.debug("Retrieving calculated isoform switches.")

		for line in io.readTable(switchesFile):
			gene = line['GeneId']
			ctrl = line['Control_transcript']
			case = line['Case_transcript']
			samples = set(line['Samples'].split(','))

			if self.valid_switch(gene, ctrl, case, tx_network):

				thisSwitch = LiteSwitch(ctrl, case, samples)
				self.update_node("switches", thisSwitch, full_name = gene)

	def valid_switch(self, gene, ctrl, case, tx_network):

		msg = ''

		if gene not in self.nodes():
			msg = "Gene {} not in the network. ".format(gene)

		if ctrl not in tx_network.nodes():
			msg += "Transcript {} not in the network. ".format(ctrl)
		if case not in tx_network.nodes():
			msg += "Transcript {} not in the network.".format(case)

		if msg:
			self.logger.warning('Switch {} - {} will not be analyzed. Reason(s): {}'.format(ctrl, case, msg))
			return False
		else:
			return True

	def genes(self):
		'''
		Iterate genes that have alternative splicing and more than one transcript expressed.
		'''

		geneAndPatients = [ (g,sum([ len(s.samples) for s in self._net.node[g]["switches"] ])) for g in self.nodes() ]
		genes = [ g for g,n in sorted(geneAndPatients, key=operator.itemgetter(1), reverse=True) ]

		for gene in genes:
			yield gene,self._net.node[gene]

	def switches(self, txs):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most
				common) will be returned for each gene.
		"""

		for gene,info in self.genes():
			if not info["switches"]: continue

			switches = sorted(info["switches"], key = lambda a: len(a.samples), reverse = True)

			for switch in switches:

				thisSwitch = self.__createSwitch(switch, txs)
				yield gene, info, thisSwitch

	def __createSwitch(self, switch, tx_network):
		"""Create a switch object from the switch dictionary.

			partialCreation(bool): if False, the heavy protein
				objects are not created.
		"""
		thisSwitch = IsoformSwitch(switch.ctrl, switch.case, switch.samples)
		nInfo = tx_network._net.node[thisSwitch.ctrl]
		tInfo = tx_network._net.node[thisSwitch.case]
		thisSwitch.addTxInfo(nInfo,tInfo)

		return thisSwitch

	def sampleSwitches(self, tx_network, numIterations = 2000):

		genesWithSwitches = [ gene for gene,info in self.genes() if info["switches"] ]
		genes = random.sample(genesWithSwitches,numIterations)

		for gene in genes:
			info = self._net.node[gene]
			switchDict = random.choice(info["switches"])
			thisSwitch = self.__createSwitch(switchDict, tx_network)

			yield gene,info,switchDict,thisSwitch
