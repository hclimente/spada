from spada.biological_entities.switch import IsoformSwitch
from spada.network.network import Network

import abc
import numpy as np
import operator
import pandas as pd
import random

class GeneNetwork(Network):
	"""docstring for GeneNetwork
	GeneNetwork contains a network of genes.

	Node information:
		id(str) 						Gene Id
		symbol(str) 					Gene Symbol
		switches(list,[]) 				List of detected isoform switches.
		driver(bool,False) 				Gene described as a driver.
		specificDriver(bool,False) 		Gene described as driver in this cancer type.
		driverType(str,"") 				Role that plays the gene in tumorigenesis.
		expressedTxsN(set,()) 			Set with transcripts with a median expression > 0.1
										in normal samples.
		expressedTxsT(set,()) 			Set with transcripts with a median expression > 0.1
										in tumor samples.
		neighborhoods(dictionary,{})	Adjusted p-value of differential expression.

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

	def add_node(self, full_name="", gene_id="?", gene_symbol="?"):
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
								switches 			= [],
								specificDriver 		= False,
								driver 				= False,
								driverType 			= None,
								expressedTxsN		= set(),
								expressedTxsT		= set(),
								neighborhoods		= {} )

			return True

	def update_node(self, key, value, full_name = "", gene_id = "", secondKey=""):
		"""Changes the value of a node attribute, specified by the key argument.
		Returns True if succesful; else, returns False."""

		geneID, geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id)
		finalValue = value

		if geneID is None:
			self.logger.error("Unable to get gene id from {} {}".format(full_name, geneID))
			return False

		return self._update_node(geneID,key,finalValue,secondKey)

	def update_nodes(self, key, values):
		for gene, value in values.items():
			if isinstance(value, set):
				for v in value:
					self.update_node(key, v, gene_id = gene)
			else:
				self.update_node(key, value, gene_id = gene)

	def add_edge(self, full_name1 = "", gene_id1 = "", symbol1 = "",
					   full_name2 = "", gene_id2 = "", symbol2 = ""):
		"""Adds an edge to the network. Return True if succesful; else, return False.
		The value of the attributes are the default, specified in GeneNetwork documentation."""

		id1 = [ x for x in [full_name1, gene_id1, symbol1] if x ]
		id2 = [ x for x in [full_name2, gene_id2, symbol2] if x ]

		if not id1 or not id2:
			self.logger.debug("Tried to add edge, but no node-information \
								 provided (Node 1[{}] and Node 2[{}]).".format(id1, id2))

		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1, gene_symbol=symbol1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2, gene_symbol=symbol2)[0]

		if (node_id1 is None or node_id1 is "") or (node_id2 is None or node_id2 is ""):
			self.logger.debug( "Cannot add edge {} - {}.".format(id1[0], id2[0]))
			return False
		elif node_id1 not in self.nodes():
			self.logger.debug("Node {} does not exist.".format(node_id1))
			return False
		elif node_id2 not in self.nodes():
			self.logger.debug("Node {} does not exist.".format(node_id2))
			return False

		return self._add_edge(node_id1, node_id2)

	def getSwitch(self, gene, nTx, tTx):

		thisSwitch = [ x for x in self.nodes()[gene]["switches"] if x.nTx == nTx and x.tTx == tTx ]

		if thisSwitch:
			return(thisSwitch[0])
		else:
		 	return(None)

	def update_edge(self, key, value, full_name1 = "", gene_id1 = "", full_name2 = "", gene_id2 = ""):
		"""Changes the value of an edge attribute, specified by the key argument.
		Returns True if succesful; else, returns False."""

		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2)[0]

		return self._update_edge(node_id1, node_id2, key, value)

	def flushSwitches(self):

		self.logger.debug("Cleaning imported network.")

		# removing isoform switches
		for gene,info in self.iterate_genes_byPatientNumber(alwaysSwitchedGenes=True):
			self._net.node[gene]["switches"] = []

	def readSwitches(self, switchesFile, tx_network):
		"""Import a set of genes with an isoform switch from candidateList.tsv.
		"""
		self.logger.debug("Retrieving calculated isoform switches.")

		switches = pd.read_csv(switchesFile, sep="\t")
		switches.columns = ["gene","normal","tumor","samples"]
		switches.samples = switches.samples.str.split(",")

		for index,row in switches.iterrows():

			if self.valid_switch(row["gene"], row["normal"], row["tumor"], tx_network):

				thisSwitch = IsoformSwitch(row["normal"], row["tumor"],row["samples"])
				nInfo = tx_network.nodes()[thisSwitch.nTx]
				tInfo = tx_network.nodes()[thisSwitch.tTx]
				thisSwitch.addTxInfo(nInfo, tInfo)

				self.update_node("switches", thisSwitch, full_name = row["gene"])

	def valid_switch(self, gene, nTx, tTx, tx_network):

		msg = ''

		if gene not in self.nodes():
			msg = "Gene {} not in the network. ".format(gene)

		if nTx not in tx_network.nodes():
			msg += "Transcript {} not in the network. ".format(nTx)
		if tTx not in tx_network.nodes():
			msg += "Transcript {} not in the network.".format(tTx)

		if msg:
			self.logger.warning('Switch {} - {} will not be analyzed. Reason(s): {}'.format(nTx, tTx, msg))
			return False
		else:
			return True

	def iterate_genes_byPatientNumber(self, onlySplicedGenes=True, onlyExpressedGenes=True, alwaysSwitchedGenes=False, bottom = 0, top = None):
		'''
		Iterate genes that have alternative splicing and more than one transcript expressed.
		'''

		if top == None:
			top = len(self.nodes()) - 1

		geneAndPatients = [ (g,sum([ len(s.samples) for s in self._net.node[g]["switches"] ])) for g in self.nodes() ][bottom:top]
		genes = [ g for g,n in sorted(geneAndPatients, key=operator.itemgetter(1), reverse=True) ]

		for gene in genes:
			info = self._net.node[gene]
			allExpressedTxs = set(info["expressedTxsN"]) | set(info["expressedTxsT"])
			self.logger.debug("Iterating gene {0}.".format(gene))
			if alwaysSwitchedGenes and info["switches"]:
				yield gene,info
			else:
				if onlySplicedGenes and len(allExpressedTxs) < 2 and not info["switches"]:
					continue
				if onlyExpressedGenes and not bool(allExpressedTxs):
					continue
				yield gene,info

	def iterate_switches_byPatientNumber(self, tx_network, only_models=False, relevance=None, removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most
				common) will be returned for each gene.
		"""

		counter = 1
		for gene,info in self.iterate_genes_byPatientNumber(alwaysSwitchedGenes=True):
			if not info["switches"]: continue

			switches = sorted(info["switches"], key = lambda a: len(a.samples),reverse=True)

			for thisSwitch in switches:
				if removeNoise and thisSwitch.isNoise:
					continue
				elif only_models and not thisSwitch.isMain:
					continue
				elif relevance is not None and relevance != thisSwitch.is_functional:
					continue

				self.logger.debug("Iterating switch number {}.".format(counter))
				counter += 1

				yield gene,info,thisSwitch

	def iterate_functionalSwitches_byPatientNumber(self, tx_network, only_models=False, removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most
				common) will be returned for each gene.
		"""

		self.iterate_switches_byPatientNumber(tx_network, only_models, True, removeNoise)

	def iterate_nonFunctionalSwitches_byPatientNumber(self, tx_network, only_models=False, removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most
				common) will be returned for each gene.
		"""

		self.iterate_switches_byPatientNumber(tx_network, only_models, False, removeNoise)

	def createSwitch(self,switchDict,tx_network,partialCreation):
		"""Create a switch object from the switch dictionary.

			partialCreation(bool): if False, the heavy protein
				objects are not created.
		"""
		thisSwitch = IsoformSwitch(switchDict["nIso"],switchDict["tIso"],switchDict["patients"])
		nInfo = tx_network._net.node[thisSwitch.nTx]
		tInfo = tx_network._net.node[thisSwitch.tTx]
		thisSwitch.addTxs(nInfo,tInfo)
		thisSwitch.addIsos(nInfo,tInfo,partialCreation)

		return thisSwitch

	def calculateCompatibilityTable(self):
		noise = []
		candidates = {}

		for gene,info in self.nodes(data=True):

			if not info["switches"]:
				continue

			nTxs = {}
			tTxs = {}
			for thisSwitch in info["switches"]:
				nTxs[thisSwitch.nTx] = nTxs.get(thisSwitch.nTx, 0) + len(thisSwitch.samples)
				tTxs[thisSwitch.tTx] = tTxs.get(thisSwitch.tTx, 0) + len(thisSwitch.samples)

			txs = set(nTxs.keys()) | set(tTxs.keys())

			scores = dict( (x, nTxs.get(x, 0) - tTxs.get(x, 0)) for x in txs )
			nTop = set([ x for x,s in scores.items() if max(scores.values()) == s ])
			tTop = set([ x for x,s in scores.items() if min(scores.values()) == s ])

			sortedSwitches = sorted(info["switches"], key=lambda a:len(a.samples), reverse=True)
			for thisSwitch in sortedSwitches:
				# only consider as top isoforms not present in the other condition
				if thisSwitch.nTx in (nTop - tTop) and thisSwitch.tTx in (tTop - nTop):
					candidates[gene] = thisSwitch
					break

			# no good switch could be found
			if not gene in candidates:
				noise.extend([ len(x.samples) for x in info["switches"] ])
			else:
				noise.extend([ len(x.samples) for x in info["switches"] if x.nTx == candidates[gene].tTx or x.tTx == candidates[gene].nTx ])

		cutoff = 0 if not noise else np.percentile(noise, 99)

		for gene,info in self.nodes(data=True):
			if not info["switches"]: continue

			for s in info["switches"]:
				s.setNoise( len(s.samples) <= cutoff )
				s.setMain( gene in candidates and s == candidates[gene] )

	def sampleSwitches(self,tx_network,partialCreation=True,numIterations=2000):

		genesWithSwitches = [ gene for gene,info in self.iterate_genes_byPatientNumber(alwaysSwitchedGenes=True) if info["switches"] ]
		genes = random.sample(genesWithSwitches,numIterations)

		for gene in genes:
			info = self._net.node[gene]
			switchDict = random.choice(info["switches"])
			thisSwitch = self.createSwitch(switchDict,tx_network,partialCreation)

			yield gene,info,switchDict,thisSwitch

	def isDriver(self, gene):

		driver = 'No'
		if self._net.node[gene]["specificDriver"]:
			driver = 'Tumor-specific_driver'
		elif self._net.node[gene]["driver"]:
			driver = 'Foreign_driver'
		elif [ x for x in self._net.neighbors(gene) if self._net.node[x]["driver"] ]:
			driver = 'Driver_interactor'

		return driver
