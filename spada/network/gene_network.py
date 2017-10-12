from spada.biological_entities import switch
from spada import utils
from spada.network import network

import abc
import numpy as np
import operator
import pandas as pd
import random

class GeneNetwork(network.Network):
	"""docstring for GeneNetwork
	GeneNetwork contains a network of genes.

	Node information:
		id(str) 						Gene Id
		symbol(str) 					Gene Symbol
		isoformSwitches(list,[]) 		List of detected isoform switches [[isoN, isoT], [isoN, isoT]]
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
		network.Network.__init__(self, name)

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
								isoformSwitches 	= [],
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

		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1, gene_symbol=symbol1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2, gene_symbol=symbol2)[0]

		if (node_id1 is None or node_id1 is "") or (node_id2 is None or node_id2 is ""):
			self.logger.warning( "Cannot add edge {} - {} ({} - {}).".format(
									full_name1, gene_id1, full_name2, gene_id2) )
			return False
		elif node_id1 not in self.nodes():
			self.logger.warning("Node {} does not exist.".format(node_id1))
			return False
		elif node_id2 not in self.nodes():
			self.logger.warning("Node {} does not exist.".format(node_id2))
			return False

		return self._add_edge(node_id1, node_id2)


	def update_edge(self, key, value, full_name1 = "", gene_id1 = "", full_name2 = "", gene_id2 = ""):
		"""Changes the value of an edge attribute, specified by the key argument.
		Returns True if succesful; else, returns False."""

		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2)[0]

		return self._update_edge(node_id1, node_id2, key, value)

	def cleanNetwork(self):

		self.logger.debug("Cleaning imported network.")

		# removing isoform switches
		for gene,info in self.iterate_genes_byPatientNumber(alwaysSwitchedGenes=True):
			self._net.node[gene]["isoformSwitches"] = []

	def importCandidates(self,candidatesFile):
		"""Import a set of genes with an isoform switch from candidateList.tsv.
		"""
		self.logger.debug("Retrieving calculated isoform switches.")

		switches = pd.DataFrame.from_csv(candidatesFile, sep="\t", header=None, index_col=None)
		switches.columns = ["Gene","Transcript_normal","Transcript_tumor","Samples"]
		switches.Samples = switches.Samples.str.split(",")

		for index,row in switches.iterrows():
			Gene = row["Gene"]

			switch = {}
			switch["nIso"] 			= row["Transcript_normal"]
			switch["tIso"] 			= row["Transcript_tumor"]
			switch["patients"] 		= row["Samples"]
			switch["functional"] 	= None
			switch["model"] 		= None
			switch["noise"] 		= None

			self.update_node("isoformSwitches",switch,full_name=Gene)

	def iterate_genes_byPatientNumber(self,onlySplicedGenes=True,onlyExpressedGenes=True,alwaysSwitchedGenes=False):
		'''
		Iterate genes that have alternative splicing and more than one transcript expressed.
		'''
		consideredGenes = self.nodes()
		if options.Options().parallelRange:
			bottom = options.Options().parallelRange - 1
			top = bottom + options.Options().step
			consideredGenes = consideredGenes[bottom:top]

		geneAndPatients = [ (x,sum([ len(z["patients"]) for z in self._net.node[x]["isoformSwitches"] ])) for x in consideredGenes ]
		genes = [ x for x,y in sorted(geneAndPatients, key=operator.itemgetter(1), reverse=True) ]

		for gene in genes:
			info = self._net.node[gene]
			allExpressedTxs = set(info["expressedTxsN"]) | set(info["expressedTxsT"])
			self.logger.debug("Iterating gene {0}.".format(gene))
			if alwaysSwitchedGenes and info["isoformSwitches"]:
				yield gene,info
			else:
				if onlySplicedGenes and len(allExpressedTxs) < 2 and not info["isoformSwitches"]:
					continue
				if onlyExpressedGenes and not bool(allExpressedTxs):
					continue
				yield gene,info

	def iterate_switches_byPatientNumber(self,tx_network,only_models=False,relevance=None,partialCreation=False,removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most
				common) will be returned for each gene.
		"""

		counter = 1
		for gene,info in self.iterate_genes_byPatientNumber(alwaysSwitchedGenes=True):
			if not info["isoformSwitches"]: continue

			switches = sorted(info["isoformSwitches"],key=lambda a:len(a['patients']),reverse=True)

			for switchDict in switches:
				if removeNoise and switchDict["noise"]:
					continue
				elif only_models and not switchDict["model"]:
					continue

				thisSwitch = self.createSwitch(switchDict,tx_network,partialCreation)

				if relevance is not None and relevance != thisSwitch.is_functional:
					continue

				self.logger.debug("Iterating switch number {0}.".format(counter))
				counter += 1

				yield gene,info,switchDict,thisSwitch

	def iterate_functionalSwitches_byPatientNumber(self,tx_network,only_models=False,partialCreation=False,removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most
				common) will be returned for each gene.
		"""

		self.iterate_switches_byPatientNumber(tx_network,only_models,True,partialCreation,removeNoise)

	def iterate_nonFunctionalSwitches_byPatientNumber(self,tx_network,only_models=False,partialCreation=False,removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most
				common) will be returned for each gene.
		"""

		self.iterate_switches_byPatientNumber(tx_network,only_models,False,partialCreation,removeNoise)

	def createSwitch(self,switchDict,tx_network,partialCreation):
		"""Create a switch object from the switch dictionary.

			partialCreation(bool): if False, the heavy protein
				objects are not created.
		"""
		thisSwitch = switch.IsoformSwitch(switchDict["nIso"],switchDict["tIso"],switchDict["patients"])
		nInfo = tx_network._net.node[thisSwitch.nTx]
		tInfo = tx_network._net.node[thisSwitch.tTx]
		thisSwitch.addTxs(nInfo,tInfo)
		thisSwitch.addIsos(nInfo,tInfo,partialCreation)

		return thisSwitch

	def calculateCompatibilityTable(self):
		incompatible = []
		bestpatible = {}

		for gene,info in self.nodes(data=True):
			if not info["isoformSwitches"]:
				continue

			txs = []
			[ txs.extend([x["nIso"],x["tIso"]]) for x in info["isoformSwitches"] ]

			d = dict.fromkeys(set(txs), {"N":0,"T":0} )

			for x in info["isoformSwitches"]:
				d[x["nIso"]] = { "N":d[x["nIso"]]["N"]+len(x["patients"]), "T":d[x["nIso"]]["T"] }
				d[x["tIso"]] = { "T":d[x["tIso"]]["T"]+len(x["patients"]), "N":d[x["tIso"]]["N"] }

			scores = [ d[x]["N"] - d[x]["T"] for x in d ]
			consensus = { 'N': [ x for x in d if max(scores) == (d[x]["N"]-d[x]["T"]) ],
						  'T': [ x for x in d if min(scores) == (d[x]["N"]-d[x]["T"]) ] }

			sortedSwitches = sorted(info["isoformSwitches"],key=lambda a:len(a['patients']),reverse=True)
			for x in sortedSwitches:
				if x["nIso"] in consensus["N"] and x["tIso"] in consensus["T"]:
					bestpatible[gene] = {"N":x["nIso"], "T":x["tIso"], "n":len(x["patients"])}
					break

			# no good switch could be found
			if not gene in bestpatible:
				continue

			incompatible.append(max([max( d[x]["N"] for x in d if x == bestpatible[gene]["T"] ),max( d[x]["T"] for x in d if x == bestpatible[gene]["N"] )]))

		threshold = np.percentile(incompatible,99)

		for gene,info in self.nodes(data=True):
			if not info["isoformSwitches"]: continue

			for s in info["isoformSwitches"]:
				if len(s["patients"]) <= threshold:
					s["noise"] = True
				else:
					s["noise"] = False

				if gene in bestpatible:
					if s["nIso"] == bestpatible[gene]["N"] and s["tIso"] == bestpatible[gene]["T"]:
						s["model"] = True
					else:
						s["model"] = False
				else:
					s["model"] = False

	def sampleSwitches(self,tx_network,partialCreation=True,numIterations=2000):

		genesWithSwitches = [ gene for gene,info in self.iterate_genes_byPatientNumber(alwaysSwitchedGenes=True) if info["isoformSwitches"] ]
		genes = random.sample(genesWithSwitches,numIterations)

		for gene in genes:
			info = self._net.node[gene]
			switchDict = random.choice(info["isoformSwitches"])
			thisSwitch = self.createSwitch(switchDict,tx_network,partialCreation)

			yield gene,info,switchDict,thisSwitch

	def getGeneAnnotation(self,gene,hallmarksDict,biologicalProcessDict):

		annotation = "Nothing"
		driverAnnotation = "Nothing"

		if self._net.node[gene]["driver"]:
			driverAnnotation = "driver"
		elif [ x for x in self._net.neighbors(gene) if self._net.node[x]["driver"] ]:
			driverAnnotation = "d1"

		if [ x for x in hallmarksDict if gene in hallmarksDict[x] ]:
			annotation = ",".join(sorted([ x for x in hallmarksDict if gene in hallmarksDict[x] ]))
		elif [ x for x in biologicalProcessDict if gene in biologicalProcessDict[x] ]:
			annotation = ",".join(sorted([ x for x in biologicalProcessDict if gene in biologicalProcessDict[x] ]))

		return (annotation,driverAnnotation)
