import network
from libs import utils
from libs import options
from biological_entities import switch

import abc
import biana
import numpy as np
import pandas as pd
import random
import subprocess

class GeneNetwork(network.Network):
	"""docstring for GeneNetwork
	GeneNetwork contains a network of genes.

	Node information:
		id(str) 						Gene Id
		symbol(str) 					Gene Symbol
		isoformSwitches(list,[]) 		List of detected isoform switches [[isoN, isoT], [isoN, isoT]]
		driver(bool,False) 				Gene described as a driver.
		asDriver(bool,False) 			Gene described as driver through alternative splicing.
		druggable(bool,False) 			Gene described as druggable.
		specificDriver(bool,False) 		Gene described as driver in this cancer type.
		driverType(str,"") 				Role that plays the gene in tumorigenesis.
		expressedTxsNormal(set,()) 		Set with transcripts with a median expression > 0.1 
										in normal samples.
		expressedTxsTumor(set,()) 		Set with transcripts with a median expression > 0.1 
										in tumor samples.
		neighborhoods(dictionary,{})	Adjusted p-value of differential expression.

	Edge information:
		id1(str) 						Gene id of interactor 1.
		id2(str) 						Gene id of interactor 2.
		experimental(bool,None) 		Interaction found through experiments.

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
								asDriver 			= False,
								druggable 			= False, 
								expressedTxsNormal	= set(),
								expressedTxsTumor	= set(),
								neighborhoods		= {} )

			return True

	def update_node(self, key, value,full_name = "",gene_id = "",secondKey=""):
		"""Changes the value of a node attribute, specified by the key argument. 
		Returns True if succesful; else, returns False."""

		geneID, geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id)
		finalValue = value
		
		if geneID is None:
			self.logger.error("Unable to get gene id from {0} {1}".format(full_name, geneID))
			return False

		return self._update_node(geneID,key,finalValue,secondKey)

	def add_edge(self, full_name1 = "", gene_id1 = "", full_name2 = "", gene_id2 = ""):
		"""Adds an edge to the network. Return True if succesful; else, return False.
		The value of the attributes are the default, specified in GeneNetwork documentation."""

		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2)[0]

		if (node_id1 is None or node_id1 is "") or (node_id2 is None or node_id2 is ""): 
			self.logger.warning( "Cannot add edge {1} - {3} ({0} - {2}).".format(
									full_name1, gene_id1, full_name2, gene_id2) )
			return False
		elif node_id1 not in self.nodes():
			self.logger.warning("Node {0} does not exist.".format(node_id1))
			return False
		elif node_id2 not in self.nodes():
			self.logger.warning("Node {0} does not exist.".format(node_id2))
			return False

		return self._add_edge(	node_id1, 
								node_id2, 
								experimental = None)

	def update_edge(self, key, value, full_name1 = "", gene_id1 = "", full_name2 = "", gene_id2 = ""):
		"""Changes the value of an edge attribute, specified by the key argument. 
		Returns True if succesful; else, returns False."""
		
		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2)[0]

		return self._update_edge(node_id1, node_id2, key, value)

	def readGeneInfo(self):
		"""Read tsv files containing characteristics of the genes. Updates the nodes.:
			- compilationTable.tsv: gene annotation of drivers, epigenetic factors and RBPs.
			- expressedGenes.lst: R-generated file containing the transcripts above
			a threshold of expression.
		"""
		
		# expression info
		for line in utils.readTable(options.Options().qout + "transcript_expression.tsv"):
			gene = line[0]
			tx = line[1]
			expressedNormal = float(line[2]) > 0.1
			expressedTumor = float(line[3]) > 0.1

			if self.add_node(full_name=gene):
				if expressedNormal:
					self.update_node("expressedTxsNormal",tx,full_name=gene)
				if expressedTumor:
					self.update_node("expressedTxsTumor",tx,full_name=gene)

		# druggability info
		for line in utils.readTable("data/Databases/dgidb_export_all_drivers_bygene_results.tsv"):
			geneSymbol = line[0]
				
			for gene,info in self.nodes(data=True):
				if info["symbol"] == geneSymbol:
					self.update_node("druggable",True,gene_id=gene )
					break

		# driver info
		## is driver: only for ucsc
		for line in utils.readTable("data/ucsc/drivers_ucsc_notation.txt",sep="|"):
			geneid = line[1]

			self.update_node("driver",True,gene_id=geneid)

		## which kind of driver
		for line in utils.readTable("data/Databases/cancer_networks_SuppTables_v7_S7.csv"):
			geneSymbol = line[0]
			role = line[1]
			
			for gene,info in self.nodes(data=True):
				if info["symbol"] == geneSymbol:
					if info["driver"]:
						self.update_node("driverType",role,gene_id=gene )
					break

		## AS driver: only for ucsc
		for line in utils.readTable("data/ucsc/asdrivers_ucsc_notation.txt",sep="|"):
			geneid = line[1]

			self.update_node("asDriver",True,gene_id=geneid)

	def importExternalCandidates(self):
		
		self.logger.debug("Retrieving externally calculated isoform switches.")

		for line in utils.readTable(options.Options().externalSwitchesFile):
			Gene = line[0]

			switch = {}
			switch["tIso"] 			= line[1]
			switch["nIso"] 			= line[2]
			switch["patients"] 		= []
			switch["functional"] 	= None
			switch["model"] 		= True
			switch["noise"] 		= False

			self.update_node("isoformSwitches",switch,full_name=Gene)

	def cleanNetwork(self):
		
		self.logger.debug("Cleaning imported network.")

		# removing isoform switches
		for gene,info in self.iterate_genes_ScoreWise():
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

	def importSpecificDrivers(self):
		self.logger.debug("Importing specific drivers.")

		for nameComponents in utils.readTable(options.Options().specificDrivers, header=False):

			geneID = self.nameFilter(gene_symbol = nameComponents[0])[0]

			if not geneID: continue

			self.logger.debug("Adding {0} as specific driver.".format(geneID))

			self.update_node("specificDriver", True, gene_id = geneID)
			self.update_node("driver", True, gene_id = geneID)
			self.update_node("score", 1, gene_id = geneID)

	def importKnownInteractions(self):

		self.logger.debug("Importing interactions from BIANA.")

		affinity_methods = { 
							'492':'in vivo', '493':'in vitro', '0':'molecular interaction', 
							'4':'affinity chromatography technology', '6':'anti bait coimmunoprecipitation', 
							'7':'anti tag coimmunoprecipitation', '8':'array technology', 
							'9':'bacterial display', '19':'coimmunoprecipitation', 
							'28':'cosedimentation in solution', '29':'cosedimentation through density gradient', 
							'30':'cross-linking', '34':'display technology', '47':'far western blotting', 
							'48':'filamentous phage display', '49':'filter binding', '66':'lambda phage display', 
							'71':'molecular sieving', '73':'mrna display', '81':'peptide array', 
							'84':'phage display', '89':'protein array', '92':'protein in situ array', 
							'95':'proteinchip(r) on a surface-enhanced laser desorption/ionization', 
							'96':'pull down', '98':'ribosome display', '108':'t7 phage display', 
							'115':'yeast display', '225':'chromatin immunoprecipitation array', 
							'400':'affinity technology', '402':'chromatin immunoprecipitation assay', 
							'405':'competition binding', '411':'enzyme linked immunosorbent assay', 
							'412':'electrophoretic mobility supershift assay', 
							'413':'electrophoretic mobility shift assay', '440':'saturation binding', 
							'657':'systematic evolution of ligands by exponential enrichment', 
							'676':'tandem affinity purification', '678':'antibody array', 
							'695':'sandwich immunoassay', '729':'luminescence based mammalian interactome mapping', 
							'813':'proximity enzyme linked immunosorbent assay', 
							'858':'immunodepleted coimmunoprecipitation', '892':'solid phase assay', 
							'899':'p3 filamentous phage display', '900':'p8 filamentous phage display', 
							'921':'surface plasmon resonance array', '946':'ping', '947':'bead aggregation assay', 
							'963':'interactome parallel affinity capture', '1017':'rna immunoprecipitation', 
							'1028':'modified chromatin immunoprecipitation', 
							'1029':'proteomics of isolated chromatin segments', '1031':'protein folding/unfolding', 
							'1087':'monoclonal antibody blockade'
						}

		complementation_methods={
    			'492':'in vivo', '493':'in vitro', '0':'molecular interaction',
    			'10':	'beta galactosidase complementation', '11':	'beta lactamase complementation', 
    			'14':	'adenylate cyclase complementation', '18':	'two hybrid', 
    			'90':	'protein complementation assay', '97':	'reverse ras recruitment system',
    			'111':	'dihydrofolate reductase reconstruction', '112':	'ubiquitin reconstruction', 
    			'228':	'cytoplasmic complementation assay', 
    			'229':	'green fluorescence protein complementation assay', 
    			'230':	'membrane bound complementation assay', 
    			'231':	'mammalian protein protein interaction trap',
    			'232':	'transcriptional complementation assay', '369':	'lex-a dimerization assay', 
    			'370':	'tox-r dimerization assay', '397':	'two hybrid array', 
    			'398':	'two hybrid pooling approach', '399':	'two hybrid fragment pooling approach',
    			'432':	'one hybrid', '437':	'protein tri hybrid', '438':	'rna tri hybrid',
    			'588':	'3 hybrid method', '655':	'lambda repressor two hybrid', 
    			'726':	'reverse two hybrid', '727':	'lexa b52 complementation', 
    			'728':	'gal4 vp16 complementation', '809':	'bimolecular fluorescence complementation',
    			'895':	'protein kinase A complementation', '916':	'lexa vp16 complementation', 
    			'1037':	'Split renilla luciferase complementation',
    		}
		
		methods  = [("Method_id",18),("Method_id",696)]
		methods.extend([ ("Method_id", x) for x in complementation_methods if x not in affinity_methods ])

		seeds = [ x for x,props in self.nodes(data=True) if props["isoformSwitches"] or props["specificDriver"] ]

		bianaInputType = "geneid"

		if options.Options().annotation == "ensembl": bianaInputType = "ensembl"

		session = biana.create_new_session(
										sessionID="SmartAS", 
										dbname="BIANA_MARCH_2013", 
										dbhost="ben-yehuda",
										dbuser="biana_user", 
										dbpassword="biana_password",
										unification_protocol="uniprot_geneID_seqtax"
									)

		proteome = session.create_new_user_entity_set(
												identifier_description_list = seeds,
												attribute_restriction_list 	= [("taxid", "9606")],
												id_type 					= bianaInputType,
												new_user_entity_set_id		= "proteome",
													  )
		session.create_network( 
								user_entity_set_id 					= "proteome", 
								level 								= 5, 
								relation_type_list 					= ["interaction"],
								relation_attribute_restriction_list = methods,
								include_relations_last_level 		= False, #Seguro?
								use_self_relations 					= False
							  )

		#Iterate through all the interactions
		c = 1
		for (userEntity_id1, userEntity_id2) in proteome.getRelations():
			self.logger.debug( "Interaction {0}/{1}".format(c, len(proteome.getRelations())) )
			c += 1
			eErIDs_list 		= proteome.get_external_entity_relation_ids(
																userEntity_id1, userEntity_id2)
			method_names 		= set()
			method_ids 			= set()
			use_method_ids 		= set()
			relationObj_dict 	= session.dbAccess.get_external_entities_dict(
												externalEntityIdsList 		= eErIDs_list, 
												attribute_list 				= [],
												relation_attribute_list 	= ["method_id","psimi_name"], 
												participant_attribute_list 	= []
																			  )
			
			if session.get_defined_node_attributes("proteome", userEntity_id1, bianaInputType, True):
				geneId_1 = session.get_defined_node_attributes("proteome", userEntity_id1, bianaInputType, True).pop()
			else:
				self.logger.warning("No {0} id for user entity {1}.".format(bianaInputType, userEntity_id1))
				continue
			
			if session.get_defined_node_attributes("proteome", userEntity_id2, bianaInputType, True):
				geneId_2 = session.get_defined_node_attributes("proteome", userEntity_id2, bianaInputType, True).pop()
			else:
				self.logger.warning("No {0} id for user entity {1}.".format(bianaInputType, userEntity_id2))
				continue
			
			for current_eErID in eErIDs_list:
				relationObj = relationObj_dict[current_eErID]
		
				if "psimi_name" in relationObj.get_attributes_dict():
					method_names.update([ str(x.value) for x in relationObj.get_attributes_dict()["psimi_name"] ])
				if "method_id" in relationObj.get_attributes_dict():
					method_ids.update([ x.value for x in relationObj.get_attributes_dict()["method_id"]])
				
				self.add_edge(gene_id1=geneId_1, gene_id2=geneId_2)
				self.update_edge("experimental", True, gene_id1=geneId_1, gene_id2=geneId_2)

	def iterate_genes_ScoreWise(self):
		'''
		Iterate genes that have alternative splicing and more than one transcript expressed.
		'''
		genes = self.nodes(data=True)

		if options.Options().parallelRange:
			bottom = options.Options().parallelRange - 1
			top = bottom + options.Options().step
			genes = genes[bottom:top]

		for gene,info in genes:
			allExpressedTxs = set(info["expressedTxsNormal"]) & set(info["expressedTxsTumor"])
			if len(allExpressedTxs) < 2:
				continue
			self.logger.debug("Iterating gene {0}.".format(gene))
			yield gene,info

	def iterate_switches_ScoreWise(self,tx_network,only_models=False,relevance=None,partialCreation=False,removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS 
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most 
				common) will be returned for each gene.
		"""

		counter = 1
		for gene,info in self.iterate_genes_ScoreWise():		
			if not info["isoformSwitches"]: continue

			for switchDict in info["isoformSwitches"]:
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

	def iterate_functionalSwitches_ScoreWise(self,tx_network,only_models=False,partialCreation=False,removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS 
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most 
				common) will be returned for each gene.
		"""

		self.iterate_switches_ScoreWise(tx_network,only_models,True,partialCreation,removeNoise)

	def iterate_nonFunctionalSwitches_ScoreWise(self,tx_network,only_models=False,partialCreation=False,removeNoise=True):
		"""Iterate through the isoform switches of a gene network, and
			generate a list of (gene,geneInformation,isoformSwitch).
			Only return those switches with an overlap between the CDS 
			of the transcripts and that have different features.

			only_models(bool): if True, only the first switch (the most 
				common) will be returned for each gene.
		"""

		self.iterate_switches_ScoreWise(tx_network,only_models,False,partialCreation,removeNoise)

	def createSwitch(self,switchDict,tx_network,partialCreation):
		"""Create a switch object from the switch dictionary.

			partialCreation(bool): if False, the heavy protein 
				objects are not created.
		"""
		thisSwitch = switch.IsoformSwitch(switchDict["nIso"],switchDict["tIso"],
											  switchDict["score"],switchDict["patients"],
											  switchDict["precision"],switchDict["sensitivity"])
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

			sortedSwitches = sorted(info["isoformSwitches"],key=lambda (a):len(a['patients']),reverse=True)
			for x in sortedSwitches:
				if x["nIso"] in consensus["N"] and x["tIso"] in consensus["T"]:
					bestpatible[gene] = [x["nIso"],x["tIso"],len(x["patients"])]
					break

			# no good switch could be found
			if not gene in bestpatible:
				continue

			incompatible.append(max([max( d[x]["N"] for x in d if x == bestpatible[gene][1] ),max( d[x]["T"] for x in d if x == bestpatible[gene][0] )]))

		threshold = np.percentile(incompatible,95)

		for gene,info in self.nodes(data=True):
			if not info["isoformSwitches"]: continue

			for s in info["isoformSwitches"]:
				if len(s["patients"]) <= threshold:
					s["noise"] = True
				else:
					s["noise"] = False

				if gene in bestpatible:
					if s["nIso"] == bestpatible[gene][0] and s["tIso"] == bestpatible[gene][1]:
						s["model"] = True
					else:
						s["model"] = False
				else:
					s["model"] = False

	def sampleSwitches(self,tx_network,partialCreation=True,numIterations=2000):

		genesWithSwitches = [ gene for gene,info in self.iterate_genes_ScoreWise() if info["isoformSwitches"] ]
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
			annotation = [ x for x in hallmarksDict if gene in hallmarksDict[x] ][0]
		elif [ x for x in biologicalProcessDict if gene in biologicalProcessDict[x] ]:
			annotation = [ x for x in biologicalProcessDict if gene in biologicalProcessDict[x] ][0]

		return (annotation,driverAnnotation)