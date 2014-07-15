#!/soft/devel/python-2.7/bin/python

import network
from libs import utils as ut
from libs import biana
from libs import options

import pandas as pd
import numpy as np
import abc
import logging

class GeneNetwork(network.Network):
	"""docstring for GeneNetwork
	GeneNetwork contains a network of genes.

	Node information:
		Id - str Gene Id
		symbol - str Gene Symbol
		switch - bool One or more isoform switches detected.
		isoformSwitches - list List of detected isoform switches [[isoN, isoT], [isoN, isoT]]
		score - float,0.01 Score between 0.2 and 0.5 dictacted by the number of patients.
		Driver - bool Gene described as a driver.
		RBP - bool Gene described as a RBP.
		EpiFactor - bool Gene described as epigenetic factor.
		ExpressedTranscripts - set Set with transcripts with a significant expression.

	Edge information:
		Id1 - str Gene id of interactor 1.
		Id2 - str Gene id of interactor 2.
		score, 0.01 - float,None Weight of the itneraction.
		#methods - str,None Methods of detection of the interaction.
		#sources - str,None Sources of the interaction.
		iLoops_prediction - bool,None Interaction predicted by iLoops.
		experimental - bool,None Interaction found through Y2H or Lumier.

	GUILD score:

		Nodes:
			Base: 		0.01
			Maximum: 	1
			Patients:	0.2 - 0.5
			Driver: 	1
		Edges:
			Base: 		0.01
			Maximum: 	1
			iLoops:		?
			KnownPPI:	1 (Lumier, Y2H) - 0.2 (others)
		
	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self):
		network.Network.__init__(self)

	@abc.abstractmethod
	def nameFilter(self, **kwds):
		raise NotImplementedError()

	def add_node(self, full_name="", gene_id="?", gene_symbol="?"):

		logging.debug("Importing node: full name {0} id {1} symbol {2}".format(
								full_name, gene_id, gene_symbol) )
		geneID,geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id, gene_symbol=gene_symbol)
		
		if geneID in self._net.nodes():
			logging.error("Node {0} already exist.".format(geneID))
		elif geneID is None:
			logging.error("Could not retrieve name from node: full name {0} id {1} symbol {2}".format(
								full_name, gene_id, gene_symbol) )
		else:
			logging.debug("Node {0} imported.".format(geneID))
			self._net.add_node( 
								geneID, 
								symbol 					= geneSymbol,
								switch 					= None, 
								isoformSwitches 		= [],
								score 					= 0.01, 
								Driver 					= False, 
								RBP 					= False, 
								EpiFactor				= False, 
								ExpressedTranscripts 	= set()
							  )

	def update_node(self, key, value, full_name = "", gene_id = ""):
		geneID, geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id)
		finalValue = value
		
		if geneID is None:
			logging.error("Unable to get gene id from {0} {1}".format(full_name, geneID))
			return False

		if key is "score":
			finalValue = min(1.0, self._net.node[geneID]["score"] + value)

		return self._update_node(geneID, key, finalValue)

	def add_edge(self, full_name1 = "", gene_id1 = "", full_name2 = "", gene_id2 = ""):

		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2)[0]

		if (node_id1 is None or node_id1 is "") or (node_id2 is None or node_id2 is ""): 
			logging.warning( "Cannot add edge {1} - {3} ({0} - {2}).".format(
									full_name1, gene_id1, full_name2, gene_id2) )
			return False
		elif node_id1 not in self.nodes():
			logging.warning("Node {0} does not exist.".format(node_id1))
			return False
		elif node_id2 not in self.nodes():
			logging.warning("Node {0} does not exist.".format(node_id2))
			return False

		return self._add_edge( 
								node_id1, 
								node_id2, 
								score 				= 0.01, 
								iLoops_prediction 	= None,
								experimental 		= None
							 )

	def update_edge(self, key, value, full_name1 = "", gene_id1 = "", full_name2 = "", gene_id2 = ""):
		
		node_id1 = self.nameFilter(full_name=full_name1, gene_id=gene_id1)[0]
		node_id2 = self.nameFilter(full_name=full_name2, gene_id=gene_id2)[0]

		return self._update_edge(node_id1, node_id2, key, value)

	def readGeneInfo(self):
		"""Read tsv files containing characteristics of the genes:
			- compilationTable.tsv: gene annotation.
			- expressedGenes.lst: transcripts above a threshold of expression.
		"""
		
		for line in ut.readTable("Data/Databases/compilationTable.tsv"):
			self.add_node(full_name=line[0])
			geneID = self.nameFilter(full_name=line[0])[0]

			apoptosis = [line[10]]
			embryo = [line[9]]
			
			if "yes" in [line[8], line[13]]:
				self.update_node( "Driver", True, gene_id = geneID )
			if "yes" in [line[3], line[4], line[5], line[11]]:
				self.update_node( "RBP", True, gene_id = geneID )
			if "yes" in [line[6]]:
				self.update_node( "EpiFactor", True, gene_id = geneID )
		
		for line in ut.readTable(options.Options().qout + "expressedGenes.lst", header=False):
			geneId 	= ""
			lst 	= [ x for x in self.nodes() if self._net.node[x]["symbol"] == line[1] ]
			
			if lst: 
				geneID = self.nameFilter( gene_id=lst.pop())[0]
				self.update_node( "ExpressedTranscripts", line[0], gene_id=geneID )

	def importCandidates(self):
		"""Import a set of genes with an isoform switch.
		"""
		logging.debug("Retrieving calculated isoform switches.")

		samples = options.Options().replicates
		min_samples = round(samples * 0.1)
		candidates = pd.DataFrame.from_csv(options.Options().qout + "candidateList.top.tsv", sep="\t")
		candidates = candidates[ candidates.Replicates >= min_samples ]
		candidates_filt = candidates[ ["Gene", "Replicates"] ]
		
		#Remove less significant switches
		candidates_filt = candidates_filt.groupby("Gene").sum()
		candidates_filt.Replicates = 0.2 + 0.3 * (candidates_filt.Replicates - min_samples)/(samples - min_samples)
				
		for gene in candidates_filt.index.values:
			nIso = candidates.loc[candidates.Gene==gene, "Transcript_normal"].iloc[0]
			tIso = candidates.loc[candidates.Gene==gene, "Transcript_tumor"].iloc[0]
			Score = candidates_filt.loc[gene, "Replicates"]
			geneID, geneSymbol = self.nameFilter(full_name=gene)

			if geneID is not None:
				self.update_node( "switch", True, gene_id=geneID )
				self.update_node( "isoformSwitches", [nIso, tIso], gene_id=geneID )
				self.update_node( "score", Score, gene_id=geneID )

	def importSpecificDrivers(self, drivers_file, otherDrivers = False):
		logging.debug("Importing specific drivers.")
		if not otherDrivers:
			for node in [ x for x in self.nodes() if self._net.node[x]["Driver"] ]:
				self.update_node("Driver", False, gene_id = node)
				self.update_node("score", 0.01, gene_id = node)
		for nameComponents in ut.readTable(drivers_file, header = False):
			geneID, geneSymbol = self.nameFilter(gene_id = nameComponents[1])
			self.update_node("Driver", True, gene_id = geneID)
			self.update_node("score", 1, gene_id = geneID)

	def importKnownInteractions(self):

		logging.debug("Importing interactions from BIANA.")

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
		
		methods  = [("Method_id", 18),("Method_id", 696)]
		methods.extend([ ("Method_id", x) for x in complementation_methods if x not in affinity_methods ])

		seeds = [ x for x in self.nodes() if self._net.node[x]["switch"] or self._net.node[x]["Driver"] ]

		bianaInputType = "geneid"

		if options.Options().inputType == "ensembl":
			bianaInputType = "ensembl" #es asi el identificador?

		session = biana.create_new_session(
										sessionID="SmartAS", 
										dbname="BIANA_MARCH_2013", 
										dbhost="ben-yehuda",
										dbuser="biana_user", 
										dbpassword="biana_password",
										unification_protocol="uniprot_geneID_seqtax"
									)

		proteome = session.create_new_user_entity_set(
														identifier_description_list 			= seeds,
														attribute_restriction_list 				= [("taxid", "9606")],
														id_type 								= bianaInputType,
														new_user_entity_set_id					= "proteome",
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
			logging.debug( "Interaction {0}/{1}".format(c, len(proteome.getRelations())) )
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
				logging.warning("No {0} id for user entity {1}.".format(bianaInputType, userEntity_id1))
				continue
			
			if session.get_defined_node_attributes("proteome", userEntity_id2, bianaInputType, True):
				geneId_2 = session.get_defined_node_attributes("proteome", userEntity_id2, bianaInputType, True).pop()
			else:
				logging.warning("No {0} id for user entity {1}.".format(bianaInputType, userEntity_id2))
				continue
			
			for current_eErID in eErIDs_list:
				relationObj = relationObj_dict[current_eErID]
		
				if "psimi_name" in relationObj.get_attributes_dict():
					method_names.update([ str(x.value) for x in relationObj.get_attributes_dict()["psimi_name"] ])
				if "method_id" in relationObj.get_attributes_dict():
					method_ids.update([ x.value for x in relationObj.get_attributes_dict()["method_id"]])
				
				self.add_edge(gene_id1=geneId_1, gene_id2=geneId_2)
				self.update_edge("experimental", True, gene_id1=geneId_1, gene_id2=geneId_2)

				if [ x for x in method_ids if x in set([18, 696]) ]:
					self.update_edge("score", 1.0, gene_id1=geneId_1, gene_id2=geneId_2)
				elif [ x for x in method_ids if x not in affinity_methods ]:
					self.update_edge("score", 0.2, gene_id1=geneId_1, gene_id2=geneId_2)