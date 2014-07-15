#!/soft/devel/python-2.7/bin/python

import pandas as pd
import numpy as np
import network
from libs import utils as ut
from libs import biana
from libs import options

class GeneNetwork(network.Network):
	"""docstring for GeneNetwork
	GeneNetwork contains a network of genes.

	Node information:
		Id - str Gene Id
		symbol - str Gene Symbol
		switch - bool One or more isoform switches detected.
		isoformSwitches - list List of detected isoform switches [[isoN, isoT], [isoN, isoT]]
		score - float Score between 0.2 and 0.5 dictacted by the number of patients.
		Driver - bool Gene described as a driver.
		RBP - bool Gene described as a RBP.
		EpiFactor - bool Gene described as epigenetic factor.
		ExpressedTranscripts - set Set with transcripts with a significant expression.

	Edge information:
		Id1 - str Gene id of interactor 1.
		Id2 - str Gene id of interactor 2.
		weight - float Weight of the itneraction.
		methods - str Methods of detection of the interaction.
		sources - str Sources of the interaction.
	"""
			
	def __init__(self):
		network.Network.__init__(self)

	def add_node(self, full_name="", gene_id="?", gene_symbol="?"):

		geneID,geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id, gene_symbol=gene_symbol)
		
		if [ x for x in self._net.nodes() if x == geneID ]:
			print("Error: node {0} exists.".format(geneID))
		else:
			self._net.add_node( 
								geneID, 
								symbol 					= geneSymbol,
								switch 					= None, 
								isoformSwitches 		= [],
								score 					= None, 
								Driver 					= False, 
								RBP 					= False, 
								EpiFactor				= False, 
								ExpressedTranscripts 	= set()
							  )

	def update_node(self, key, value, full_name = "", gene_id = ""):
		geneID, geneSymbol = self.nameFilter(full_name=full_name, gene_id=gene_id)
		
		if geneID is None:
			print("Unable to get gene id from {0} {1}".format(full_name, geneID))
			return False

		self._update_node(geneID, key, value)
		return True

	def update_edge(self, key, value, full_name1 = "", gene_id1 = "", full_name2 = "", gene_id2 = ""):
		geneIDs = []
		lst1 = [full_name1, full_name2]
		lst2 = [gene_id1, gene_id2]
		
		for i in [1,2]:
			geneID = ""
			if lst1[i]:
				geneID, geneSymbol = self.nameFilter(full_name=lst1[i])
			elif lst2[i]:
				geneID, geneSymbol = self.nameFilter(gene_id=lst2[i])

			if geneID is None or geneID is "": 
				print( "Cannot add edge {0}.".format(str([full_name1, gene_id1, full_name2, gene_id2])))
				return False

			geneIDs.append(geneID)

		self._update_edge(geneIDs[0], geneIDs[1], key, value)
		return True

	def add_edge(self, node_id1, node_id2, weight=0.0, methods="", sources=""):
		Weight 	= weight
		Methods = methods
		Sources = sources

		for w in [ self._net[x][y]["weight"] for (x,y) in self._net.edges() if (node_id1,node_id2) == (x,y) or (node_id1,node_id2) == (y,x) ]:
			Weight = min(1.0, Weight + w)

		self._net.add_edge( 
							node_id1, 
							node_id2, 
							weight=Weight, 
							methods=Methods, 
							sources=Sources 
						  )

	def nameFilter(self, full_name="", gene_id="", gene_symbol=""):
		geneSymbol 	= None
		geneID 		= None

		if full_name:
			nameComponents 	= full_name.split("|")

			if len(nameComponents) > 1:		
				geneSymbol 	= nameComponents[0]
				geneID 		= nameComponents[1]

		elif "locus" not in gene_id and "locus" not in gene_symbol:
			if gene_id:
				geneID 		= gene_id
			if gene_symbol:
				geneSymbol 	= gene_symbol

		return (geneID, geneSymbol)

	def readGeneInfo(self):
		
		for line in ut.readTable("Data/Databases/compilationTable.tsv"):
			self.add_node(full_name=line[0])
			geneID, geneSymbol = self.nameFilter(full_name=line[0])

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
				geneID, geneSymbol = self.nameFilter( gene_id=lst.pop())
				self.update_node( "ExpressedTranscripts", line[0], gene_id=geneID )

	def importCandidates(self):
		"""Import a set of genes with an isoform switch.
		"""

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
				if self._net.node[geneID]["score"] is None:
					self.update_node( "score", Score, gene_id=geneID )
				else:
					self.update_node( "score", min(1.0, self._net.node[geneID]["score"] + Score), gene_id=geneID )

	def importSpecificDrivers(self, path, otherDrivers = False):
		if not otherDrivers:
			for node in self.nodes():
				self.update_node("Driver", False, gene_id = node)
		for nameComponents in ut.readTable(path, header = False):
			geneID, geneSymbol = self.nameFilter(gene_id = nameComponents[1])
			self.update_node("Driver", True, gene_id = geneID)

	def importiLoopsInteractions(self, path, gn, nTx, tTx):
		edges = pd.DataFrame.from_csv(path + "iLoops/InteraXChanges_" + gn + "_" + nTx + "_" + tTx + "_full1.tsv", sep="\t")
		
		#Filter out "locus" genes, as they are not understandable by GUILD
		locus_filter = [ "locus" not in x for x in [ x for x in edges.Partner_gene ] ]
		edges = edges[ locus_filter ] 
		
		#Select genes for which there is only a prediction for one isoform
		dRC_filter = np.abs(edges.dRC) <= 50
		
		for aGene in edges.Partner_gene.unique():
			genenameMask = edges.Partner_gene == aGene
			geneInfo 	 = edges[ genenameMask ]
			geneScore 	 = 0.0
		
			if any(geneInfo.Annotation == "Driver"):
				geneScore += 0.5
		
			if ( genenameMask & dRC_filter ).any():
				geneScore += float( max(np.abs(geneInfo.dRC)) ) / 50
			else:
				geneScore += 0.5
		
			self.add_edge(gn, aGene, weight=geneScore)

	def importKnownInteractions(self):

		affinity_dict = { 
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
		
		affinity = set( affinity_dict.keys() )
		methods  = [("Method_id", 18),("Method_id", 696)]
		methods.extend([ ("Method_id", x) for x in affinity ])

		#seeds = [ x for x in self.nodes() if self._net.node[x]["switch"] or self._net.node[x]["driver"] ]
		#For testing purposes
		seeds = ["7157"]
		level = 2
		########

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
														#negative_attribute_restriction_list 	= [("taxid", "9606")],
														id_type 								= "geneid",
														new_user_entity_set_id					= "proteome",
													  )

		session.create_network( 
								user_entity_set_id 					= "proteome", 
								level 								= level, 
								relation_type_list 					= ["interaction"],
								relation_attribute_restriction_list = methods,
								include_relations_last_level 		= False, #Seguro?
								use_self_relations 					= False
							  )

		for (uE_id1, uE_id2) in proteome.getRelations():
			eErIDs_list 		= proteome.get_external_entity_relation_ids(uE_id1, uE_id2)
			method_names 		= set()
			method_ids 			= set()
			source_databases 	= set()
			use_method_ids 		= set()
			relationObj_dict 	= session.dbAccess.get_external_entities_dict(
																				externalEntityIdsList = eErIDs_list, 
																				attribute_list = [],
																				relation_attribute_list = ["method_id","psimi_name"], 
																				participant_attribute_list = []
																			  )
			
			
			#geneId_1 = [ x for x in session.get_user_entity(uE_id1).get_attribute(attribute_identifier="geneid") ].pop(0)
			#geneId_2 = [ x for x in session.get_user_entity(uE_id2).get_attribute(attribute_identifier="geneid") ].pop(0)
			if session.get_defined_node_attributes("proteome", uE_id1, "geneid", True):
				geneId_1 = session.get_defined_node_attributes("proteome", uE_id1, "geneid", True).pop()
			else:
				self._rejectedNodes.add(uE_id1)
				print("No gene id for " + str(uE_id1))
			
			if session.get_defined_node_attributes("proteome", uE_id2, "geneid", True):
				geneId_2 = session.get_defined_node_attributes("proteome", uE_id2, "geneid", True).pop()
			else:
				self._rejectedNodes.add(uE_id2)
				print("No gene id for " + str(uE_id2))
			
			for current_eErID in eErIDs_list:
				relationObj = relationObj_dict[current_eErID]
		
				if "psimi_name" in relationObj.get_attributes_dict():
					method_names.update([ str(x.value) for x in relationObj.get_attributes_dict()["psimi_name"] ])
				if "method_id" in relationObj.get_attributes_dict():
					method_ids.update([ x.value for x in relationObj.get_attributes_dict()["method_id"]])
				
				source_databases.add( str(session.dbAccess.get_external_database(database_id = relationObj.get_source_database() )) )
				
				nodeScore = 0.0
				edgeWeight = 0.0
				if [ x for x in method_ids if x in set([18, 696]) ]:
					nodeScore = 0.02
					edgeWeight = 1.0
				elif [ x for x in method_ids if x not in affinity ]:
					nodeScore = 0.01
					edgeWeight = 0.2
				else:
					continue

				self.add_node(full_name=geneId_1, score=nodeScore)
				self.add_node(full_name=geneId_2, score=nodeScore)
				self.add_edge(geneId_1, geneId_2, weight=edgeWeight, methods=method_names, sources=source_databases)

		nx.write_dot(self.n(), "kk.dot")

	def printGUILDInput(self, path):
		with open(path + "GUILD_nodes", "w") as NODES:
			for node, info in self.nodes(data=True):
				NODES.write("{0}\t{1}\n".format(node, info["score"]))

		with open(path + "GUILD_edges", "w") as EDGES:
			for node1, node2, info in self.edges(data=True):
				EDGES.write("{0}\t{1}\t{2}\n".format(node1, info["weight"], node2))