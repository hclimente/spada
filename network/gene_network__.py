#!/soft/devel/python-2.7/bin/python

import pandas as pd
import numpy as np
import network
from libs import utils as ut
from libs import biana

class GeneNetwork(network.Network):
	def __init__(self):
		network.Network.__init__(self)

	def add_node(self, full_name="", normalIso="", tumorIso="", score=0.0, gene_id="?", gene_symbol="?", switch=False, driver=False):
		debug = False

		geneSymbol 	= gene_symbol
		geneID 		= gene_id
		Score 		= score
		
		if full_name:
			nameComponents 	= full_name.split("|")

			if len(nameComponents) == 1:
				if debug: print(full_name + " not added. Unable to get gene id and symbol.")
				return
		
			geneSymbol 	= nameComponents[0]
			geneID 		= nameComponents[1]

		elif "locus" in gene_id:
			if debug: print( "Unknown gene " + gene_id + " not added.")
			return
		elif "locus" in gene_symbol:
			if debug: print( "Unknown gene " + gene_symbol + " not added.")
			return

		for s in [ self._net.node[x]["score"] for x in self._net.nodes() if x == geneID ]:
			Score = min(1.0, s + Score)

		self._net.add_node( geneID, normal_iso=normalIso, tumor_iso=tumorIso, score=Score, symbol=geneSymbol, switch=switch, driver=driver )

	def add_edge(self, node_id1, node_id2, weight=0.0, methods="", sources=""):
		Weight 	= weight
		Methods = methods
		Sources = sources

		for w in [ self._net[x][y]["weight"] for (x,y) in self._net.edges() if (node_id1,node_id2) == (x,y) or (node_id1,node_id2) == (y,x) ]:
			Weight = min(1.0, Weight + w)

		self._net.add_edge( node_id1, node_id2, weight=Weight, methods=Methods, sources=Sources )

	def importCandidates(self, path, samples):
		"""Import a set of genes with an isoform switch.
		Keyword arguments:
    	path -- path containing the candidateList file.
    	samples -- total number of samples, for filtering purposes.
		"""
		min_samples = round(samples * 0.1)
		candidates = pd.DataFrame.from_csv(path + "candidateList.top.tsv", sep="\t")
		candidates = candidates[ candidates.Replicates >= min_samples ]
		candidates_filt = candidates[ ["Gene", "Replicates"] ]
		
		#Remove less significant switches
		candidates_filt = candidates_filt.groupby("Gene").sum()
		candidates_filt.Replicates = 0.2 + 0.3 * (candidates_filt.Replicates - min_samples)/(samples - min_samples)
				
		for gene in candidates_filt.index.values:
			nIso = candidates.loc[candidates.Gene==gene, "Transcript_normal"].iloc[0]
			tIso = candidates.loc[candidates.Gene==gene, "Transcript_tumor"].iloc[0]
			Score = candidates_filt.loc[gene, "Replicates"]
			self.add_node( full_name=gene, score=Score, normalIso=nIso, tumorIso=tIso, switch=True )

	def importDrivers(self, path):
		for nameComponents in ut.readTable(path):
			self.add_node( gene_symbol=nameComponents[0], gene_id=nameComponents[1], score=1.0, driver=True )

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