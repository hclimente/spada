#!/soft/devel/python-2.7/bin/python

from libs import utils
from network import gene_network, isoform_network
from methods import method
from libs import options
from interface import iLoops_parser
from interface import out_network

import os
import tarfile
import time

import pdb

class AnalyzeInteractions(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)
		self._expressed_transcripts = set()

		for txSet in [ self._gene_network._net.node[x]["ExpressedTranscripts"] for x in self._gene_network.nodes() ]:
			for tx in txSet:
				self._expressed_transcripts.add(tx)
		self._analyzable_candidates	= self.getAnalyzableCandidates()

	def run(self):
		self.logger.info("Examining iLoops results.")
		
		#Sort by score: drivers will be first to be analyzed. Then, genes with an isoform switch with 
		#decreasing number of patients

		sortedNodes = sorted(self._gene_subnetwork.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene,properties in sortedNodes:
			if not properties["isoformSwitches"]: continue
			for switch,score,patients in properties["isoformSwitches"]:
				nIso = switch[0]
				tIso = switch[1]

				if self.getPredictedInteractions(gene, nIso, tIso):
					self.calculateGeneInteractions(gene, nIso, tIso)
					out_network.interactorsInfo(self._gene_network, self._transcript_network, 
												gene, nIso, tIso )
					self.detectChanges(gene,nIso,tIso)

					self._gene_subnetwork.saveNetwork("geneSubnetwork.pkl")
					self._transcript_network.saveNetwork("txNetwork.pkl")
					break

	def getAnalyzableCandidates(self):
		candidatesGaudi = {}
		for line in utils.readTable(options.Options().qout + "candidatesGaudi.lst"):
			index = int(line[1])
			if index >= 0 and index <= 1:
				candidatesGaudi[line[0]] = line[0] + ".tar.gz"
			if index == 2:
				relative = line[2][-1].split(" ").pop()
				candidatesGaudi[line[0]] = relative + ".tar.gz"

		return candidatesGaudi

	def getPredictedInteractions(self, gene, nIso, tIso):
		parser = iLoops_parser.iLoopsParser()
		symbol = self._gene_subnetwork._net.node[gene]["symbol"]
		
		if not nIso in self._analyzable_candidates and not tIso in self._analyzable_candidates:
 			self.logger.debug("Gene {0} ({1}): lacking {2} and {3} predictions.".format(gene, symbol, nIso, tIso) )
			return False
 		
 		elif not nIso in self._analyzable_candidates:
 			self.logger.debug("Gene {0} ({1}): lacking {2} predictions.".format(gene, symbol, nIso) )
 			return False

 		elif not tIso in self._analyzable_candidates:
			self.logger.debug("Gene {0} ({1}): lacking {2} predictions.".format(gene, symbol, tIso) )
 			return False

 		elif self._transcript_network._net.edges(nIso) and self._transcript_network._net.edges(tIso):
 			self.logger.debug("Gene {0} ({1}): {2} and {3} predictions already analyzed.".format(gene, symbol, nIso, tIso) )
 			return True
 		
 		self.logger.info("Gene {0} ({1}): analyzing {2} and {3} predictions.".format(gene, symbol, nIso, tIso) )
 		
 		for iso, ori in zip([nIso, tIso], ["Normal","Tumor"]):
			tarFile = "iLoops/{0}/{1}/{2}".format(
													options.Options().inputType, 
													options.Options().iLoopsVersion, 
													self._analyzable_candidates[iso]
												 )
			
			self.logger.info("Analyzing {0} predictions ({1} transcript).".format(iso, ori) )

			while not os.path.isfile(tarFile):
				time.sleep(900)
			
			tar = tarfile.open(tarFile)
			xmlFile = tar.getmembers()[0].name
			tar.extract(xmlFile, path="iLoops/{0}/{1}/".format(options.Options().inputType, 
																  options.Options().iLoopsVersion ))
			tar.close()

			xmlFile = "iLoops/{0}/{1}/{2}".format(options.Options().inputType, 
												  options.Options().iLoopsVersion,
												  xmlFile )

			interactions = parser.parseInteractions(
									thisCandidate				  = iso,
									xmlOutput					  = xmlFile,
									expressedIsoforms			  = self._expressed_transcripts,
									output_proteins               = False, 
									output_alignments             = False,
									output_domain_mappings        = False,
									output_protein_features       = False,
									output_domain_assignations    = False,
									output_interactions           = True,
									output_interaction_signatures = False,
									output_RF_results             = True,
									output_RF_precisions          = True
													)
			os.remove(xmlFile)

			for partner in interactions:
				partnerSameLoops = [ x[0] for x in self._transcript_network.nodes(data=True) if x[1]["iLoopsFamily"] == self._transcript_network._net.node[partner]["iLoopsFamily"] ]
				for iso2 in partnerSameLoops:
					self._transcript_network.add_edge( iso, iso2 )
					self._transcript_network.update_edge( iso, iso2, "iLoops_prediction", True )
					self._transcript_network.update_edge( iso, iso2, "RC", interactions[partner] )

		return True

	def calculateGeneInteractions(self, gene, nIso, tIso):

		if [ x for x in self._gene_subnetwork._net.edges(gene,data=True) if x[2]["iLoops_prediction"] ]:
 			self.logger.debug("Gene {0}: {1} and {2} predictions already summarized.".format(gene, nIso, tIso) )
 			return

		self.logger.info("Extrapolating predictions for {0} and {1} to gene level.".format(nIso,tIso))
		nInfo = self._transcript_network._net.node[nIso]
		tInfo = self._transcript_network._net.node[tIso]

		for partner in self._transcript_network._net.neighbors(nIso):
			gene 	= nInfo["gene_id"]
			gene2 	= self._transcript_network._net.node[partner]["gene_id"]
			deltaRC = 100

			if partner in self._transcript_network._net.neighbors(tIso):
				nRC 	= self._transcript_network._net.edge[nIso][partner]["RC"]
				tRC 	= self._transcript_network._net.edge[tIso][partner]["RC"]
				deltaRC = nRC - tRC

			self._gene_subnetwork.add_edge(gene_id1=gene, gene_id2=gene2)
			self._gene_subnetwork.update_edge("deltaRC",deltaRC,gene_id1=gene,gene_id2=gene2)
			self._gene_subnetwork.update_edge("iLoops_prediction",True,gene_id1=gene,gene_id2=gene2)

		for partner in [ x for x in self._transcript_network._net.neighbors(tIso) if x not in self._transcript_network._net.neighbors(nIso) ]:
			gene 	= tInfo["gene_id"]
			gene2 	= self._transcript_network._net.node[partner]["gene_id"]
			
			self._gene_subnetwork.add_edge(gene_id1=gene, gene_id2=gene2)
			self._gene_subnetwork.update_edge("deltaRC",-100,gene_id1=gene,gene_id2=gene2)
			self._gene_subnetwork.update_edge("iLoops_prediction",True,gene_id1=gene,gene_id2=gene2)

	def detectChanges(self,gene,nIso,tIso):
		nInfo = self._transcript_network._net.node[nIso]
		tInfo = self._transcript_network._net.node[tIso]

		for gene1,gene2 in self._gene_subnetwork._net.edges(gene):
			edgeInfo = self._gene_subnetwork._net.edge[gene1][gene2]
			if edgeInfo["iLoops_prediction"] and edgeInfo["experimental"]:
				if edgeInfo["deltaRC"] <= -20 and edgeInfo["score"] >= 0.6:
					self.logger.info("{0} - {1} interaction lost.".format(gene1, gene2))
				elif edgeInfo["deltaRC"] >= 20 and edgeInfo["score"] < 0.2:
					if nInfo["RC"] == 50 or tInfo["RC"] == 500:
						self.logger.info("{0} - {1} interaction gained.".format(gene1, gene2))

if __name__ == '__main__':
	pass