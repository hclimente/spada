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
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)
		self._expressed_transcripts = set()

		for txSet in [ self._gene_network._net.node[x]["ExpressedTranscripts"] for x in self._gene_network.nodes() ]:
			for tx in txSet:
				self._expressed_transcripts.add(tx)
		self._analyzable_candidates	= self.getAnalyzableCandidates()

	def run(self):
		self.logger.info("Examining iLoops results.")
		
		#Sort by score: drivers will be first to be analyzed. Then, genes with an isoform switch with 
		#decreasing patients

		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene,properties in sortedNodes:
			if not self._gene_network._net.node[gene]["isoformSwitches"]: continue
			switch = self._gene_network._net.node[gene]["isoformSwitches"][0][0]
			nIso = switch[0]
			tIso = switch[1]

			if self.getPredictedInteractions(gene, nIso, tIso):
				self.calculateGeneInteractions(gene, nIso, tIso)
				out_network.interactorsInfo(self._gene_network, self._transcript_network, 
											gene, nIso, tIso )

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
		symbol = self._gene_network._net.node[gene]["symbol"]
		
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
				self._transcript_network.add_edge( iso, partner )
				self._transcript_network.update_edge( iso, partner, "iLoops_prediction", True )
				self._transcript_network.update_edge( iso, partner, "RC", interactions[partner] )
			
			self.logger.info("Predictions for {0} analyzed.".format(iso) )
			self._gene_network.saveNetwork("geneNetwork_2.pkl")
			self._transcript_network.saveNetwork("txNetwork_2.pkl")

		return True

	def calculateGeneInteractions(self, gene, nIso, tIso):

		for partner in self._transcript_network._net.neighbors(nIso):
			gene 	= self._transcript_network._net.node[nIso]["gene_id"]
			nRC		= self._transcript_network._net.edge[nIso][partner]["RC"]
			deltaRC = 100
			tRC 	= ""
			if partner not in self._transcript_network._net.neighbors(tIso):
				tRC = self._transcript_network._net.edge[tIso][partner]["RC"]
				deltaRC = nRC - tRC

			for isoSameLoops in [ x for x in self._transcript_network.nodes() if self._transcript_network._net.node[x]["iLoopsFamily"] == self._transcript_network._net.node[partner]["iLoopsFamily"] ]:
				self._gene_network.add_edge(
											 gene_id1=gene, 
											 gene_id2=self._transcript_network._net.node[isoSameLoops]["gene_id"]
										   )
				self._gene_network.update_edge( "score", deltaRC, gene_id1=gene, gene_id2=self._transcript_network._net.node[isoSameLoops]["gene_id"])

		for partner in [ x for x in self._transcript_network._net.neighbors(tIso) if x not in self._transcript_network._net.neighbors(nIso) ]:
			gene 	= self._transcript_network._net.node[tIso]["gene_id"]
			tRC		= self._transcript_network._net.edge[tIso][partner]["RC"]
			
			for isoSameLoops in [ x for x in self._transcript_network.nodes() if self._transcript_network._net.node[x]["iLoopsFamily"] == self._transcript_network._net.node[partner]["iLoopsFamily"] ]:
				self._gene_network.add_edge(
											 gene_id1=gene, 
											 gene_id2=self._transcript_network._net.node[isoSameLoops]["gene_id"]
										   )
				self._gene_network.update_edge( "score", -100, gene_id1=gene, gene_id2=self._transcript_network._net.node[isoSameLoops]["gene_id"])

		self._gene_network.saveNetwork("geneNetwork_3.pkl")
		self._transcript_network.saveNetwork("txNetwork_3.pkl")

if __name__ == '__main__':
	pass