#!/soft/devel/python-2.7/bin/python

from interface import iLoops_parser
from interface import out_network
from libs import options
from libs import utils
from network import gene_network, isoform_network
from methods import method

import os
import tarfile
import time

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
		analyzedTxs = set()
		
		#Sort by score: drivers will be first to be analyzed. Then, genes with an isoform switch with 
		#decreasing number of patients
		sortedNodes = sorted(self._gene_subnetwork.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene,properties in sortedNodes:
			if not properties["isoformSwitches"]: continue

			nIso = properties["isoformSwitches"][0].nTx
			tIso = properties["isoformSwitches"][0].tTx

			status = self.getPredictedInteractions(gene, nIso, tIso, analyzedTxs)

			if status == "analysis_finished":
				self.calculateGeneInteractions(gene, nIso, tIso)
				self._gene_network.saveNetwork("geneNetwork.pkl")
				self._transcript_network.saveNetwork("txNetwork.pkl")

			if status == "analysis_finished" or status == "already_analyzed":
				out_network.interactorsInfo(self._gene_network, self._transcript_network, 
											gene, nIso, tIso )
				self.detectChanges(gene,nIso,tIso)
				analyzedTxs.add(nIso)
				analyzedTxs.add(tIso)

	def getAnalyzableCandidates(self):
		candidatesGaudi = {}
		for line in utils.readTable(options.Options().qout + "candidatesGaudi.lst"):
			index = int(line[1])
			if index >= 0 and index <= 1:
				candidatesGaudi[line[0]] = line[0]
			if index == 2:
				relative = line[2][-1].split(" ").pop()
				candidatesGaudi[line[0]] = relative

		return candidatesGaudi

	def getPredictedInteractions(self, gene, nIso, tIso, analyzed):
		parser = iLoops_parser.iLoopsParser()
		symbol = self._gene_network._net.node[gene]["symbol"]

		nIsoNeighbours = [ x for x,y in self._transcript_network._net.edges(nIso) if x != nIso]
		nIsoNeighbours.extend([ y for x,y in self._transcript_network._net.edges(nIso) if y != nIso])
		nIsoNeighbours = set(nIsoNeighbours)

		tIsoNeighbours = [ x for x,y in self._transcript_network._net.edges(tIso) if x != tIso]
		tIsoNeighbours.extend([ y for x,y in self._transcript_network._net.edges(tIso) if y != tIso])
		tIsoNeighbours = set(tIsoNeighbours)
		
		if not nIso in self._analyzable_candidates and not tIso in self._analyzable_candidates:
 			self.logger.debug("Gene {0} ({1}): lacking {2} and {3} predictions.".format(gene, symbol, nIso, tIso) )
			return "predictions_lacking"
 		
 		elif not nIso in self._analyzable_candidates:
 			self.logger.debug("Gene {0} ({1}): lacking {2} predictions.".format(gene, symbol, nIso) )
 			return "predictions_lacking"

 		elif not tIso in self._analyzable_candidates:
			self.logger.debug("Gene {0} ({1}): lacking {2} predictions.".format(gene, symbol, tIso) )
 			return "predictions_lacking"

 		elif len(nIsoNeighbours - analyzed) > 0 and len(tIsoNeighbours - analyzed) > 0:
 			self.logger.debug("Gene {0} ({1}): {2} and {3} predictions already analyzed.".format(gene, symbol, nIso, tIso) )
 			return "already_analyzed"
 		
 		self.logger.info("Gene {0} ({1}): analyzing {2} and {3} predictions.".format(gene, symbol, nIso, tIso) )
 		
 		for iso, ori in zip([nIso, tIso], ["Normal","Tumor"]):
			tarFile = "iLoops/{0}/{1}/{2}.tar.gz".format(
													options.Options().inputType, 
													options.Options().iLoopsVersion, 
													self._analyzable_candidates[iso]
												 )
			
			self.logger.info("Analyzing {0} predictions ({1} transcript).".format(iso, ori) )

			while not os.path.isfile(tarFile):
				time.sleep(900)
			
			interactions = parser.parseInteractions(self._analyzable_candidates[iso],self._expressed_transcripts)

			for partner in interactions:
				partnerSameLoops = [ x[0] for x in self._transcript_network.nodes(data=True) if x[1]["iLoopsFamily"] == self._transcript_network._net.node[partner]["iLoopsFamily"] and x[0] in self._expressed_transcripts ]
				for iso2 in partnerSameLoops:
					self._transcript_network.add_edge( iso, iso2 )
					self._transcript_network.update_edge( iso, iso2, "iLoops_prediction", True )
					self._transcript_network.update_edge( iso, iso2, "RC", interactions[partner] )

		return "analysis_finished"

	def calculateGeneInteractions(self, gene, nIso, tIso):

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

			self._gene_network.add_edge(gene_id1=gene, gene_id2=gene2)
			self._gene_network.update_edge("deltaRC",deltaRC,gene_id1=gene,gene_id2=gene2)
			self._gene_network.update_edge("iLoops_prediction",True,gene_id1=gene,gene_id2=gene2)

		for partner in [ x for x in self._transcript_network._net.neighbors(tIso) if x not in self._transcript_network._net.neighbors(nIso) ]:
			gene 	= tInfo["gene_id"]
			gene2 	= self._transcript_network._net.node[partner]["gene_id"]
			
			self._gene_network.add_edge(gene_id1=gene, gene_id2=gene2)
			self._gene_network.update_edge("deltaRC",-100,gene_id1=gene,gene_id2=gene2)
			self._gene_network.update_edge("iLoops_prediction",True,gene_id1=gene,gene_id2=gene2)

	def detectChanges(self,gene,nIso,tIso):
		nInfo = self._transcript_network._net.node[nIso]
		tInfo = self._transcript_network._net.node[tIso]

		for gene1,gene2 in self._gene_network._net.edges(gene):
			edgeInfo = self._gene_network._net.edge[gene1][gene2]
			if edgeInfo["iLoops_prediction"] and edgeInfo["experimental"]:
				if edgeInfo["deltaRC"] <= -20 and edgeInfo["score"] >= 0.6:
					self.logger.info("{0} - {1} interaction lost.".format(gene1, gene2))
				elif edgeInfo["deltaRC"] >= 20 and edgeInfo["score"] < 0.2:
					if nInfo["RC"] == 50 or tInfo["RC"] == 500:
						self.logger.info("{0} - {1} interaction gained.".format(gene1, gene2))

if __name__ == '__main__':
	pass