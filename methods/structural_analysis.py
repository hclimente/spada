#!/soft/devel/python-2.7/bin/python

from network import biological_entities
from network import gene_network, isoform_network
from methods import method

from collections import Counter
from Bio import pairwise2

import pdb

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)

	def run(self):
		self.logger.info("Structural analysis.")

		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene, properties in sortedNodes:
			if not self._gene_network._net.node[gene]["isoformSwitches"]: continue
			self.logger.debug("Searching structural information for gene {0}.".format(gene))
			switch = self._gene_network._net.node[gene]["isoformSwitches"][0]
			nIso = switch[0]
			tIso = switch[1]

			normalProtein = biological_entities.Protein( 
											nIso, 
											self._transcript_network._net.node[nIso]["Uniprot"], 
											self._transcript_network._net.node[nIso]["proteinSequence"], 
											self._transcript_network._net.node[nIso]["exonStructure"]
														)
			tumorProtein = biological_entities.Protein( 
											tIso, 
											self._transcript_network._net.node[tIso]["Uniprot"], 
											self._transcript_network._net.node[tIso]["proteinSequence"], 
											self._transcript_network._net.node[tIso]["exonStructure"]
													  )

			if self._transcript_network._net.node[nIso]["Uniprot"] is not None:
				self.analyzeStructuralImpact(normalProtein)
			if self._transcript_network._net.node[tIso]["Uniprot"] is not None:
				self.analyzeStructuralImpact(tumorProtein)
			
			if not normalProtein.hasPdbs and not tumorProtein.hasPdbs:
				self.logger.debug("No Uniprot for {0} or {1}".format(nIso, tIso))
				continue
			
			self.logger.info("Structural information found for gene {0}.".format(gene))
			pdb.set_trace()

			if normalProtein.hasPdbs:
				normalProtein.getAlteredRegions(tumorProtein)
				normalProtein.report()
			if tumorProtein.hasPdbs:
				tumorProtein.getAlteredRegions(normalProtein)
				tumorProtein.report()

	def analyzeStructuralImpact(self, protein):

		noInteractions = True

		for line in utils.readTable("Data/Databases/Interactome3D/2014_01/interactions.dat"):
			interactionPdb 	= "Data/Databases/Interactome3D/2014_01/interactions/" + line[21]

			if protein.uniprot == line[0]:
				self.logger.debug("Relevant interaction for {0} at {1}.".format(
															protein.tx, interactionPdb))
				protein.calculateVolumes(interactionPdb, "A", self._transcript_network._net.node[protein.tx]["exonStructure"])
				noInteractions = False
			elif protein.uniprot == line[1]:
				self.logger.debug("Relevant interaction for {0} at {1}.".format(
															protein.tx, interactionPdb))
				protein.calculateVolumes(interactionPdb, "B")
				noInteractions = False

		if noInteractions:
			self.logger.debug("No relevant structures found for {0}.".format(protein.tx))
			return False

		return True

if __name__ == '__main__':
	pass