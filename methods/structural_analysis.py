#!/soft/devel/python-2.7/bin/python

from network import biological_entities
from network import gene_network, isoform_network
from methods import method

from collections import Counter
from Bio import pairwise2
import logging

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network):
		self._gene_network			= gn_network
		self._transcript_network	= tx_network
		
		method.Method.__init__(self)

	def run(self):
		logging.info("Structural analysis.")

		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene, properties in sortedNodes:
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

			if self._transcript_network._net[nIso]["Uniprot"] is not None:
				self.analyzeStructuralImpact(normalProtein)
			if self._transcript_network._net[tIso]["Uniprot"] is not None:
				self.analyzeStructuralImpact(tumorProtein)
			
			if normalProtein.hasPdbs and tumorProtein.hasPdbs:
				logging.debug("No Uniprot for {0} or {1}".format(nIso, tIso))
				continue
			
			logging.info("Structural information found for gene {0}.".format(gene))

			if normalProtein.hasPdbs:
				normalProtein.getAlteredRegions(tumorProtein)
				normalProtein.report()
			if tumorProtein.hasPdbs:
				tumorProtein.getAlteredRegions(normalProtein)
				tumorProtein.report()

	def analyzeStructuralImpact(self, tx, uniprot, seq, exons):

		protein = 
		noInteractions = True

		for line in utils.readTable("Data/Databases/Interactome3D/2014_01/interactions.dat"):
			interactionPdb 	= "Data/Databases/Interactome3D/2014_01/interactions/" + line[21]

			if protein.uniprot == line[0]:
				logging.debug("Relevant interaction for {0} at {1}.".format(
															protein.tx, interactionPdb))
				protein.calculateVolumes(interactionPdb, "A", self._transcript_network._net.node[protein.tx]["exonStructure"])
				noInteractions = False
			elif protein.uniprot == line[1]:
				logging.debug("Relevant interaction for {0} at {1}.".format(
															protein.tx, interactionPdb))
				protein.calculateVolumes(interactionPdb, "B")
				noInteractions = False

		if noInteractions:
			logging.debug("No relevant structures found for {0}.".format(protein.tx))
			return False

		return True

if __name__ == '__main__':
	pass