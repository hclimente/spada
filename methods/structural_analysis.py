#!/soft/devel/python-2.7/bin/python

from libs import utils
from biological_entities import protein
from methods import method

from collections import Counter
from Bio import pairwise2
import fisher

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self.logger.info("Structural analysis.")

		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene, properties in sortedNodes:
			if not self._gene_network._net.node[gene]["isoformSwitches"]: continue
			self.logger.debug("Searching structural information for gene {0}.".format(gene))
			switch = self._gene_network._net.node[gene]["isoformSwitches"][0][0]
			nIso = switch[0]
			tIso = switch[1]

			normalProtein = protein.Protein(nIso, 
											self._transcript_network._net.node[nIso]["Uniprot"], 
											self._transcript_network._net.node[nIso]["proteinSequence"], 
											self._transcript_network._net.node[nIso]["exonStructure"],
											self._transcript_network._net.node[nIso]["cdsCoords"],
											self._transcript_network._net.node[nIso]["strand"] )
			tumorProtein = protein.Protein( tIso, 
											self._transcript_network._net.node[tIso]["Uniprot"], 
											self._transcript_network._net.node[tIso]["proteinSequence"], 
											self._transcript_network._net.node[tIso]["exonStructure"],
											self._transcript_network._net.node[tIso]["cdsCoords"],
											self._transcript_network._net.node[tIso]["strand"] )

			if self._transcript_network._net.node[nIso]["Uniprot"] is not None:
				self.analyzeStructuralImpact(normalProtein)
			if self._transcript_network._net.node[tIso]["Uniprot"] is not None:
				self.analyzeStructuralImpact(tumorProtein)
			
			if not normalProtein.hasPdbs and not tumorProtein.hasPdbs:
				self.logger.debug("No Uniprot for {0} or {1}".format(nIso, tIso))
				continue
			
			self.logger.debug("Structural information found for gene {0}.".format(gene))
			normalProtein.getAlteredRegions(tumorProtein)
			tumorProtein.getAlteredRegions(normalProtein)

			nIsoSpecific = bool([ x for x in normalProtein._structure if x.isoformSpecific ])
			tIsoSpecific = bool([ x for x in tumorProtein._structure if x.isoformSpecific ])

			if nIsoSpecific == tIsoSpecific:
				self.logger.debug("Isoform specific residues were not found exclusively in one isoform.")
				continue

			if normalProtein.hasPdbs and nIsoSpecific:
				pval,percent = self.getStatistics(normalProtein)
				if percent >= 10:
					self.logger.info("{0}% of interaction alteration at gene {1}, normal isoform {2} (Uniprot {3}).".format(percent,gene,nIso,normalProtein.uniprot))
				self.logger.debug("{0}% of interaction alteration at gene {1}, nIso {2},  Uniprot {3}: pval {4}.".format(percent,gene,nIso,normalProtein.uniprot,pval))
				normalProtein.report()
				normalProtein.printPDBInfo()
			if tumorProtein.hasPdbs and tIsoSpecific:
				pval,percent = self.getStatistics(tumorProtein) 
				if percent >= 10:
					self.logger.info("{0}% of interaction alteration at gene {1}, tormal isoform {2} (Uniprot {3}).".format(percent,gene,tIso,normalProtein.uniprot))
				self.logger.debug("{0}% of interaction alteration at gene {1}, tIso {2},  Uniprot {3}: pval {4}.".format(percent,gene,tIso,normalProtein.uniprot,pval))
				tumorProtein.report()
				tumorProtein.printPDBInfo()

	def analyzeStructuralImpact(self, protein):

		noInteractions = True

		for line in utils.readTable("Data/Databases/Interactome3D/2014_01/interactions.dat"):
			interactionPdb 	= "Data/Databases/Interactome3D/2014_01/interactions/" + line[21]

			if protein.uniprot == line[0]:
				self.logger.debug("Relevant interaction for {0} at {1}.".format(
															protein.tx, interactionPdb))
				if protein.calculateVolumes(interactionPdb, "A", self._transcript_network._net.node[protein.tx]["exonStructure"]):
					noInteractions = False
			elif protein.uniprot == line[1]:
				self.logger.debug("Relevant interaction for {0} at {1}.".format(
															protein.tx, interactionPdb))
				if protein.calculateVolumes(interactionPdb, "B", self._transcript_network._net.node[protein.tx]["exonStructure"]):
					noInteractions = False

		if noInteractions:
			self.logger.debug("No relevant structures found for {0}, {1}.".format(protein.tx,protein.uniprot))
			return False

		return True

	def getStatistics(self, protein):

		stats = { 	"isoSp": {"B": 0, "I": 0, "S": 0, "u": 0}, 
					"nIsoSp": {"B": 0, "I": 0, "S": 0, "u": 0} }

		for residue in protein._structure: 
			if residue.isoformSpecific:
				if not residue.tag: 	 	stats["isoSp"]["u"] += 1
				elif residue.tag == "IS":	stats["isoSp"]["I"] += 1
				elif residue.tag == "NIS": 	stats["isoSp"]["S"] += 1
				elif residue.tag == "B":  	stats["isoSp"]["B"]	+= 1

			else:
				if not residue.tag: 	 	stats["nIsoSp"]["u"] += 1
				elif residue.tag == "IS":	stats["nIsoSp"]["I"] += 1
				elif residue.tag == "NIS": 	stats["nIsoSp"]["S"] += 1
				elif residue.tag == "B":  	stats["nIsoSp"]["B"] += 1

		self.logger.debug("{0}, Interacting surface:{1}\tNon-interacting surface:{2}\tBuried:{3}\tUnknown location:{4}".format(protein.tx,stats["isoSp"]["I"],stats["isoSp"]["S"],stats["isoSp"]["B"],stats["isoSp"]["u"]) )
		self.logger.debug("{0}, Interacting surface:{1}\tNon-interacting surface:{2}\tBuried:{3}\tUnknown location:{4}".format(protein.tx,stats["nIsoSp"]["I"],stats["nIsoSp"]["S"],stats["nIsoSp"]["B"],stats["nIsoSp"]["u"]) )

		try: percent = stats["isoSp"]["I"]/(stats["isoSp"]["I"]+stats["nIsoSp"]["I"])*100
		except ZeroDivisionError: percent = 0

		pval = 0

		for tag in ["S","B","I"]:
			lst = [ x for x in ["S","B","I"] if x != tag ]
			negIsoSp = 0
			negNIsoSp = 0
			
			for tag2 in lst: negIsoSp += stats["isoSp"][tag2]
			for tag2 in lst: negIsoSp += stats["nIsoSp"][tag2]

			p = fisher.pvalue(stats["isoSp"][tag], negIsoSp, stats["nIsoSp"][tag], negNIsoSp)
			self.logger.debug("{0} - {1} p-values: left:{2}\tright:{3}.".format(protein.tx,tag,p.left_tail,p.right_tail))

			if tag == "I": pval = p.right_tail

		return (pval,percent)

if __name__ == '__main__':
	pass