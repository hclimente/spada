#!/soft/devel/python-2.7/bin/python

from interface import interpro_analysis
from libs import options
from libs import utils
from methods import method

from collections import Counter
from Bio import pairwise2
import fisher

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self.logger.info("Structural analysis.")

		self.findBrokenSurfaces()
		self.interProAnalysis()
		self.disorderAnalysis()

	def findBrokenSurfaces(self):

		self.logger.info("Searching Interactome3D broken surfaces.")

		I3D_REPORT = open(options.Options().qout+"structural_analysis/I3D_analysis.tsv","w")
		I3D_REPORT.write("Gene\tTranscript\tUniprot\tPercent_affected_by_AS\t")
		I3D_REPORT.write("Fisher_p-value\tResidue_information\tIsoform_specific_residues\n")

		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene, properties in sortedNodes:

			if not properties["isoformSwitches"]: continue
			self.logger.debug("Searching structural information for gene {0}.".format(gene))
			switch = properties["isoformSwitches"][0]
			normalProtein = switch.nIsoform
			tumorProtein = switch.tIsoform
			
			if not (normalProtein and tumorProtein): continue
			elif not normalProtein.hasPdbs and not tumorProtein.hasPdbs: continue
			
			self.logger.debug("I3D information found for gene {0}.".format(gene))
		
			nIsoSpecific = bool([ x for x in normalProtein._structure if x.isoformSpecific ])
			tIsoSpecific = bool([ x for x in tumorProtein._structure if x.isoformSpecific ])

			if nIsoSpecific == tIsoSpecific:
				self.logger.debug("Isoform specific residues were not found exclusively in one isoform.")
				continue

			for protein,hasIsoSp in zip([normalProtein,tumorProtein],[nIsoSpecific,tIsoSpecific]):
				if protein.hasPdbs and hasIsoSp:
					pval,percent = self.getStatistics(protein)
					isoInfo,isoSpec = protein.report()
					I3D_REPORT.write("{0}\t{1}\t{2}\t".format(gene,protein.tx,protein.uniprot))
					I3D_REPORT.write("{0}\t{1}\t{2}\t{3}\n".format(percent,pval,isoInfo,isoSpec))
					protein.printPDBInfo()

		I3D_REPORT.close()

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

	def interProAnalysis(self):

		self.logger.info("Studying InterPro domain information.")

		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene, properties in sortedNodes:
			if not properties["isoformSwitches"]: continue

			self.logger.debug("Searching structural information for gene {0}.".format(gene))
			switch = properties["isoformSwitches"][0]
			normalProtein = switch.nIsoform
			tumorProtein = switch.tIsoform

			for protein in [normalProtein,tumorProtein]:
				if not protein: continue

				out = interpro_analysis.InterproAnalysis().launchAnalysis(protein.tx, protein.seq)
				protein.readInterpro(out)

			if not normalProtein or not tumorProtein: continue
			if not normalProtein._features and not tumorProtein._features: continue

			nIsoFeatures = Counter([ x["accession"] for x in normalProtein._features ])
			tIsoFeatures = Counter([ x["accession"] for x in tumorProtein._features ])

			uniqNIsoFeats = [ (x,nIsoFeatures[x]-tIsoFeatures.get(x,0)) for x in nIsoFeatures if nIsoFeatures[x]-tIsoFeatures.get(x,0) > 0 ]
			uniqTIsoFeats = [ (x,tIsoFeatures[x]-nIsoFeatures.get(x,0)) for x in tIsoFeatures if tIsoFeatures[x]-nIsoFeatures.get(x,0) > 0 ]

			filename = "{0}structural_analysis/IP_{1}_{2}.tsv".format(options.Options().qout,gene,properties["symbol"])
			with open(filename,"a") as INTERPRO_REPORT:
				if uniqNIsoFeats:
					INTERPRO_REPORT.write("#Normal_isoform_specific_features\n")
					for feat,reps in uniqNIsoFeats:
						INTERPRO_REPORT.write("\t#{0}\t{1}\n".format(feat,reps))

				if uniqTIsoFeats:
					INTERPRO_REPORT.write("#Tumor_isoform_specific_features\n")
					for feat,reps in uniqTIsoFeats:
						INTERPRO_REPORT.write("\t#{0}\t{1}\n".format(feat,reps))

				for protein in [normalProtein,tumorProtein]:
					for featInfo in protein._features:
		
						INTERPRO_REPORT.write("{0}\t{1}\t".format(protein.tx,featInfo["accession"]))
						INTERPRO_REPORT.write("{0}\t{1}\t".format(featInfo["analysis"],featInfo["description"]))
						INTERPRO_REPORT.write("{0}\n".format(featInfo["percentAffected"]))

	def disorderAnalysis(self):
		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		for gene, properties in sortedNodes:
			if not properties["isoformSwitches"]: continue
		
if __name__ == '__main__':
	pass