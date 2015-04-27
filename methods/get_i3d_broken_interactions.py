from libs import options
from libs import utils
from methods import method

import fisher

class GetI3DBrokenInteractions(method.Method):
	def __init__(self,gn_network,tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

		self.hallmarks = self.readGeneset("h.all.v5.0.entrez.gmt")
		self.biologicalProcess = self.readGeneset("c5.bp.v4.0.entrez.gmt")

	def run(self):

		self.logger.info("Searching I3D broken surfaces.")

		self.OUT = open("{0}i3d/i3d_broken.tsv".format(options.Options().qout),'w')
		self.OUT.write("Gene\tSymbol\tnTx\ttTx\tCancer\tUniprot\tTag\tWhatsHappening\t")
		self.OUT.write("Percent\tp\tOR\tSequenceCover\tIsoformStructure\tIsoformSpecific\n")
		
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,partialCreation=False,removeNoise=True,only_models=True):
			self.findBrokenSurfaces(thisSwitch,gene,info)

		self.OUT.close()

	def clean(self):
		utils.cmd("rm","-r","{0}i3d".format(options.Options().qout))
		utils.cmd("mkdir","{0}i3d".format(options.Options().qout))
	
	def findBrokenSurfaces(self,thisSwitch,gene,info):

		self.logger.debug("I3D: searching broken surfaces for gene {0}.".format(gene) )
		
		nIso = thisSwitch.nIsoform
		tIso = thisSwitch.tIsoform

		if nIso is None or tIso is None: return False
		elif not nIso.uniprot and not tIso.uniprot: return False
		elif not nIso.hasPdbs and not tIso.hasPdbs: return False
		
		nIsoSpecific = bool([ x for x in nIso._structure if x.isoformSpecific ])
		tIsoSpecific = bool([ x for x in tIso._structure if x.isoformSpecific ])

		if sum([nIsoSpecific,tIsoSpecific]) == 2:
			self.logger.debug("Isoform specific residues were not found exclusively in one isoform.")
			return False

		self.logger.debug("I3D: information found for gene {0}.".format(gene))

		tag = "Nothing"
		if info["Driver"]:
			tag = "Driver"
		elif [ x for x in self._gene_network._net.neighbors(gene) if self._gene_network._net.node[x]["Driver"] ]:
			tag = "d1"
		elif [ x for x in self.hallmarks if gene in self.hallmarks[x] ]:
			hallmark = [ x for x in self.hallmarks if gene in self.hallmarks[x] ][0]
			tag = hallmark
		elif [ x for x in self.biologicalProcess if gene in self.biologicalProcess[x] ]:
			bp = [ x for x in self.biologicalProcess if gene in self.biologicalProcess[x] ][0]
			tag = bp

		for protein,hasIsoSpecificResidues,what in zip([nIso,tIso],[nIsoSpecific,tIsoSpecific],["nIso_interaction","tIso_interaction"]):
			if protein.hasPdbs and hasIsoSpecificResidues:
				pval,OR,percent = self.getStatistics(protein)
				isoInfo,isoSpec = protein.report()
				seqCover = float(len(isoInfo)-isoInfo.count("*"))/len(isoInfo)
				self.OUT.write("{0}\t{1}\t{2}\t".format(gene,info["symbol"],nIso.tx))
				self.OUT.write("{0}\t{1}\t{2}\t".format(tIso.tx,options.Options().tag,protein.uniprot))
				self.OUT.write("{0}\t{1}\t{2}\t".format(tag,what,percent))
				self.OUT.write("{0}\t{1}\t{2}\n".format(pval,OR,seqCover))
				self.OUT.write("{0}\t{1}\n".format(isoInfo,isoSpec))
				protein.printPDBInfo()

				return True

		return False

	def getStatistics(self,protein):

		stats = { "isoSp": {"B": 0, "I": 0, "S": 0, "u": 0}, 
				  "nIsoSp": {"B": 0, "I": 0, "S": 0, "u": 0} }

		for residue in protein._structure: 
			if residue.isoformSpecific:
				if residue.tag is None: 	stats["isoSp"]["u"] += 1
				elif residue.tag == "IS":	stats["isoSp"]["I"] += 1
				elif residue.tag == "NIS": 	stats["isoSp"]["S"] += 1
				elif residue.tag == "B":  	stats["isoSp"]["B"]	+= 1

			else:
				if residue.tag is None: 	stats["nIsoSp"]["u"] += 1
				elif residue.tag == "IS":	stats["nIsoSp"]["I"] += 1
				elif residue.tag == "NIS": 	stats["nIsoSp"]["S"] += 1
				elif residue.tag == "B":  	stats["nIsoSp"]["B"] += 1

		self.logger.debug("{0}, Interacting surface:{1}\tNon-interacting surface:{2}\tBuried:{3}\tUnknown location:{4}".format(protein.tx,stats["isoSp"]["I"],stats["isoSp"]["S"],stats["isoSp"]["B"],stats["isoSp"]["u"]) )
		self.logger.debug("{0}, Interacting surface:{1}\tNon-interacting surface:{2}\tBuried:{3}\tUnknown location:{4}".format(protein.tx,stats["nIsoSp"]["I"],stats["nIsoSp"]["S"],stats["nIsoSp"]["B"],stats["nIsoSp"]["u"]) )

		try: percent = float(stats["isoSp"]["I"])/(stats["isoSp"]["I"]+stats["nIsoSp"]["I"])*100
		except ZeroDivisionError: percent = 0

		pval = "NA"
		OR = "NA"

		for tag in ["S","B","I"]:
			lst = [ x for x in ["S","B","I"] if x != tag ]
			negIsoSp = 0
			negNIsoSp = 0
			
			for tag2 in lst: negIsoSp += stats["isoSp"][tag2]
			for tag2 in lst: negIsoSp += stats["nIsoSp"][tag2]

			p = fisher.pvalue(stats["isoSp"][tag],negIsoSp,stats["nIsoSp"][tag],negNIsoSp)
			try:
				oddsRatio = stats["isoSp"][tag]*stats["nIsoSp"][tag]/(negIsoSp*negNIsoSp)
			except ZeroDivisionError: 
				oddsRatio = "NA"
			self.logger.debug("{0} - {1} p-values: left:{2}\tright:{3}. OR: {4}".format(protein.tx,tag,p.left_tail,p.right_tail,oddsRatio))

			if tag == "I": 
				pval = p.right_tail
				OR = oddsRatio

		return (pval,OR,percent)

	def readGeneset(self,sSetFile):

		geneSets = {}
		geneSetFile = "{0}Data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = line[2:]
			geneSets[geneSet] = genes

		return geneSets