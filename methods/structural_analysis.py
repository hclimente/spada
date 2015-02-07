#!/soft/devel/python-2.7/bin/python

from interface import interpro_analysis
from libs import options
from libs import utils
from methods import method

from collections import Counter
from Bio import pairwise2
import fisher
import os

import pdb

class StructuralAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self.logger.info("Structural analysis.")

		isoInfo = {}
		for line in open(options.Options().wd+"Data/TCGA/UnifiedFasta_iLoops13.fa"):
			if ">" in line:
				elements = line.strip().split("#")
				isoInfo[elements[0][1:]] = {}
				isoInfo[elements[0][1:]]["UniProt"] = elements[2]
				isoInfo[elements[0][1:]]["iLoopsFamily"] = elements[3]

		I3D_REPORT = open(options.Options().qout+"structural_analysis/I3D_analysis.tsv","w")
		I3D_REPORT.write("Gene\tSymbol\tTranscript\tUniprot\tPercent_affected_by_AS\t")
		I3D_REPORT.write("Fisher_p-value\tResidue_information\tIsoform_specific_residues\n")

		IP_REPORT = open("{0}structural_analysis/InterPro_report.tsv".format(options.Options().qout),"w")
		IP_REPORT.write("Symbol\tGene\tnTx\ttTx\tAnalysis\tFeature_accesion\t")
		IP_REPORT.write("Feature\t(Additional) repetitions\n")

		IUPRED_REPORT = open(options.Options().qout+"structural_analysis/iupred_analysis.tsv","w")
		IUPRED_REPORT.write("Gene\tTranscript\tAnalysis\tPercent_affected_by_AS\t")
		IUPRED_REPORT.write("Motifs\t#aa_Isoform_specific_ordered\t")
		IUPRED_REPORT.write("#aa_Non_isoform_specific_ordered\t")
		IUPRED_REPORT.write("#aa_Isoform_specific_disordered\t")
		IUPRED_REPORT.write("#aa_Non_isoform_specific_disordered\n")

		RELEVANCE_INFO = open(options.Options().qout+"structural_analysis/structural_summary.tsv","w")
		RELEVANCE_INFO.write("Gene\tNormalTranscript\tTumorTranscript\tiLoopsChange\t")
		RELEVANCE_INFO.write("BrokenSurface\tFunctionalChange\tDisorderChange\n")

		for gene,info,switchDict,switch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network):
			switch._iloops_change 	  = self.findDiffLoops(switch,gene,info,isoInfo)
			switch._broken_surfaces   = self.findBrokenSurfaces(switch,gene,info,I3D_REPORT)
			switch._functional_change = self.interProAnalysis(switch,gene,info,IP_REPORT)
			switch._disorder_change   = self.disorderAnalysis(switch,gene,info,IUPRED_REPORT)

			RELEVANCE_INFO.write("{0}\t{1}\t{2}\t{3}\t".format(gene,switch.nTx,switch.tTx,switch._iloops_change))
			RELEVANCE_INFO.write("{0}\t{1}\t{2}\n".format(switch._broken_surfaces,switch._functional_change,switch._disorder_change))

		I3D_REPORT.close()
		IP_REPORT.close()
		IUPRED_REPORT.close()
		RELEVANCE_INFO.close()

	def findDiffLoops(self,switch,gene,info,isoInfo):

		self.logger.debug("iLoops: looking for loop changes for gene {0}.".format(gene) )
		
		if switch.nTx in isoInfo and switch.tTx in isoInfo:
			nLoops = isoInfo[switch.nTx]["iLoopsFamily"]
			tLoops = isoInfo[switch.tTx]["iLoopsFamily"]

			if nLoops != tLoops:
				self.logger.debug("iLoops: information found for gene {0}.".format(gene) )
				return True
			else:
				return False

		return None

	def findBrokenSurfaces(self,switch,gene,info,I3D_REPORT):

		self.logger.debug("I3D: searching broken surfaces for gene {0}.".format(gene) )
		
		nIso = switch.nIsoform
		tIso = switch.tIsoform

		if nIso is None or tIso is None: return False
		elif not nIso.uniprot and not tIso.uniprot: return False
		elif not nIso.hasPdbs and not tIso.hasPdbs: return False
		
		nIsoSpecific = bool([ x for x in nIso._structure if x.isoformSpecific ])
		tIsoSpecific = bool([ x for x in tIso._structure if x.isoformSpecific ])

		if sum([nIsoSpecific,tIsoSpecific]) == 2:
			self.logger.debug("Isoform specific residues were not found exclusively in one isoform.")
			return False

		self.logger.debug("I3D: information found for gene {0}.".format(gene))

		for protein,hasIsoSpecificResidues in zip([nIso,tIso],[nIsoSpecific,tIsoSpecific]):
			if protein.hasPdbs and hasIsoSpecificResidues:
				pval,percent = self.getStatistics(protein)
				isoInfo,isoSpec = protein.report()
				I3D_REPORT.write("{0}\t{1}\t{2}\t{3}\t".format(gene,info["symbol"],protein.tx,protein.uniprot))
				I3D_REPORT.write("{0}\t{1}\t{2}\t{3}\n".format(percent,pval,isoInfo,isoSpec))
				protein.printPDBInfo()

				return True

		return False

	def getStatistics(self, protein):

		stats = { "isoSp": {"B": 0, "I": 0, "S": 0, "u": 0}, 
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

		try: percent = float(stats["isoSp"]["I"])/(stats["isoSp"]["I"]+stats["nIsoSp"]["I"])*100
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

	def interProAnalysis(self,switch,gene,info,IP_REPORT):

		self.logger.debug("InterPro: searching changes for gene {0}.".format(gene) )
		
		normalProtein = switch.nIsoform
		tumorProtein = switch.tIsoform

		for protein in [normalProtein,tumorProtein]:
			if not protein: continue

			interproFile = interpro_analysis.InterproAnalysis().launchAnalysis(protein.tx, protein.seq)
			if not interproFile: continue
			protein.getFeatures(interproFile)

		if not normalProtein or not tumorProtein: return False
		if not normalProtein._features and not tumorProtein._features: return False
		
		nIsoFeatures = Counter([ x["accession"] for x in normalProtein._features ])
		tIsoFeatures = Counter([ x["accession"] for x in tumorProtein._features ])

		nIso_uniqFeats = [ (x,nIsoFeatures[x]-tIsoFeatures.get(x,0)) for x in nIsoFeatures if nIsoFeatures[x]-tIsoFeatures.get(x,0) > 0 ]
		tIso_uniqFeats = [ (x,tIsoFeatures[x]-nIsoFeatures.get(x,0)) for x in tIsoFeatures if tIsoFeatures[x]-nIsoFeatures.get(x,0) > 0 ]

		iterationZip = zip([nIso_uniqFeats,tIso_uniqFeats],[normalProtein,tumorProtein],["Lost in tumor","Gained in tumor"])

		for uniqFeat,protein,whatsHappening in iterationZip:
			for feat,reps in uniqFeat:
				featInfo = [ x for x in protein._features if x["accession"]==feat ][0]

				IP_REPORT.write("{0}\t{1}\t{2}\t".format(info["symbol"],gene,switch.nTx))
				IP_REPORT.write("{0}\t{1}\t".format(switch.tTx,whatsHappening))
				IP_REPORT.write("{0}\t{1}\t".format(featInfo["analysis"],featInfo["accession"]))
				IP_REPORT.write("{0}\t{1}\n".format(featInfo["description"],reps))
				# IP_REPORT.write("{0}\n".format(featInfo["percentAffected"]))

				if switch._functional_change is None:
					switch._functional_change = set()

				toSave = [featInfo["analysis"],featInfo["accession"],featInfo["description"],reps]
				if not switch._functional_change:
					switch._functional_change = []
				switch._functional_change.append(toSave)

		if switch._functional_change: 
			self.logger.debug("InterPro: information found for gene {0}.".format(gene))
			return True
		else: return False

	def disorderAnalysis(self,switch,gene,info,IUPRED_REPORT):

		self.logger.debug("IUPRED: Searching disorder for gene {0}.".format(gene))
		
		normalProtein = switch.nIsoform
		tumorProtein = switch.tIsoform

		for protein in [normalProtein,tumorProtein]:
			if not protein: continue

			with open(options.Options().qout+"protein.fa","w") as FASTA:
				FASTA.write(">{0}\n{1}\n".format(protein.tx,protein.seq))

			for mode in ["short","long"]:

				counter = {	"specific" 	: { "ordered": 0.0,"disordered":0.0}, 
							"unspecific": { "ordered": 0.0,"disordered":0.0} }

				motif 	= ""
				motifs 	= []
				gap 	= ""

				proc = utils.cmdOut(options.Options().wd+"Pipeline/libs/bin/iupred/iupred",options.Options().qout+"protein.fa",mode)

				#Parse iupred output
				for line in [ x.strip().split(" ") for x in proc.stdout if "#" not in x ]:
					resNum  = int(line[0])
					residue	= line[1]
					score 	= float(line[-1])

					thisRes = protein._structure[resNum-1]
					thisRes.set_iuPredScore(score)

					res = residue.upper() if thisRes.isoformSpecific else residue.lower()

					#Extract motifs in isoform specific, disordered regions.
					#A gap of 2 can be allowed.

					if thisRes.isDisordered:
						if gap: #Check if there is a preexisting motif
							motif += gap
							gap	   = ""
						
						motif += res
					elif motif:
						if len(gap) < 2: #Max allowed gap
							gap += res
						else: 
							if len(motif) > 1:
								motifs.append(motif)
							motif 	= ""
							gap 	= ""

					#Count disordered and ordered residues, based on isoform specificity
					if thisRes.isDisordered:
						if thisRes.isoformSpecific:
							counter["specific"]["disordered"] 	+= 1
						else: 
							counter["unspecific"]["disordered"] += 1
					else:
						if thisRes.isoformSpecific:
							counter["specific"]["ordered"] 		+= 1
						else: 
							counter["unspecific"]["ordered"] 	+= 1

				if len(motif) > 1:
					motifs.append(motif)

				usefulMotifs = set()
				for motif in motifs:
					isoSpecResidues 	= float(sum(1 for c in motif if c.isupper()))
					nonIsoSpecResidues 	= float(sum(1 for c in motif if c.islower()))

					ratio = isoSpecResidues/(isoSpecResidues+nonIsoSpecResidues)	

					if ratio >= 0.2:
						usefulMotifs.add(motif)
						if switch._disorder_change is None:
							switch._disorder_change = []
						switch._disorder_change.append(motif)
			
				IUPRED_REPORT.write("{0}\t{1}\t".format(gene,info["symbol"]))
				IUPRED_REPORT.write("{0}\t{1}\t".format(protein.tx,mode))
				IUPRED_REPORT.write("{0}\t".format(",".join(usefulMotifs)))
				IUPRED_REPORT.write("{0}\t".format(int(counter["specific"]["ordered"])) )
				IUPRED_REPORT.write("{0}\t".format(int(counter["unspecific"]["ordered"])) )
				IUPRED_REPORT.write("{0}\t".format(int(counter["specific"]["disordered"])) )
				IUPRED_REPORT.write("{0}\n".format(int(counter["unspecific"]["disordered"])) )

			os.remove(options.Options().qout+"protein.fa")

		if switch._disorder_change is None or not switch._disorder_change: 
			return False
		else: 
			self.logger.debug("IUPRED: information found for gene {0}.".format(gene))
			return True
		
if __name__ == '__main__':
	pass