#!/soft/devel/python-2.7/bin/python

from libs import options
from libs import utils
from methods import method

from itertools import groupby
from operator import itemgetter
import pandas as pd

import pdb

class ResultSummary(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self.logger.info("Summarizing results.")

		self.generalSwitchInfo()
		#self.structuralInfo()
		#self.functionalChange()
		#self.printRelevantGene()

	def generalSwitchInfo(self):
		#self.proteinOverview()
		self.switchAndExonOverview()

	def structuralInfo(self):
		self.loopChange()
		self.disorderChange()
		self.functionalChange()

	def neighborhoodInfo(self):
		self.topNeighborhoods()

	def proteinOverview(self):
		sortedNodes = sorted(self._gene_network.nodes(data=True), 
							 key=lambda (a, dct): dct['score'], reverse=True)
		txDict = self._transcript_network.nodes(data=True)

		centrality = []
		nIsoLength = []
		tIsoLength = []
		
		for gene,info in sortedNodes:
			if not info["isoformSwitches"]: continue
			#elif abs(info["diffExpression_logFC"]) > 0.5 or info["diffExpression_p"] < 0.05: continue
			centrality.append(self._gene_network._net.degree(gene))
			for switchDict in info["isoformSwitches"]:
				if [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["nIso"] ]!=[None]:
					seq = [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["nIso"] ][0]
					nIsoLength.append(len(seq))
				if [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["tIso"] ]!=[None]:
					seq = [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["tIso"] ][0]
					tIsoLength.append(len(seq))

		with open("{0}result_summary/protein_centrality.tsv".format(options.Options().qout), "w" ) as F:
			for degree in centrality:
				F.write("{0}\n".format(degree))

		with open("{0}result_summary/nIso_length.tsv".format(options.Options().qout), "w" ) as F:
			for length in nIsoLength:
				F.write("{0}\n".format(length))
		
		with open("{0}result_summary/tIso_length.tsv".format(options.Options().qout), "w" ) as F:
			for length in tIsoLength:
				F.write("{0}\n".format(length))

	def switchAndExonOverview(self):
		stats = {}

		stats["hasCds"] = { "both": 0, "onlyN": 0,"onlyT": 0,"none": 0 }
		stats["hasCdsChange"] = { "y": 0, "n": 0 }
		stats["hasUtrChange"] = { "y": 0, "n": 0 }
		stats["isDriver"] = { "y": 0, "n": 0 }
		stats["isRelevant"] = { "y": 0, "n": 0 }

		total = 0
		
		for gene,info,switchDict,switch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network):
			nIso = switch.nTranscript
			tIso = switch.tTranscript

			if nIso.cds and tIso.cds: stats["hasCds"]["both"]  += 1
			elif nIso.cds: 			  stats["hasCds"]["onlyN"] += 1
			elif tIso.cds: 			  stats["hasCds"]["onlyT"] += 1
			else: 					  stats["hasCds"]["none"]  += 1

			if switch.cds_diff: stats["hasCdsChange"]["y"] += 1
			else: 				stats["hasCdsChange"]["n"] += 1

			if switch.utr_diff: stats["hasUtrChange"]["y"] += 1
			else: 				stats["hasUtrChange"]["n"] += 1

			if info["Driver"]:  stats["isDriver"]["y"] +=1
			else: 				stats["isDriver"]["n"] +=1

			if switch.is_relevant:  stats["isRelevant"]["y"] +=1
			else: 					stats["isRelevant"]["n"] +=1

			total += 1

			nSpecificCds = [ x for x in switch._normal_transcript.cds if x not in switch._tumor_transcript.cds ]
			nSpecificUtr = [ x for x in switch._normal_transcript.utr if x not in switch._tumor_transcript.utr ]
			tSpecificCds = [ x for x in switch._tumor_transcript.cds if x not in switch._normal_transcript.cds ]
			tSpecificUtr = [ x for x in switch._tumor_transcript.utr if x not in switch._normal_transcript.utr ]
			
			for specificCds,specificUtr,cds in zip([nSpecificCds,tSpecificCds],[nSpecificUtr,tSpecificUtr],[switch._normal_transcript.cds,switch._tumor_transcript.cds]):
			
				specific = specificUtr
				specific.extend(specificCds)
				specific = sorted(specific)
				
				setSpecCds = set(specificCds)
				setSpecUtr = set(specificUtr)
				setSpec = set(specific)

				exons = [ map(itemgetter(1),g) for k,g in groupby(enumerate(specific), lambda (i,x): i-x) ]
				
				for exon in exons:
					exonInfo = {}
					exonInfo["switch"] = "{0}_{1}_{2}".format(gene,switchDict["nIso"],switchDict["tIso"])
					exonInfo["length"] = len(exon)

					if setSpec & setSpecCds and setSpec & setSpecUtr:
						exonicCds = sorted(setSpec & setSpecCds)
						exonInfo["role"] = "CDS-UTR"
						exonInfo["keepORF"] = len(exonicCds)%3
						medianPos = exonicCds[len(exonicCds)/2]
						pos = [ i for i,x in enumerate(cds) if x==medianPos ][0]	 
						exonInfo["position"] = float(pos)/len(cds)
					elif setSpec & setSpecCds:
						exonInfo["role"] = "CDS"
						exonInfo["keepORF"] = len(exon)%3
						medianPos = exon[len(exon)/2]
						pos = [ i for i,x in enumerate(cds) if x==medianPos ][0]
						exonInfo["position"] = float(pos)/len(cds)
					elif setSpec & setSpecUtr:
						exonInfo["role"] = "UTR"
						exonInfo["keepORF"] = None
						exonInfo["position"] = None
					else:
						exonInfo["role"] = "wut"
						exonInfo["keepORF"] = "wut"
						exonInfo["position"] = "wut"

		with open("{0}result_summary/switches.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tAnalysis\tBoth\tOnly_nIso\tOnly_tIso\tNone\tTotal\n")
			F.write("{0}\tCDS_study\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(stats["hasCds"]["both"],stats["hasCds"]["onlyN"]) )
			F.write("{0}\t{1}\t".format(stats["hasCds"]["onlyT"],stats["hasCds"]["none"]) )
			F.write("{0}\n".format(total))

			F.write("Cancer\tAnalysis\tYes\tNo\tTotal\n")
			F.write("{0}\tCDS_change\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(stats["hasCdsChange"]["y"],stats["hasCdsChange"]["n"]) )
			F.write("{0}\n".format(total))

			F.write("{0}\tUTR_change\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(stats["hasUtrChange"]["y"],stats["hasUtrChange"]["n"]) )
			F.write("{0}\n".format(total))

			F.write("{0}\tDriver_affection\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(stats["isDriver"]["y"],stats["isDriver"]["n"]) )
			F.write("{0}\n".format(total))

			F.write("{0}\tRelevant\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(stats["isRelevant"]["y"],stats["isRelevant"]["n"]) )
			F.write("{0}\n".format(total))

		with open("{0}result_summary/exons.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tSwitch\tLength\tType\tKeepOrf\tPosition\n")
			for length in exonInfo:
				F.write("{0}\t{1}\t".format(options.Options().tag,exonInfo["switch"]))
				F.write("{0}\t{1}\t".format(exonInfo["length"],exonInfo["role"]))
				F.write("{0}\t{1}\n".format(exonInfo["keepORF"],exonInfo["position"]))
		
	def loopChange(self):
		stats = {}

		stats["loopsChange"] = {}
		stats["loopsChange"].setdefault("different", 0)
		stats["loopsChange"].setdefault("same", 0)
		stats["loopsChange"].setdefault("onlyN", 0)
		stats["loopsChange"].setdefault("onlyT", 0)
		stats["loopsChange"].setdefault("noLoops", 0)

		total = 0

		for gene,info,switchDict,switch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network):
			nIso = switch.nTranscript
			tIso = switch.tTranscript

			nLoops = self._transcript_network._net.node[nIso.name]["iLoopsFamily"]
			tLoops = self._transcript_network._net.node[tIso.name]["iLoopsFamily"]

			if nLoops or tLoops:
				if nLoops == tLoops: 
					stats["loopsChange"]["same"] += 1
				elif nLoops:
					if tLoops:
						stats["loopsChange"]["different"] += 1
					else: 
						stats["loopsChange"]["onlyN"]	 += 1
				elif tLoops:
					stats["loopsChange"]["onlyT"]+= 1
			else:
				stats["loopsChange"]["noLoops"]	+=1

			total += 1

		with open("{0}result_summary/structural_loops.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tMapped_loops\tDifferent\t")
			F.write("Same\tOnly_nIso\tOnly_tIso\tNone\tTotal\n")
			F.write("{0}\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(stats["loopsChange"]["different"],stats["loopsChange"]["same"]) )
			F.write("{0}\t{1}\t".format(stats["loopsChange"]["onlyN"],stats["loopsChange"]["onlyT"]) )
			F.write("{0}\t{1}\n".format(stats["loopsChange"]["noLoops"],total ))

	def functionalChange(self):

		total = 0
		withFunctionalChanges = 0
		changeList = []
		
		PRINTS = 0
		ProSiteProfiles = 0
		ProSitePatterns = 0
		Pfam = 0
		totalType = 0

		mostCommonChanges = {}

		changes = [ z.functionalChange for w,x,y,z in self._gene_network.iterate_switches_ScoreWise(self._transcript_network) ]
		
		for x in changes:
			if x:
				withFunctionalChanges += 1
				for y in x:
					changeList.extend([y])
			total += 1

		for x in changeList:
			if x[0] == "PRINTS":
				PRINTS += x[3]
			elif x[0] == "ProSiteProfiles":
				ProSiteProfiles += x[3]
			elif x[0] == "Pfam":
				Pfam += x[3]
			elif x[0] == "ProSitePatterns":
				Pfam += x[3]
			totalType += x[3]

			mostCommonChanges.setdefault(x[2],0)
			mostCommonChanges[x[2]] += x[3]

		with open("{0}result_summary/structural_function.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tAnalysis\tYes\tNo\tTotal\n")
			F.write("{0}\tfunctional_change\t".format(options.Options().tag))
			F.write("{0}\t{1}\t".format(withFunctionalChanges,total-withFunctionalChanges))
			F.write("{0}\n".format(total))

			F.write("Cancer\tAnalysis\tPRINTS\tProSiteProfiles\tPfam\tProSitePatterns\tTotal\n")
			F.write("{0}\tchange_type\t".format(options.Options().tag))
			F.write("{0}\t{1}\t".format(PRINTS,ProSiteProfiles))
			F.write("{0}\t{1}\t".format(Pfam,ProSitePatterns))
			F.write("{0}\n".format(totalType))

			F.write("Cancer\tAnalysis\tFeature\tCount\n")
			for x in mostCommonChanges:
				F.write("{0}\tmotifs\t".format(options.Options().tag))
				F.write("{0}\t{1}\n".format(x,mostCommonChanges[x]))


	def disorderChange(self):
		allMotifs = []

		withDisorder = 0
		total = 0

		motifs = [ z.disorderChange for w,x,y,z in self._gene_network.iterate_switches_ScoreWise(self._transcript_network) ]

		for x in motifs:
			if x:
				withDisorder += 1
				for y in x:
					allMotifs.extend([y])
			total += 1

		meanLength = float(sum([ len(x) for x in allMotifs ])/len(allMotifs)) if len(allMotifs) > 0 else float('nan')

		with open("{0}result_summary/structural_disorder.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tAnalysis\tYes\tNo\tTotal\n")
			F.write("{0}\tdisordered_change\t".format(options.Options().tag))
			F.write("{0}\t{1}\t".format(withDisorder,total-withDisorder))
			F.write("{0}\n".format(total))

			F.write("Cancer\tAnalysis\tMean_length\n")
			F.write("{0}\tmean_length\t".format(options.Options().tag))
			F.write("{0}\n".format(meanLength))

	def printRelevantGene(self):

		affectedGenes = {}

		patientsNGenes = [ [x["symbol"],z.patients] for w,x,y,z in self._gene_network.iterate_switches_ScoreWise(self._transcript_network) ]

		patients = []
		[ patients.extend(z) for x,z in patientsNGenes ]
		patients = set(patients)
		
		genes = set([x for x,y in patientsNGenes])

		for thatGene in genes:
			affectedGenes.setdefault(thatGene,{})
			for thisGene,pats in patientsNGenes:
				if thisGene != thatGene:
					continue

				for pat in pats:
					affectedGenes[thatGene][pat] = 1

		with open("{0}result_summary/relevantSwitches.tsv".format(options.Options().qout), "w" ) as F:
			colnames = "Genes"
			for patient in patients:
				colnames += "\t" + options.Options().tag
			
			F.write("{0}\n".format(colnames))
			
			for gene in affectedGenes:
				F.write(gene)
				for patient in patients:
					F.write("\t"+str(affectedGenes[gene].get(patient,0)))
				F.write("\n")

	def topNeighborhoods(self):
		import pdb
		pdb.set_trace()
		neighborhoodAnalysis = []
		[ neighborhoodAnalysis.extend(z.neighborhoodChange) for w,x,y,z in self._gene_network.iterate_switches_ScoreWise(self._transcript_network) ]
		neighborhoodAnalysis = set(neighborhoodAnalysis)
		
		with open("{0}result_summary/neighborhoods_genesets.tsv".format(options.Options().qout,analysis), "w" ) as F:
			F.write("Cancer\tAnalysis\tGeneset\n")
			[ F.write("{0}\tenriched_genesets\t{1}\n".format(options.Options().tag,x)) for x in neighborhoodAnalysis ]

		return
		expressedGenes = [ (x,y) for x,y in self._gene_network.nodes(data=True) if y["ExpressedTranscripts"]]

		patients = options.Options().replicates
		if options.Options().unpairedReplicates:
			patients = options.Options().unpairedReplicates

		for analysis in neighborhoodAnalysis:

			geneSets = []
			[ geneSets.extend(x["neighborhoods"][analysis]) for w,x,y,z in self._gene_network.iterate_switches_ScoreWise(self._transcript_network) ]
			geneSets = set(geneSets)

			with open("{0}result_summary/neighborhoods_genes_{1}.tsv".format(options.Options().qout,analysis), "w" ) as F:
				F.write("Cancer\tAnalysis")
				[ F.write("\t{0}".format(gene)) for gene in expressedGenes ]
				F.write("\n")

				for geneSet in geneSets:
					F.write("{0}\tgene_set_enrichmen\t{2}".format(options.Options().tag,geneSet))
					for x,y in expressedGenes:
						if geneSet in geneSets["neighborhoods"][analysis]:
							F.write("\t1")
						else:
							F.write("\t0")
					F.write("\n")

			with open("{0}result_summary/neighborhoods_patients_{1}.tsv".format(options.Options().qout,analysis), "w" ) as F:
				F.write("GeneSet")
				[ F.write("\t{0}".format(patient)) for patient in patients ]
				F.write("\n")

				for geneSet in geneSets:
					F.write("{0}".format(geneSet))
					for patient in patients:
						for gene,info in self._gene_network.nodes(data=True):
							if patient in [ x.patients for x in info["isoformSwitches"] ]:
								F.write("\t1")
							else:
								F.write("\t0")
						F.write("\n")
