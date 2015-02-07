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

		self.sizes = {}

		self.proteinStats = {}
		self.proteinStats["centrality"] = []
		self.proteinStats["nIsoLength"] = []
		self.proteinStats["tIsoLength"] = []

		self.switchStats = {}
		self.switchStats["hasCds"] = { "both": 0, "onlyN": 0,"onlyT": 0,"none": 0 }
		self.switchStats["hasCdsChange"] = { "y": 0, "n": 0 }
		self.switchStats["hasUtrChange"] = { "y": 0, "n": 0 }
		self.switchStats["isDriver"] = { "y": 0, "n": 0 }
		self.switchStats["isRelevant"] = { "y": 0, "n": 0 }

		self.exonStats = []

		self.loops = {}
		self.loops["loopsChange"] = {}
		self.loops["loopsChange"].setdefault("different", 0)
		self.loops["loopsChange"].setdefault("same", 0)
		self.loops["loopsChange"].setdefault("onlyN", 0)
		self.loops["loopsChange"].setdefault("onlyT", 0)
		self.loops["loopsChange"].setdefault("noLoops", 0)

		self.featuresTable = {}

	def run(self):
		self.logger.info("Summarizing results.")

		txDict = self._transcript_network.nodes(data=True)
		total = 0

		for gene,info,switchDict,switch in self._gene_network.iterate_nonRelevantSwitches_ScoreWise(self._transcript_network,partialCreation=True):
			self.logger.debug("Getting statistics for switch {0}_{1}_{2}.".format(gene,switch.nTx,switch.tTx))

			# general protein, switch and gene info
			self.proteinOverview(switchDict,txDict)
			self.switchAndExonOverview(gene,info,switchDict,switch)

			# structural info
			self.changedFeatures(gene,info,switch)
			self.loopChange(switch)
			#self.disorderChange(switch)
			#self.functionalChange(switch)
			
			total += 1

		sortedNodes = sorted(self._gene_network.nodes(data=True), 
							 key=lambda (a, dct): dct['score'], reverse=True)
		for gene,info in sortedNodes:
			self.proteinStats["centrality"].append(self._gene_network._net.degree(gene))

		self.printSwitchInfo(total)
		self.printStructutalInfo(total)
		
		#self.functionalChange()
		#self.printRelevantGene()	

	def neighborhoodInfo(self):

		self.logger.info("Neighborhood information.")

		self.topNeighborhoods()

	def proteinOverview(self,switchDict,txDict):
		
		if [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["nIso"] ]!=[None]:
			seq = [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["nIso"] ][0]
			self.proteinStats["nIsoLength"].append(len(seq))
		if [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["tIso"] ]!=[None]:
			seq = [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["tIso"] ][0]
			self.proteinStats["tIsoLength"].append(len(seq))

	def switchAndExonOverview(self,gene,info,switchDict,thisSwitch):

		nIso = thisSwitch.nTranscript
		tIso = thisSwitch.tTranscript

		if nIso.cds and tIso.cds: self.switchStats["hasCds"]["both"]  += 1
		elif nIso.cds: 			  self.switchStats["hasCds"]["onlyN"] += 1
		elif tIso.cds: 			  self.switchStats["hasCds"]["onlyT"] += 1
		else: 					  self.switchStats["hasCds"]["none"]  += 1

		if thisSwitch.cds_diff: self.switchStats["hasCdsChange"]["y"] += 1
		else:					self.switchStats["hasCdsChange"]["n"] += 1

		if thisSwitch.utr_diff: self.switchStats["hasUtrChange"]["y"] += 1
		else: 				self.switchStats["hasUtrChange"]["n"] += 1

		if info["Driver"]:  self.switchStats["isDriver"]["y"] +=1
		else: 				self.switchStats["isDriver"]["n"] +=1


		if thisSwitch.is_relevant:  self.switchStats["isRelevant"]["y"] +=1
		else: 						self.switchStats["isRelevant"]["n"] +=1

		nSpecificCds = [ x for x in thisSwitch._normal_transcript.cds if x not in thisSwitch._tumor_transcript.cds ]
		nSpecificUtr = [ x for x in thisSwitch._normal_transcript.utr if x not in thisSwitch._tumor_transcript.utr ]
		tSpecificCds = [ x for x in thisSwitch._tumor_transcript.cds if x not in thisSwitch._normal_transcript.cds ]
		tSpecificUtr = [ x for x in thisSwitch._tumor_transcript.utr if x not in thisSwitch._normal_transcript.utr ]
		
		for specificCds,specificUtr,cds in zip([nSpecificCds,tSpecificCds],[nSpecificUtr,tSpecificUtr],[thisSwitch._normal_transcript.cds,thisSwitch._tumor_transcript.cds]):
		
			specific = specificUtr
			specific.extend(specificCds)
			specific = sorted(specific)
			
			setSpecCds = set(specificCds)
			setSpecUtr = set(specificUtr)
			setSpec = set(specific)

			exons = [ map(itemgetter(1),g) for k,g in groupby(enumerate(specific), lambda (i,x): i-x) ]
			
			for exon in exons:
				setExon = set(exon)
				exonInfo = {}
				exonInfo["switch"] = "{0}_{1}_{2}".format(gene,switchDict["nIso"],switchDict["tIso"])
				exonInfo["length"] = len(exon)

				if setExon & setSpecCds and setExon & setSpecUtr:
					exonicCds = sorted(setExon & setSpecCds)
					exonInfo["role"] = "CDS-UTR"
					exonInfo["keepORF"] = len(exonicCds)%3
					medianPos = exonicCds[len(exonicCds)/2]
					pos = [ i for i,x in enumerate(cds) if x==medianPos ][0]
					exonInfo["position"] = float(pos)/len(cds)
				elif setExon & setSpecCds:
					exonInfo["role"] = "CDS"
					exonInfo["keepORF"] = len(switchAndExonOverview)%3
					medianPos = exon[len(exon)/2]
					pos = [ i for i,x in enumerate(cds) if x==medianPos ][0]
					exonInfo["position"] = float(pos)/len(cds)
				elif setExon & setSpecUtr:
					exonInfo["role"] = "UTR"
					exonInfo["keepORF"] = None
					exonInfo["position"] = None
				else:
					exonInfo["role"] = "wut"
					exonInfo["keepORF"] = "wut"
					exonInfo["position"] = "wut"

				self.exonStats.append(exonInfo)
		
	def printSwitchInfo(self,total):
		with open("{0}result_summary/protein_centrality.tsv".format(options.Options().qout), "w" ) as F:
			for degree in self.proteinStats["centrality"]:
				F.write("{0}\t{1}\n".format(options.Options().tag,degree))

		with open("{0}result_summary/nIso_length.tsv".format(options.Options().qout), "w" ) as F:
			for length in self.proteinStats["nIsoLength"]:
				F.write("{0}\n".format(length))
		
		with open("{0}result_summary/tIso_length.tsv".format(options.Options().qout), "w" ) as F:
			for length in self.proteinStats["tIsoLength"]:
				F.write("{0}\n".format(length))

		with open("{0}result_summary/switches.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tAnalysis\tBoth\tOnly_nIso\tOnly_tIso\tNone\tTotal\n")
			F.write("{0}\tCDS_study\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(self.switchStats["hasCds"]["both"],self.switchStats["hasCds"]["onlyN"]) )
			F.write("{0}\t{1}\t".format(self.switchStats["hasCds"]["onlyT"],self.switchStats["hasCds"]["none"]) )
			F.write("{0}\n".format(total))

			F.write("Cancer\tAnalysis\tYes\tNo\tTotal\n")
			F.write("{0}\tCDS_change\t".format(options.Options().tag))
			F.write("{0}\t{1}\t".format(self.switchStats["hasCdsChange"]["y"],self.switchStats["hasCdsChange"]["n"]) )
			F.write("{0}\n".format(total))

			F.write("{0}\tUTR_change\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(self.switchStats["hasUtrChange"]["y"],self.switchStats["hasUtrChange"]["n"]) )
			F.write("{0}\n".format(total))

			F.write("{0}\tDriver_affection\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(self.switchStats["isDriver"]["y"],self.switchStats["isDriver"]["n"]) )
			F.write("{0}\n".format(total))

			F.write("{0}\tRelevant\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(self.switchStats["isRelevant"]["y"],self.switchStats["isRelevant"]["n"]) )
			F.write("{0}\n".format(total))

		with open("{0}result_summary/exons.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tSwitch\tLength\tType\tKeepOrf\tPosition\n")
			for exon in self.exonStats:
				F.write("{0}\t{1}\t".format(options.Options().tag,exon["switch"]))
				F.write("{0}\t{1}\t".format(exon["length"],exon["role"]))
				F.write("{0}\t{1}\n".format(exon["keepORF"],exon["position"]))

	def printStructutalInfo(self,total):
		with open("{0}result_summary/structural_summary.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tSwitch\tPfam\tPRINTS\tProSitePatterns\t")
			F.write("IUPREDLong\tIUPREDShort\tI3D\tDriver\n")
			for s in self.featuresTable:
				F.write("{0}\t{1}\t".format(options.Options().tag,s))
				F.write("{0}\t{1}\t".format(len(self.featuresTable[s][0]),len(self.featuresTable[s][1])))
				F.write("{0}\t{1}\t".format(len(self.featuresTable[s][2]),len(self.featuresTable[s][3])))
				F.write("{0}\t{1}\t".format(len(self.featuresTable[s][4]),len(self.featuresTable[s][5])))
				F.write("{0}\n".format(self.featuresTable[s][6]))

		with open("{0}result_summary/structural_features.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tSwitch\tAnalysis\tAction\tFeature\tDriver\n")
			
			for s in self.featuresTable:
				Pfam = self.featuresTable[s][0]
				PRINTS = self.featuresTable[s][1]
				ProSitePatterns = self.featuresTable[s][2]
				IUPREDLong = self.featuresTable[s][3]
				IUPREDShort = self.featuresTable[s][4]
				I3D = self.featuresTable[s][5]
				Driver = self.featuresTable[s][6]

				for pfamDom in Pfam:
					F.write("{0}\t{1}\tPfam\t".format(options.Options().tag,s))
					F.write("{0}\t{1}\t{2}\n".format(pfamDom[1],pfamDom[0],Driver))

				for prints in PRINTS:
					F.write("{0}\t{1}\tPRINTS\t".format(options.Options().tag,s))
					F.write("{0}\t{1}\t{2}\n".format(prints[1],prints[0],Driver))

				for patt in ProSitePatterns:
					F.write("{0}\t{1}\tProSitePatterns\t".format(options.Options().tag,s))
					F.write("{0}\t{1}\t{2}\n".format(patt[1],patt[0],Driver))

				if IUPREDLong:
					for iulong in IUPREDLong[0]:
						F.write("{0}\t{1}\tIUPREDLong\t".format(options.Options().tag,s))
						F.write("{0}\t{1}\t{2}\n".format(IUPREDLong[1],iulong,Driver))

				if IUPREDShort:
					for iushort in IUPREDShort[0]:
						F.write("{0}\t{1}\tIUPREDShort\t".format(options.Options().tag,s))
						F.write("{0}\t{1}\t{2}\n".format(IUPREDShort[1],iushort,Driver))

				for i3d in I3D:
					F.write("{0}\t{1}\tI3D\t".format(options.Options().tag,s))
					F.write("{0}\t{1}\t{2}\n".format(i3d[1],i3d[0],Driver))

		with open("{0}result_summary/structural_loops.tsv".format(options.Options().qout), "w" ) as F:
			F.write("Cancer\tDifferent\t")
			F.write("Same\tOnly_nIso\tOnly_tIso\tNone\tTotal\n")
			F.write("{0}\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["different"],self.loops["loopsChange"]["same"]) )
			F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["onlyN"],self.loops["loopsChange"]["onlyT"]) )
			F.write("{0}\t{1}\n".format(self.loops["loopsChange"]["noLoops"],total ))

	def changedFeatures(self,gene,info,switch):

		tag = "{0}_{1}_{2}".format(gene,switch.nTx,switch.tTx)

		self.featuresTable.setdefault(tag,[])

		if switch.functionalChange:
			Pfam = []
			PRINTS = []
			ProSitePatterns = []

			functionalChanges = ""
			for element in utils.readTable("{0}structural_analysis/InterPro_report.tsv".format(options.Options().qout)):
				if element[0]==gene and element[2] in [switch.nTx,switch.tTx]:
					if element[4]=="Pfam":
						Pfam.append((element[6],element[3]))
					elif element[4]=="PRINTS":
						PRINTS.append((element[6],element[3]))
					elif element[4]=="ProSitePatterns":
						ProSitePatterns.append((element[6],element[3]))

			self.featuresTable[tag].extend([Pfam,PRINTS,ProSitePatterns])
			
		else:
			self.featuresTable[tag].extend([[],[],[]])
		
		if switch.disorderChange:
			longDisorder = []
			shortDisorder = []

			for element in utils.readTable("{0}structural_analysis/iupred_analysis.tsv".format(options.Options().qout)):
				if element[0]==gene and element[2] in [switch.nTx,switch.tTx] and element[4]:
					if element[3]=='short':
						if element[2] == switch.nTx:
							shortDisorder = (element[4].split(","),'Lost in tumor')
						else:
							shortDisorder = (element[4].split(","),'Gained in tumor')
					elif element[3]=='long':
						if element[2] == switch.nTx:
							longDisorder = (element[4].split(","),'Lost in tumor')
						else:
							longDisorder = (element[4].split(","),'Gained in tumor')
			
			self.featuresTable[tag].extend([shortDisorder,longDisorder])
		else:
			self.featuresTable[tag].extend([[],[]])
		
		if switch.brokenSurfaces:
			surfaceChanges = []

			for element in utils.readTable("{0}structural_analysis/I3D_analysis.tsv".format(options.Options().qout)):
				if element[0]==gene and element[2] in [switch.nTx,switch.tTx] and element[4]:
					if element[2] == switch.nTx:
						surfaceChanges.append((element[4],"Lost in tumor"))
					else:
						surfaceChanges.append((element[4],"Gained in tumor"))
			
			self.featuresTable[tag].append(surfaceChanges)
		else:
			self.featuresTable[tag].append([])

		self.featuresTable[tag].append(info["Driver"])

	def loopChange(self,switch):
		nIso = switch.nTranscript
		tIso = switch.tTranscript

		nLoops = self._transcript_network._net.node[nIso.name]["iLoopsFamily"]
		tLoops = self._transcript_network._net.node[tIso.name]["iLoopsFamily"]

		if nLoops or tLoops:
			if nLoops == tLoops: 
				self.loops["loopsChange"]["same"] += 1
			elif nLoops:
				if tLoops:
					self.loops["loopsChange"]["different"] += 1
				else: 
					self.loops["loopsChange"]["onlyN"] += 1
			elif tLoops:
				self.loops["loopsChange"]["onlyT"] += 1
		else:
			self.loops["loopsChange"]["noLoops"] +=1

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
