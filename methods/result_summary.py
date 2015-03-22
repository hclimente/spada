#!/soft/devel/python-2.7/bin/python

from libs import options
from libs import utils
from methods import method

import cPickle
import fisher
from itertools import groupby
from operator import itemgetter

import pdb

class ResultSummary(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self._random_gene_network = cPickle.load(open("{0}randomGeneNetwork.pkl".format(options.Options().qout)))

		self.proteinStats = { "Random": { "nIsoLength": [], "tIsoLength": [] },
							  "NonRandom": { "nIsoLength": [], "tIsoLength": [] }}

		self.switchStats = {}
		
		# driver/driver+relevance tests
		self.switchStats["driverTest"] = { "Driver": {"NoSwitch": 0, "Switch": 0}, 
										   "NonDriver": {"NoSwitch": 0, "Switch": 0} }
		self.switchStats["driverRTest"] = { "Driver": {"NonRelevant": 0, "Relevant": 0}, 
										   "NonDriver": {"NonRelevant": 0, "Relevant": 0} }

		# relevance test
		self.switchStats["relevanceTest"] = { "Random": {"NonRelevant": 0, "Relevant": 0}, 
											  "NonRandom": {"NonRelevant": 0, "Relevant": 0} }

		# cds and utr tests
		self.switchStats["cdsTest"] = { "Random": {"Change": 0, "NoChange": 0}, 
										"NonRandom": {"Change": 0, "NoChange": 0} }
		self.switchStats["utrTest"] = { "Random": {"Change": 0, "NoChange": 0}, 
										"NonRandom": {"Change": 0, "NoChange": 0} }
		self.switchStats["cdsPresence"] = { "Random": 
											{ "Both": 0, "OnlyN": 0, "OnlyT": 0, "None": 0 }, 
											"NonRandom": 
											{ "Both": 0, "OnlyN": 0, "OnlyT": 0, "None": 0 } 
										  }

		self.exonStats  = { "Random": [], "NonRandom": [] }

		self.loops = {}
		self.loops["loopsChange"] = {}
		self.loops["loopsChange"].setdefault("different", 0)
		self.loops["loopsChange"].setdefault("same", 0)
		self.loops["loopsChange"].setdefault("onlyN", 0)
		self.loops["loopsChange"].setdefault("onlyT", 0)
		self.loops["loopsChange"].setdefault("noLoops", 0)

		self.featuresTable = []

	def run(self):
		self.logger.info("Summarizing results.")

		txDict = self._transcript_network.nodes(data=True)

		for gene,info,switchDict,thisSwitch in self._random_gene_network.sampleSwitches(self._transcript_network,partialCreation=True):
			self.switchAndExonOverview(True,gene,info,switchDict,thisSwitch)
			self.changedFeatures(True,gene,info,switchDict,thisSwitch)

		# tests at switch level
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,options.Options().onlyModels,partialCreation=True,removeNoise=True):
			self.logger.debug("Getting statistics for switch {0}_{1}_{2}.".format(gene,thisSwitch.nTx,thisSwitch.tTx))

			self.switchAndExonOverview(False,gene,info,switchDict,thisSwitch)
			self.changedFeatures(False,gene,info,switchDict,thisSwitch)

			# general protein, switch and gene info
			# self.proteinOverview(switchDict,txDict)
			
			# structural info
			
			# self.loopChange(switch)
			#self.disorderChange(switch)
			#self.functionalChange(switch)

		# tests at gene level
		for gene,info in self._gene_network.iterate_genes_ScoreWise():
			dicts = [ x for x in info["isoformSwitches"] if x["model"] ]

			switchDict = None
			thisSwitch = None
			if dicts:
				switchDict = dicts[0]
				thisSwitch = self._gene_network.createSwitch(switchDict,self._transcript_network,True)

			self.driverTests(gene,info,switchDict,thisSwitch)

		self.printSwitchInfo()
		#self.printStructutalInfo()
		
		#self.functionalChange()
		#self.printRelevantGene()	

	def neighborhoodInfo(self):

		self.logger.info("Neighborhood information.")

		self.topNeighborhoods()

	def proteinOverview(self,random,switchDict,txDict):

		if random: randomTag = "Random"
		else: randomTag = "NonRandom"
		
		if [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["nIso"] ]!=[None]:
			seq = [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["nIso"] ][0]
			self.proteinStats[randomTag]["nIsoLength"].append(len(seq))
		if [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["tIso"] ]!=[None]:
			seq = [ x[1]["proteinSequence"] for x in txDict if x[0]==switchDict["tIso"] ][0]
			self.proteinStats[randomTag]["tIsoLength"].append(len(seq))

	def driverTests(self,gene,info,switchDict,thisSwitch):
		if info["Driver"]: 	driverTag = "Driver"
		else:  				driverTag = "NonDriver"

		if switchDict is None: 	
			switchTag = "NoSwitch"
		else:					
			switchTag = "Switch"
			if thisSwitch.is_relevant:
				relevanceTag = "Relevant"
			else:
				relevanceTag = "NonRelevant"
			
		self.switchStats["driverTest"][driverTag][switchTag] += 1
		if switchDict is not None:
			self.switchStats["driverRTest"][driverTag][relevanceTag] += 1

	def switchAndExonOverview(self,random,gene,info,switchDict,thisSwitch):
		nIso = thisSwitch.nTranscript
		tIso = thisSwitch.tTranscript

		if random: 	randomTag = "Random"
		else:		randomTag = "NonRandom"

		if thisSwitch.is_relevant:	relevanceTag = "Relevant"
		else:						relevanceTag = "NonRelevant"

		if thisSwitch.cds_diff: cdsTag = "Change"
		else:					cdsTag = "NoChange"

		if thisSwitch.utr_diff: utrTag = "Change"
		else: 					utrTag = "NoChange"

		if nIso.cds and tIso.cds: presenceTag = "Both"
		elif nIso.cds: 			  presenceTag = "OnlyN"
		elif tIso.cds: 			  presenceTag = "OnlyT"
		else: 					  presenceTag = "None"

		self.switchStats["relevanceTest"][randomTag][relevanceTag] += 1
		self.switchStats["cdsTest"][randomTag][cdsTag] += 1
		self.switchStats["utrTest"][randomTag][utrTag] += 1
		self.switchStats["cdsPresence"][randomTag][presenceTag] += 1

		'''
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
				exonInfo["random"] = randomTag
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

				self.exonStats[randomTag].append(exonInfo)
		'''
		
	def printSwitchInfo(self):

		with open("{0}result_summary/nIso_length_nonRandom{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			for length in self.proteinStats["NonRandom"]["nIsoLength"]:
				F.write("{0}\n".format(length))
		
		with open("{0}result_summary/tIso_length_nonRandom{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			for length in self.proteinStats["NonRandom"]["tIsoLength"]:
				F.write("{0}\n".format(length))

		with open("{0}result_summary/nIso_length_random{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			for length in self.proteinStats["Random"]["nIsoLength"]:
				F.write("{0}\n".format(length))
		
		with open("{0}result_summary/tIso_length_random{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			for length in self.proteinStats["Random"]["tIsoLength"]:
				F.write("{0}\n".format(length))

		with open("{0}result_summary/switches{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tAnalysis\tBoth\tOnly_nIso\tOnly_tIso\tNone\t")
			F.write("Random_Both\tRandom_Only_nIso\tRandom_Only_tIso\tRandom_None\n")
			F.write("{0}\tCDS_study\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["Both"]))
			F.write("{0}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["OnlyN"]))
			F.write("{0}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["OnlyT"]))
			F.write("{0}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["None"]))
			F.write("{0}\t".format(self.switchStats["cdsPresence"]["Random"]["Both"]))
			F.write("{0}\t".format(self.switchStats["cdsPresence"]["Random"]["OnlyN"]))
			F.write("{0}\t".format(self.switchStats["cdsPresence"]["Random"]["OnlyT"]))
			F.write("{0}\n".format(self.switchStats["cdsPresence"]["Random"]["None"]))

			F.write("Cancer\tAnalysis\tYes\tNo\tRandom_Yes\tRandom_No\tp\n")

			F.write("{0}\tCDS_change\t".format(options.Options().tag))
			F.write("{0}\t".format(self.switchStats["cdsTest"]["NonRandom"]["Change"]))
			F.write("{0}\t".format(self.switchStats["cdsTest"]["NonRandom"]["NoChange"]))
			F.write("{0}\t".format(self.switchStats["cdsTest"]["Random"]["Change"]))
			F.write("{0}\n".format(self.switchStats["cdsTest"]["Random"]["NoChange"]))

			p = fisher.pvalue(self.switchStats["cdsTest"]["NonRandom"]["Change"],
							self.switchStats["cdsTest"]["NonRandom"]["NoChange"],
							self.switchStats["cdsTest"]["Random"]["Change"],
							self.switchStats["cdsTest"]["Random"]["NoChange"])

			F.write("{0}\n".format(p.two_tail) )

			F.write("{0}\tUTR_change\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["utrTest"]["NonRandom"]["Change"]))
			F.write("{0}\t".format(self.switchStats["utrTest"]["NonRandom"]["NoChange"]))
			F.write("{0}\t".format(self.switchStats["utrTest"]["Random"]["Change"]))
			F.write("{0}\n".format(self.switchStats["utrTest"]["Random"]["NoChange"]))

			p = fisher.pvalue(self.switchStats["utrTest"]["NonRandom"]["Change"],
							self.switchStats["utrTest"]["NonRandom"]["NoChange"],
							self.switchStats["utrTest"]["Random"]["Change"],
							self.switchStats["utrTest"]["Random"]["NoChange"])

			F.write("{0}\n".format(p.two_tail) )

			##### GENE LEVEL #####
			F.write("Cancer\tAnalysis\tFeat-Switch\tFeat-NoSwitch\tNoFeat-Switch\tNoFeat-NoSwitch\tp\n")
			# driver enrichment in switches
			F.write("{0}\tDriver_enrichment\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["driverTest"]["Driver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["driverTest"]["Driver"]["NoSwitch"]))
			F.write("{0}\t".format(self.switchStats["driverTest"]["NonDriver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["driverTest"]["NonDriver"]["NoSwitch"]))
			
			p = fisher.pvalue(self.switchStats["driverTest"]["Driver"]["Switch"],
							self.switchStats["driverTest"]["Driver"]["NoSwitch"],
							self.switchStats["driverTest"]["NonDriver"]["Switch"],
							self.switchStats["driverTest"]["NonDriver"]["NoSwitch"])

			F.write("{0}\n".format(p.two_tail) )

			# driver enrichment in relevant
			F.write("{0}\tDriver_relevance_enrichment\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["driverRTest"]["Driver"]["Relevant"]) )
			F.write("{0}\t".format(self.switchStats["driverRTest"]["Driver"]["NonRelevant"]) )
			F.write("{0}\t".format(self.switchStats["driverRTest"]["NonDriver"]["Relevant"]) )
			F.write("{0}\t".format(self.switchStats["driverRTest"]["NonDriver"]["NonRelevant"]) )
			
			p = fisher.pvalue(self.switchStats["driverRTest"]["Driver"]["Relevant"],
							self.switchStats["driverRTest"]["Driver"]["NonRelevant"],
							self.switchStats["driverRTest"]["NonDriver"]["Relevant"],
							self.switchStats["driverRTest"]["NonDriver"]["NonRelevant"])

			F.write("{0}\n".format(p.two_tail) )
		'''
		with open("{0}result_summary/exons{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tSwitch\tLength\tType\tKeepOrf\tPosition\n")
			for exon in self.exonStats:
				F.write("{0}\t{1}\t".format(options.Options().tag,exon["switch"]))
				F.write("{0}\t{1}\t".format(exon["length"],exon["role"]))
				F.write("{0}\t{1}\n".format(exon["keepORF"],exon["position"]))
		'''

	def printStructutalInfo(self):

		with open("{0}result_summary/structural_summary{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tGene\tSymbol\tnTx\ttTx\tPfam\t")
			F.write("IUPREDLong\tIUPREDShort\tRelevant\tModel\t")
			F.write("Noise\tDriver\tASDriver\tDriverType\n")
			for tag,random,featureDict in self.featuresTable:
				switchElements = tag.split("_")
				gene = switchElements[0]
				symbol = switchElements[1]
				nTx = switchElements[2]
				tTx = switchElements[3]

				F.write("{0}\t{1}\t".format(options.Options().tag,gene))
				F.write("{0}\t{1}\t{2}\t".format(symbol,nTx,tTx))
				F.write("{0}\t{1}\t".format(len(featureDict["Pfam"]),len(featureDict["IUPREDLong"])))
				F.write("{0}\t".format(len(featureDict["IUPREDShort"])))
				F.write("{0}\t{1}\t".format(featureDict["Relevant"],featureDict["Model"]))
				F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["Driver"]))
				F.write("{0}\t{1}\n".format(featureDict["ASDriver"],featureDict["DriverType"]))

		with open("{0}result_summary/structural_features{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tGene\tSymbol\tnTx\ttTx\tAnalysis\tWhatsHappenning\t")
			F.write("Feature\tRelevant\tModel\tNoise\tDriver\tASDriver\tDriverType\n")
			
			for s in self.featuresTable:
				switchElements = s.split("_")
				gene = switchElements[0]
				symbol = switchElements[1]
				nTx = switchElements[2]
				tTx = switchElements[3]

				for pfamDom in featureDict["Pfam"]:
					F.write("{0}\t{1}\tPfam\t".format(options.Options().tag,gene))
					F.write("{0}\t{1}\t{2}\t".format(symbol,nTx,tTx))
					F.write("{0}\t{1}\t".format(pfamDom[1],pfamDom[0]))
					F.write("{0}\t{1}\t".format(featureDict["Relevant"],featureDict["Model"]))
					F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["Driver"]))
					F.write("{0}\t{1}\n".format(featureDict["ASDriver"],featureDict["DriverType"]))

				if featureDict["IUPREDLong"]:
					for iulong in featureDict["IUPREDLong"][0]:
						F.write("{0}\t{1}\tIUPREDLong\t".format(options.Options().tag,gene))
						F.write("{0}\t{1}\t{2}\t".format(symbol,nTx,tTx))

						F.write("{0}\t{1}\t".format(featureDict["IUPREDLong"][1],iulong))
						F.write("{0}\t{1}\t".format(featureDict["Relevant"],featureDict["Model"]))
						F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["Driver"]))
						F.write("{0}\t{1}\n".format(featureDict["ASDriver"],featureDict["DriverType"]))

				if featureDict["IUPREDShort"]:
					for iushort in featureDict["IUPREDShort"][0]:
						F.write("{0}\t{1}\tIUPREDShort\t".format(options.Options().tag,gene))
						F.write("{0}\t{1}\t{2}\t".format(symbol,nTx,tTx))
						F.write("{0}\t{1}\t".format(featureDict["IUPREDShort"][1],iushort))
						F.write("{0}\t{1}\t".format(featureDict["Relevant"],featureDict["Model"]))
						F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["Driver"]))
						F.write("{0}\t{1}\n".format(featureDict["ASDriver"],featureDict["DriverType"]))

		with open("{0}result_summary/structural_loops{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tDifferent\t")
			F.write("Same\tOnly_nIso\tOnly_tIso\tNone\n")
			F.write("{0}\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["different"],self.loops["loopsChange"]["same"]) )
			F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["onlyN"],self.loops["loopsChange"]["onlyT"]) )
			F.write("{0}\t{1}\n".format(self.loops["loopsChange"]["noLoops"]))

	def changedFeatures(self,random,gene,info,switchDict,switch):	

		tag = "{0}_{1}_{2}_{3}".format(gene,info["symbol"],switch.nTx,switch.tTx)

		if random:	randomTag = "Random"
		else:	randomTag = "NonRandom"

		switchFeatures = {}

		if switch.functionalChange:
			Pfam = []

			functionalChanges = ""
			for element in utils.readTable("{0}structural_analysis/interpro_analysis.tsv".format(options.Options().qout)):
				if element[0]==gene and element[2] in [switch.nTx,switch.tTx]:
					if element[4]=="Pfam":
						Pfam.append((element[6],element[3]))

			switchFeatures["Pfam"] = Pfam
			
		else:
			switchFeatures["Pfam"] = []
		
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
			
			switchFeatures["IUPREDShort"] = shortDisorder
			switchFeatures["IUPREDLong"] = longDisorder
		else:
			switchFeatures["IUPREDShort"] = []
			switchFeatures["IUPREDLong"] = []
		
		switchFeatures["Driver"] = int(info["Driver"])
		switchFeatures["ASDriver"] = int(info["ASDriver"])
		switchFeatures["DriverType"] = info["DriverType"]
		switchFeatures["Relevant"] = int(switch.is_relevant)
		switchFeatures["Model"] = int(switchDict["model"])
		switchFeatures["Noise"] = int(switchDict["noise"])

		self.featuresTable.append((tag,randomTag,switchFeatures))

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
