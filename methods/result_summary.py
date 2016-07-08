from interface import out_network
from libs import options
from libs import utils
from methods import method

import pickle
from scipy.stats import fisher_exact
from itertools import groupby
import numpy as np
from operator import itemgetter
from scipy import stats

class ResultSummary(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self._random_gene_network = pickle.load(open("{0}randomGeneNetwork_fixNormal.pkl".format(options.Options().qout),"rb"))
		self._random_gene_network.createLogger()

		self.proteinStats = { "Random": [], "NonRandom": []}

		self.switchStats = {}
		
		# driver/driver+relevance tests
		self.switchStats["d0Enrichment"] = { "driver": {"NoSwitch": 0.0, "Switch": 0.0}, 
										   	 "nonDriver": {"NoSwitch": 0.0, "Switch": 0.0} }
		self.switchStats["d0Patients"] = { "driver": [], "nonDriver": [] }
		self.switchStats["d0PatientsFunctional"] = { "driver": [], "nonDriver": [] }
		self.switchStats["d1Enrichment"] = { "D1Driver": {"NoSwitch": 0.0, "Switch": 0.0}, 
										   	 "NonD1Driver": {"NoSwitch": 0.0, "Switch": 0.0} }
		self.switchStats["d1Patients"] = { "D1Driver": [], "NonD1Driver": [] }
		self.switchStats["d1PatientsFunctional"] = { "D1Driver": [], "NonD1Driver": [] }
		self.switchStats["d0Functional"] = { "driver": {"NonFunctional": 0.0, "Functional": 0.0}, 
										   	"nonDriver": {"NonFunctional": 0.0, "Functional": 0.0} }
		self.switchStats["d1Functional"] = { "D1Driver": {"NonFunctional": 0.0, "Functional": 0.0}, 
										   	"NonD1Driver": {"NonFunctional": 0.0, "Functional": 0.0} }

		# relevance test
		self.switchStats["functionality"] = { "Random": {"NonFunctional": 0.0, "Functional": 0.0}, 
											  "NonRandom": {"NonFunctional": 0.0, "Functional": 0.0} }

		# cds and utr tests
		self.switchStats["cdsTest"] = { "Random": {"Change": 0.0, "NoChange": 0.0}, 
										"NonRandom": {"Change": 0.0, "NoChange": 0.0} }
		self.switchStats["utrTest"] = { "Random": {"Change": 0.0, "NoChange": 0.0}, 
										"NonRandom": {"Change": 0.0, "NoChange": 0.0} }
		self.switchStats["cdsPresence"] = { "Random": 
											{ "Both": 0.0, "OnlyN": 0.0, "OnlyT": 0.0, "None": 0.0 }, 
											"NonRandom": 
											{ "Both": 0.0, "OnlyN": 0.0, "OnlyT": 0.0, "None": 0.0 } 
										  }

		self.exonStats  = { "Random": [], "NonRandom": [] }
		self.alternativeSplicingStats  = { "Random": [], "NonRandom": [] }

		self.loops = {}
		self.loops["loopsChange"] = {}
		self.loops["loopsChange"].setdefault("different", 0.0)
		self.loops["loopsChange"].setdefault("same", 0.0)
		self.loops["loopsChange"].setdefault("onlyN", 0.0)
		self.loops["loopsChange"].setdefault("onlyT", 0.0)
		self.loops["loopsChange"].setdefault("noLoops", 0.0)

		self.featuresTable = []

	def clean(self):
		o = options.Options()

		#utils.cmd("rm","-r",".testOld2/{0}/result_summary".format(o.out))
		#utils.cmd("mv",".testOld/{0}/result_summary".format(o.out),".testOld2/{0}/result_summary".format(o.out))
		#utils.cmd("mv","{0}/result_summary".format(o.qout),".testOld/{0}/result_summary".format(o.out))
		utils.cmd("mkdir", "{}result_summary".format(o.qout))
	
	def run(self):
		self.logger.info("Summarizing results.")

		txDict = self._transcript_network.nodes(data=True)

		# print out switches and random switches
		out_network.outCandidateList(self._gene_network,self._transcript_network)
		out_network.outCandidateList(self._random_gene_network,self._transcript_network,
			filename="random.candidateList_info.tsv")

		testedSwitches = 0

		# tests at switch level
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network,only_models=True,partialCreation=True,removeNoise=True):
			self.logger.debug("Getting statistics for switch {0}_{1}_{2}.".format(gene,thisSwitch.nTx,thisSwitch.tTx))

			# general protein, switch and gene info
			self.switchAndExonOverview(False,gene,info,switchDict,thisSwitch)
			self.proteinOverview(False,switchDict,txDict,thisSwitch)

			# structural info
			self.changedStructuralFeatures(False,gene,info,switchDict,thisSwitch)
			# self.loopChange(switch)
			# self.disorderChange(switch)

			testedSwitches += 1
			
		for gene,info,switchDict,thisSwitch in self._random_gene_network.sampleSwitches(self._transcript_network,numIterations=testedSwitches):
			self.switchAndExonOverview(True,gene,info,switchDict,thisSwitch)
			self.proteinOverview(True,switchDict,txDict,thisSwitch)
			self.changedStructuralFeatures(True,gene,info,switchDict,thisSwitch)

		# tests at gene level
		for gene,info in self._gene_network.iterate_genes_byPatientNumber(alwaysSwitchedGenes=True):
			dicts = [ x for x in info["isoformSwitches"] if x["model"] and not x["noise"] ]

			modelDict = None
			if dicts:
				modelDict = dicts[0]

			self.driverTests(gene,info,modelDict)

		self.printSwitchInfo()
		self.printStructutalInfo()
		
		#self.printFunctionalGene()	

	def neighborhoodInfo(self):

		self.logger.info("Neighborhood information.")

		self.topNeighborhoods()

	def proteinOverview(self,random,switchDict,txDict,thisSwitch):

		if random: randomTag = "Random"
		else: randomTag = "NonRandom"
		
		nLength = 0
		tLength = 0
		nSpLength = 0
		tSpLength = 0

		if thisSwitch.nIsoform:
			nLength = len(thisSwitch.nIsoform.seq)
			nIso = thisSwitch.nIsoform.getSegments('isoform-specific')
			nSp = []
			[nSp.extend(x) for x in nIso]
			nSpLength = len(nSp)
		if thisSwitch.tIsoform:
			tLength = len(thisSwitch.tIsoform.seq)
			tIso = thisSwitch.tIsoform.getSegments('isoform-specific')
			tSp = []
			[tSp.extend(x) for x in tIso]
			tSpLength = len(tSp)

		self.proteinStats[randomTag].append((nLength,tLength,nSpLength,tSpLength))

	def driverTests(self,gene,info,switchDict):
		if info["driver"]: 	driverTag = "driver"
		else:  				driverTag = "nonDriver"

		if switchDict is None: 	
			switchTag = "NoSwitch"
			functionalTag = "NonFunctional"
		else:					
			switchTag = "Switch"
			if switchDict["functional"]:
				functionalTag = "Functional"
			else:
				functionalTag = "NonFunctional"
		
		# driver enrichment
		self.switchStats["d0Enrichment"][driverTag][switchTag] += 1
		self.switchStats["d0Functional"][driverTag][functionalTag] += 1

		if not info["driver"]:
			d1 = [ x for x in self._gene_network._net.neighbors(gene) if self._gene_network._net.node[x]["driver"] ]
			if d1: 	
				d1DriverTag="D1Driver"
			else: 	
				d1DriverTag="NonD1Driver"

			# d1 enrichment
			self.switchStats["d1Enrichment"][d1DriverTag][switchTag] += 1
			self.switchStats["d1Functional"][d1DriverTag][functionalTag] += 1

		# number of patients enrichment
		if switchDict:
			self.switchStats["d0Patients"][driverTag].append(len(switchDict["patients"]))
			if switchDict["functional"]:
				self.switchStats["d0PatientsFunctional"][driverTag].append(len(switchDict["patients"]))
			
			if not info["driver"]:
				self.switchStats["d1Patients"][d1DriverTag].append(len(switchDict["patients"]))
				if switchDict["functional"]:
					self.switchStats["d1PatientsFunctional"][d1DriverTag].append(len(switchDict["patients"]))

	def switchAndExonOverview(self,random,gene,info,switchDict,thisSwitch):
		nTx = thisSwitch.nTranscript
		tTx = thisSwitch.tTranscript

		if random: 	randomTag = "Random"
		else:		randomTag = "NonRandom"

		if thisSwitch.is_functional:	functionalTag = "Functional"
		else:						functionalTag = "NonFunctional"

		if thisSwitch.cds_diff: cdsTag = "Change"
		else:					cdsTag = "NoChange"

		if thisSwitch.utr_diff: utrTag = "Change"
		else: 					utrTag = "NoChange"

		if nTx.cds and tTx.cds: presenceTag = "Both"
		elif nTx.cds: 			  presenceTag = "OnlyN"
		elif tTx.cds: 			  presenceTag = "OnlyT"
		else: 					  presenceTag = "None"

		self.switchStats["functionality"][randomTag][functionalTag] += 1
		self.switchStats["cdsTest"][randomTag][cdsTag] += 1
		self.switchStats["utrTest"][randomTag][utrTag] += 1
		self.switchStats["cdsPresence"][randomTag][presenceTag] += 1

		nSpecificCds = set([ x for x in nTx.cds if x not in tTx.cds ])
		nSpecificUtr = set([ x for x in nTx.utr if x not in tTx.utr ])
		tSpecificCds = set([ x for x in tTx.cds if x not in nTx.cds ])
		tSpecificUtr = set([ x for x in tTx.utr if x not in nTx.utr ])
		
		for specificCds,specificUtr,cds,origin in zip([nSpecificCds,tSpecificCds],[nSpecificUtr,tSpecificUtr],[nTx.cds,tTx.cds],["nIso","tIso"]):
		
			specificRegions = sorted(list(specificUtr | specificCds))
			exons = [ set(map(itemgetter(1),g)) for k,g in groupby(enumerate(specificRegions), lambda x: x[0]-x[1]) ]
			
			for exon in exons:

				exonInfo = {}
				exonInfo["switch"] = "{0}_{1}_{2}".format(gene,switchDict["nIso"],switchDict["tIso"])
				exonInfo["length"] = len(exon)
				exonInfo["origin"] = origin

				if exon & specificCds:
					exonicCds = sorted(list(exon & specificCds))
					
					exonInfo["cdsLength"] = len(exonicCds)
					exonInfo["keepORF"] = True if len(exonicCds)%3==0 else False
					
					firstPos = exonicCds[0]
					pos = [ i for i,x in enumerate(cds) if x==firstPos ][0]
					exonInfo["position"] = float(pos)/len(cds)
					exonInfo["relativeSize"] = float(exonInfo["cdsLength"])/len(cds)

					if exon & specificUtr:
						exonInfo["role"] = "CDS-UTR"
					else:
						exonInfo["role"] = "CDS"

				elif exon & specificUtr:
					exonInfo["role"] = "UTR"
					exonInfo["cdsLength"] = 0
					exonInfo["keepORF"] = "NA"
					exonInfo["position"] = "NA"
					exonInfo["relativeSize"] = "NA"

				self.exonStats[randomTag].append(exonInfo)

		txCorresp = thisSwitch.analyzeSplicing()
		orfChange = 0

		for i in range(len(txCorresp)):
			nVersion = txCorresp[i][0]
			tVersion = txCorresp[i][1]

			if nVersion == tVersion:
				tag = "COMMON"
			elif i == 0:
				tag = "BEGINNING"
			elif i == len(txCorresp) - 1:
				tag = "ENDING"
			else:
				tag = "MIDDLE"
			
			exon = {}

			exon["gene"] = gene
			exon["symbol"] = info["symbol"]
			exon["nTranscript"] = thisSwitch.nTx
			exon["tTranscript"] = thisSwitch.tTx
			exon["nVersion"] = 0 if nVersion is None else len(nVersion)
			exon["tVersion"] = 0 if tVersion is None else len(tVersion)
			exon["tag"] = tag

			orfChange = orfChange + exon["nVersion"]%3 - exon["tVersion"]%3
			if abs(orfChange) >= 3:
				orfChange = orfChange - 3*np.sign(orfChange)
			exon["orfChange"] = orfChange

			self.alternativeSplicingStats[randomTag].append(exon)
		
	def printSwitchInfo(self):

		with open("{0}result_summary/isoform_length{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tRandom\tnIsoLength\ttIsoLength\tnIsoSpecificLength\ttIsoSpecificLength\n")
			for nlen,tlen,nsplen,tsplen in self.proteinStats["NonRandom"]:
				F.write("{0}\tNonRandom\t".format(options.Options().tag))
				F.write("{0}\t{1}\t{2}\t{3}\n".format(nlen,tlen,nsplen,tsplen))

			for nlen,tlen,nsplen,tsplen in self.proteinStats["Random"]:
				F.write("{0}\tRandom\t".format(options.Options().tag))
				F.write("{0}\t{1}\t{2}\t{3}\n".format(nlen,tlen,nsplen,tsplen))
		
		with open("{0}result_summary/switches{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tAnalysis\tBoth\tOnly_nIso\tOnly_tIso\tNone\t")
			F.write("Random_Both\tRandom_Only_nIso\tRandom_Only_tIso\tRandom_None\n")
			F.write("{}\tCDS_study\t".format(options.Options().tag ))
			F.write("{}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["Both"]))
			F.write("{}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["OnlyN"]))
			F.write("{}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["OnlyT"]))
			F.write("{}\t".format(self.switchStats["cdsPresence"]["NonRandom"]["None"]))
			F.write("{}\t".format(self.switchStats["cdsPresence"]["Random"]["Both"]))
			F.write("{}\t".format(self.switchStats["cdsPresence"]["Random"]["OnlyN"]))
			F.write("{}\t".format(self.switchStats["cdsPresence"]["Random"]["OnlyT"]))
			F.write("{}\n".format(self.switchStats["cdsPresence"]["Random"]["None"]))

			F.write("Cancer\tAnalysis\tYes\tNo\tRandom_Yes\tRandom_No\tp\tOR\n")

			F.write("{}\tCDS_change\t".format(options.Options().tag))
			F.write("{}\t".format(self.switchStats["cdsTest"]["NonRandom"]["Change"]))
			F.write("{}\t".format(self.switchStats["cdsTest"]["NonRandom"]["NoChange"]))
			F.write("{}\t".format(self.switchStats["cdsTest"]["Random"]["Change"]))
			F.write("{}\t".format(self.switchStats["cdsTest"]["Random"]["NoChange"]))

			lContingencyTable = [[self.switchStats["cdsTest"]["NonRandom"]["Change"],
								  self.switchStats["cdsTest"]["NonRandom"]["NoChange"]],
								 [self.switchStats["cdsTest"]["Random"]["Change"],
								  self.switchStats["cdsTest"]["Random"]["NoChange"]]]
			OR,pval = fisher_exact(lContingencyTable)

			F.write("{}\t{}\n".format(pval,OR) )

			F.write("{}\tUTR_change\t".format(options.Options().tag ))
			F.write("{}\t".format(self.switchStats["utrTest"]["NonRandom"]["Change"]))
			F.write("{}\t".format(self.switchStats["utrTest"]["NonRandom"]["NoChange"]))
			F.write("{}\t".format(self.switchStats["utrTest"]["Random"]["Change"]))
			F.write("{}\t".format(self.switchStats["utrTest"]["Random"]["NoChange"]))

			lContingencyTable = [[self.switchStats["utrTest"]["NonRandom"]["Change"],
								  self.switchStats["utrTest"]["NonRandom"]["NoChange"]],
								 [self.switchStats["utrTest"]["Random"]["Change"],
								  self.switchStats["utrTest"]["Random"]["NoChange"]]]
			OR,pval = fisher_exact(lContingencyTable)

			F.write("{}\t{}\n".format(pval,OR) )

			##### GENE LEVEL #####
			F.write("Cancer\tAnalysis\tFeat-Switch\tFeat-NoSwitch\tNoFeat-Switch\tNoFeat-NoSwitch\tp\tOR\n")
			# driver enrichment in switches
			F.write("{}\td0_enrichment\t".format(options.Options().tag ))
			F.write("{}\t".format(self.switchStats["d0Enrichment"]["driver"]["Switch"]))
			F.write("{}\t".format(self.switchStats["d0Enrichment"]["driver"]["NoSwitch"]))
			F.write("{}\t".format(self.switchStats["d0Enrichment"]["nonDriver"]["Switch"]))
			F.write("{}\t".format(self.switchStats["d0Enrichment"]["nonDriver"]["NoSwitch"]))
			
			lContingencyTable = [[self.switchStats["d0Enrichment"]["driver"]["Switch"],
								  self.switchStats["d0Enrichment"]["driver"]["NoSwitch"]],
								[self.switchStats["d0Enrichment"]["nonDriver"]["Switch"],
								self.switchStats["d0Enrichment"]["nonDriver"]["NoSwitch"]]]
			OR,pval = fisher_exact(lContingencyTable)

			F.write("{}\t{}\n".format(pval,OR) )

			F.write("{0}\td1_enrichment\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["d1Enrichment"]["D1Driver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["d1Enrichment"]["D1Driver"]["NoSwitch"]))
			F.write("{0}\t".format(self.switchStats["d1Enrichment"]["NonD1Driver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["d1Enrichment"]["NonD1Driver"]["NoSwitch"]))
			
			lContingencyTable = [[self.switchStats["d1Enrichment"]["D1Driver"]["Switch"],
								  self.switchStats["d1Enrichment"]["D1Driver"]["NoSwitch"]],
								[self.switchStats["d1Enrichment"]["NonD1Driver"]["Switch"],
								self.switchStats["d1Enrichment"]["NonD1Driver"]["NoSwitch"]]]
			OR,pval = fisher_exact(lContingencyTable)

			F.write("{}\t{}\n".format(pval,OR) )

			# driver enrichment in functional
			F.write("{0}\td0_functional_enrichment\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["d0Functional"]["driver"]["Functional"]) )
			F.write("{0}\t".format(self.switchStats["d0Functional"]["driver"]["NonFunctional"]) )
			F.write("{0}\t".format(self.switchStats["d0Functional"]["nonDriver"]["Functional"]) )
			F.write("{0}\t".format(self.switchStats["d0Functional"]["nonDriver"]["NonFunctional"]) )
			
			lContingencyTable = [[self.switchStats["d0Functional"]["driver"]["Functional"],
								  self.switchStats["d0Functional"]["driver"]["NonFunctional"]],
								[self.switchStats["d0Functional"]["nonDriver"]["Functional"],
								self.switchStats["d0Functional"]["nonDriver"]["NonFunctional"]]]
			OR,pval = fisher_exact(lContingencyTable)

			F.write("{}\t{}\n".format(pval,OR) )

			# driver d1 enrichment in functional
			F.write("{0}\td1_functional_enrichment\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["d1Functional"]["D1Driver"]["Functional"]) )
			F.write("{0}\t".format(self.switchStats["d1Functional"]["D1Driver"]["NonFunctional"]) )
			F.write("{0}\t".format(self.switchStats["d1Functional"]["NonD1Driver"]["Functional"]) )
			F.write("{0}\t".format(self.switchStats["d1Functional"]["NonD1Driver"]["NonFunctional"]) )
			
			lContingencyTable = [[self.switchStats["d1Functional"]["D1Driver"]["Functional"],
								  self.switchStats["d1Functional"]["D1Driver"]["NonFunctional"]],
								[self.switchStats["d1Functional"]["NonD1Driver"]["Functional"],
								self.switchStats["d1Functional"]["NonD1Driver"]["NonFunctional"]]]
			OR,pval = fisher_exact(lContingencyTable)

			F.write("{}\t{}\n".format(pval,OR) )

			F.write("Cancer\tAnalysis\tMedian case\tMedian control\tp\n")

			# patient difference in switches
			m1 = np.median(np.array(self.switchStats["d0Patients"]["driver"]))
			m2 = np.median(np.array(self.switchStats["d0Patients"]["nonDriver"]))
			F.write("{0}\tDriver_D0_patients\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(m1,m2))
			
			z,p = stats.ranksums(self.switchStats["d0Patients"]["driver"],
								  self.switchStats["d0Patients"]["nonDriver"])

			F.write("{0}\n".format(p) )

			m1 = np.median(np.array(self.switchStats["d0PatientsFunctional"]["driver"]))
			m2 = np.median(np.array(self.switchStats["d0PatientsFunctional"]["nonDriver"]))
			F.write("{0}\tDriver_D0_patients_functional\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(m1,m2))
			
			z,p = stats.ranksums(self.switchStats["d0PatientsFunctional"]["driver"],
								  self.switchStats["d0PatientsFunctional"]["nonDriver"])

			F.write("{0}\n".format(p) )

			m1 = np.median(np.array(self.switchStats["d1Patients"]["D1Driver"]))
			m2 = np.median(np.array(self.switchStats["d1Patients"]["NonD1Driver"]))

			F.write("{0}\tDriver_D1_patients\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(m1,m2))
			
			z,p = stats.ranksums(self.switchStats["d1Patients"]["D1Driver"],
								  self.switchStats["d1Patients"]["NonD1Driver"])

			F.write("{0}\n".format(p) )

			m1 = np.median(np.array(self.switchStats["d1PatientsFunctional"]["D1Driver"]))
			m2 = np.median(np.array(self.switchStats["d1PatientsFunctional"]["NonD1Driver"]))

			F.write("{0}\tDriver_D1_patients_functional\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(m1,m2))
			
			z,p = stats.ranksums(self.switchStats["d1PatientsFunctional"]["D1Driver"],
								  self.switchStats["d1PatientsFunctional"]["NonD1Driver"])

			F.write("{0}\n".format(p) )

		with open("{0}result_summary/exons{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tRandom\tSwitch\tOrigin\tType\tLength\tCDSLength\t")
			F.write("CDSRelativeSize\tPosition\tKeepOrf\n")
			for random in self.exonStats:
				for exon in self.exonStats[random]:
					F.write("{0}\t{1}\t".format(options.Options().tag,random))
					F.write("{0}\t{1}\t".format(exon["switch"],exon["role"]))
					F.write("{0}\t".format(exon["origin"]))
					F.write("{0}\t{1}\t".format(exon["length"],exon["cdsLength"]))
					F.write("{0}\t{1}\t".format(exon["relativeSize"],exon["position"]))
					F.write("{0}\n".format(exon["keepORF"]))

		with open("{0}result_summary/exons_new{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tRandom\tGene\tSymbol\tnTranscript\ttTranscript\tTag\tOrfChange\tnormalSegment\ttumorSegment\n");
			for random in self.alternativeSplicingStats:
				for exon in self.alternativeSplicingStats[random]:
					F.write("{0}\t{1}\t".format(options.Options().tag,random))
					F.write("{0}\t{1}\t".format(exon["gene"],exon["symbol"]))
					F.write("{0}\t{1}\t".format(exon["nTranscript"],exon["tTranscript"]))
					F.write("{0}\t{1}\t".format(exon["tag"],exon["orfChange"]))
					F.write("{0}\t{1}\n".format(exon["nVersion"],exon["tVersion"]))

	def printStructutalInfo(self):

		# with open("{0}result_summary/structural_summary{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
		# 	F.write("Cancer\tGene\tSymbol\tnTx\ttTx\tPfam\t")
		# 	F.write("IUPREDLong\tIUPREDShort\tFunctional\tModel\t")
		# 	F.write("Noise\tDriver\tASDriver\tDriverType\n")
		# 	for tag,random,featureDict in self.featuresTable:
		# 		switchElements = tag.split("_")
		# 		gene = switchElements[0]
		# 		symbol = switchElements[1]
		# 		nTx = switchElements[2]
		# 		tTx = switchElements[3]

		# 		F.write("{0}\t{1}\t".format(options.Options().tag,gene))
		# 		F.write("{0}\t{1}\t{2}\t".format(symbol,nTx,tTx))
		# 		F.write("{0}\t{1}\t".format(len(featureDict["Pfam"]),len(featureDict["IUPREDLong"])))
		# 		F.write("{0}\t".format(len(featureDict["IUPREDShort"])))
		# 		F.write("{0}\t{1}\t".format(featureDict["Functional"],featureDict["Model"]))
		# 		F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["driver"]))
		# 		F.write("{0}\t{1}\n".format(featureDict["asDriver"],featureDict["driverType"]))

		with open("{0}result_summary/structural_features{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tGene\tSymbol\tnTx\ttTx\tRandom\tAnalysis\tWhatsHappenning\t")
			F.write("Feature\tDriver\tASDriver\tDriverType\n")
			
			for tag,randomTag,featureDict in self.featuresTable:
				switchElements = tag.split("_")
				gene = switchElements[0]
				symbol = switchElements[1]
				nTx = switchElements[2]
				tTx = switchElements[3]

				for analysis in ["Pfam","iupred","anchor","prosite"]:
					for data in featureDict[analysis]:
						F.write("{0}\t{1}\t".format(options.Options().tag,gene))
						F.write("{0}\t{1}\t".format(symbol,nTx))
						F.write("{0}\t{1}\t{2}\t".format(tTx,randomTag,analysis))
						F.write("{0}\t{1}\t".format(data[1],data[0].replace(" ","_")))
						F.write("{0}\t{1}\t".format(featureDict["driver"],featureDict["asDriver"]))
						F.write("{0}\n".format(featureDict["driverType"]))

		# with open("{0}result_summary/structural_loops{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
		# 	F.write("Cancer\tDifferent\t")
		# 	F.write("Same\tOnly_nIso\tOnly_tIso\tNone\n")
		# 	F.write("{0}\t".format(options.Options().tag ))
		# 	F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["different"],self.loops["loopsChange"]["same"]) )
		# 	F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["onlyN"],self.loops["loopsChange"]["onlyT"]) )
		# 	F.write("{0}\t{1}\n".format(self.loops["loopsChange"]["noLoops"]))

	def changedStructuralFeatures(self,random,gene,info,switchDict,thisSwitch):	

		tag = "{0}_{1}_{2}_{3}".format(gene,info["symbol"],thisSwitch.nTx,thisSwitch.tTx)

		if random:
			randomTag = "Random"
			filetag = "_random"
		else:
			randomTag = "NonRandom"
			filetag = ""

		switchFeatures = {}

		pfam = []
		prosite = []
		disorder = []
		anchor = []

		if switchDict["functional"]:
			if thisSwitch.domainChange:
				for element in utils.readTable("{}structural_analysis/interpro_analysis{}.tsv".format(options.Options().qout,filetag)):
					if element[2]==thisSwitch.nTx and element[3]==thisSwitch.tTx:
						pfam.append((element[5],element[4]))

			if thisSwitch.disorderChange:
				for element in utils.readTable("{}structural_analysis/iupred_analysis{}.tsv".format(options.Options().qout,filetag)):
					if element[2]==thisSwitch.nTx and element[3]==thisSwitch.tTx:
						if float(element[-1]):
							disorder.append((element[5],element[4]))

			if thisSwitch.anchorChange:
				for element in utils.readTable("{}structural_analysis/anchor_analysis{}.tsv".format(options.Options().qout,filetag)):
					if element[2]==thisSwitch.nTx and element[3]==thisSwitch.tTx:
						if float(element[-1]):
							anchor.append((element[5],element[4]))

			if thisSwitch.ptmChange:
				for element in utils.readTable("{}structural_analysis/prosite_analysis{}.tsv".format(options.Options().qout,filetag)):
					if element[2]==thisSwitch.nTx and element[3]==thisSwitch.tTx:
						prosite.append((element[5],element[4]))
			
		switchFeatures["Pfam"] = pfam
		switchFeatures["iupred"] = disorder
		switchFeatures["anchor"] = anchor
		switchFeatures["prosite"] = prosite
		
		switchFeatures["driver"] = int(info["driver"])
		switchFeatures["asDriver"] = int(info["asDriver"])
		switchFeatures["driverType"] = info["driverType"]
		switchFeatures["Functional"] = int(thisSwitch.is_functional)
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

	def printFunctionalGene(self):

		affectedGenes = {}

		patientsNGenes = [ [x["symbol"],z.patients] for w,x,y,z in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network) ]

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

		with open("{0}result_summary/functionalSwitches.tsv".format(options.Options().qout), "w" ) as F:
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
		[ neighborhoodAnalysis.extend(z.neighborhoodChange) for w,x,y,z in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network) ]
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
			[ geneSets.extend(x["neighborhoods"][analysis]) for w,x,y,z in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network) ]
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
