#!/soft/devel/python-2.7/bin/python

from libs import options
from libs import utils
from methods import method

import cPickle
import fisher
from scipy import stats
from itertools import groupby
from operator import itemgetter

import pdb

class ResultSummary(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self._random_gene_network = cPickle.load(open("{0}randomGeneNetwork.pkl".format(options.Options().qout)))
		self._random_gene_network.createLogger()

		self.proteinStats = { "Random": { "nIsoLength": [], "tIsoLength": [] },
							  "NonRandom": { "nIsoLength": [], "tIsoLength": [] }}

		self.switchStats = {}
		
		# driver/driver+relevance tests
		self.switchStats["driverD0Test"] = { "Driver": {"NoSwitch": 0, "Switch": 0}, 
										   	 "NonDriver": {"NoSwitch": 0, "Switch": 0} }
		self.switchStats["driverD0Patients"] = { "Driver": [], "NonDriver": [] }
		self.switchStats["driverD1Test"] = { "D1Driver": {"NoSwitch": 0, "Switch": 0}, 
										   	 "NonD1Driver": {"NoSwitch": 0, "Switch": 0} }
		self.switchStats["driverD1Patients"] = { "D1Driver": [], "NonD1Driver": [] }
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

		testedSwitches = 0

		# tests at switch level
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,only_models=True,partialCreation=True,removeNoise=True):
			self.logger.debug("Getting statistics for switch {0}_{1}_{2}.".format(gene,thisSwitch.nTx,thisSwitch.tTx))

			# general protein, switch and gene info
			self.switchAndExonOverview(False,gene,info,switchDict,thisSwitch)
			self.proteinOverview(False,switchDict,txDict)

			# structural info
			self.changedStructuralFeatures(False,gene,info,switchDict,thisSwitch)
			# self.loopChange(switch)
			# self.disorderChange(switch)

			testedSwitches += 1
			
		for gene,info,switchDict,thisSwitch in self._random_gene_network.sampleSwitches(self._transcript_network,numIterations=testedSwitches):
			self.switchAndExonOverview(True,gene,info,switchDict,thisSwitch)
			self.proteinOverview(True,switchDict,txDict)

			# self.changedStructuralFeatures(True,gene,info,switchDict,thisSwitch)

		# tests at gene level
		for gene,info in self._gene_network.iterate_genes_ScoreWise():
			dicts = [ x for x in info["isoformSwitches"] if x["model"] and not x["noise"] ]

			switchDict = None
			thisSwitch = None
			if dicts:
				switchDict = dicts[0]
				thisSwitch = self._gene_network.createSwitch(switchDict,self._transcript_network,True)

			self.driverTests(gene,info,switchDict,thisSwitch)

		self.printSwitchInfo()
		self.printStructutalInfo()
		
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
			
		self.switchStats["driverD0Test"][driverTag][switchTag] += 1

		if not info["Driver"]:
			d1 = [ x for x in self._gene_network._net.neighbors(gene) if self._gene_network._net.node[x]["Driver"] ]
			if d1: 	
				d1DriverTag="D1Driver"
			else: 	
				d1DriverTag="NonD1Driver"

			self.switchStats["driverD1Test"][d1DriverTag][switchTag] += 1
		
		if switchDict:
			self.switchStats["driverD0Patients"][driverTag].append(len(switchDict["patients"]))
			if not info["Driver"]:
				self.switchStats["driverD1Patients"][d1DriverTag].append(len(switchDict["patients"]))
		if switchDict is not None:
			self.switchStats["driverRTest"][driverTag][relevanceTag] += 1

	def switchAndExonOverview(self,random,gene,info,switchDict,thisSwitch):
		nTx = thisSwitch.nTranscript
		tTx = thisSwitch.tTranscript

		if random: 	randomTag = "Random"
		else:		randomTag = "NonRandom"

		if thisSwitch.is_relevant:	relevanceTag = "Relevant"
		else:						relevanceTag = "NonRelevant"

		if thisSwitch.cds_diff: cdsTag = "Change"
		else:					cdsTag = "NoChange"

		if thisSwitch.utr_diff: utrTag = "Change"
		else: 					utrTag = "NoChange"

		if nTx.cds and tTx.cds: presenceTag = "Both"
		elif nTx.cds: 			  presenceTag = "OnlyN"
		elif tTx.cds: 			  presenceTag = "OnlyT"
		else: 					  presenceTag = "None"

		self.switchStats["relevanceTest"][randomTag][relevanceTag] += 1
		self.switchStats["cdsTest"][randomTag][cdsTag] += 1
		self.switchStats["utrTest"][randomTag][utrTag] += 1
		self.switchStats["cdsPresence"][randomTag][presenceTag] += 1

		nSpecificCds = set([ x for x in nTx.cds if x not in tTx.cds ])
		nSpecificUtr = set([ x for x in nTx.utr if x not in tTx.utr ])
		tSpecificCds = set([ x for x in tTx.cds if x not in nTx.cds ])
		tSpecificUtr = set([ x for x in tTx.utr if x not in nTx.utr ])
		
		for specificCds,specificUtr,cds in zip([nSpecificCds,tSpecificCds],[nSpecificUtr,tSpecificUtr],[nTx.cds,tTx.cds]):
		
			specificRegions = sorted(list(specificUtr | specificCds))
			exons = [ set(map(itemgetter(1),g)) for k,g in groupby(enumerate(specificRegions), lambda (i,x): i-x) ]
			
			for exon in exons:

				exonInfo = {}
				exonInfo["switch"] = "{0}_{1}_{2}".format(gene,switchDict["nIso"],switchDict["tIso"])
				exonInfo["length"] = len(exon)

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
					exonInfo["keepORF"] = None
					exonInfo["position"] = None
					exonInfo["relativeSize"] = None

				self.exonStats[randomTag].append(exonInfo)
		
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
			F.write("{0}\tDriver_D0_enrichment\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["driverD0Test"]["Driver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["driverD0Test"]["Driver"]["NoSwitch"]))
			F.write("{0}\t".format(self.switchStats["driverD0Test"]["NonDriver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["driverD0Test"]["NonDriver"]["NoSwitch"]))
			
			p = fisher.pvalue(self.switchStats["driverD0Test"]["Driver"]["Switch"],
							self.switchStats["driverD0Test"]["Driver"]["NoSwitch"],
							self.switchStats["driverD0Test"]["NonDriver"]["Switch"],
							self.switchStats["driverD0Test"]["NonDriver"]["NoSwitch"])

			F.write("{0}\n".format(p.two_tail) )

			F.write("{0}\tDriver_D1_enrichment\t".format(options.Options().tag ))
			F.write("{0}\t".format(self.switchStats["driverD1Test"]["D1Driver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["driverD1Test"]["D1Driver"]["NoSwitch"]))
			F.write("{0}\t".format(self.switchStats["driverD1Test"]["NonD1Driver"]["Switch"]))
			F.write("{0}\t".format(self.switchStats["driverD1Test"]["NonD1Driver"]["NoSwitch"]))
			
			p = fisher.pvalue(self.switchStats["driverD1Test"]["D1Driver"]["Switch"],
							self.switchStats["driverD1Test"]["D1Driver"]["NoSwitch"],
							self.switchStats["driverD1Test"]["NonD1Driver"]["Switch"],
							self.switchStats["driverD1Test"]["NonD1Driver"]["NoSwitch"])

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

			F.write("Cancer\tAnalysis\tMean cond1\tMean cond2\tp\n")
			# patient difference in switches
			F.write("{0}\tDriver_D0_patients\t".format(options.Options().tag ))
			F.write("{0}\t".format(sum(self.switchStats["driverD0Patients"]["Driver"])/len(self.switchStats["driverD0Patients"]["Driver"])))
			F.write("{0}\t".format(sum(self.switchStats["driverD0Patients"]["NonDriver"])/len(self.switchStats["driverD0Patients"]["NonDriver"])))
			
			t,p = stats.ttest_ind(self.switchStats["driverD0Patients"]["Driver"],
								  self.switchStats["driverD0Patients"]["NonDriver"])

			F.write("{0}\n".format(p) )

			F.write("{0}\tDriver_D1_patients\t".format(options.Options().tag ))
			F.write("{0}\t".format(sum(self.switchStats["driverD1Patients"]["D1Driver"])/len(self.switchStats["driverD1Patients"]["D1Driver"])))
			F.write("{0}\t".format(sum(self.switchStats["driverD1Patients"]["NonD1Driver"])/len(self.switchStats["driverD1Patients"]["NonD1Driver"])))
			
			t,p = stats.ttest_ind(self.switchStats["driverD1Patients"]["D1Driver"],
								  self.switchStats["driverD1Patients"]["NonD1Driver"])

			F.write("{0}\n".format(p) )

		with open("{0}result_summary/exons{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tRandom\tSwitch\tType\tLength\tCDSLength\t")
			F.write("CDSRelativeSize\tPosition\tKeepOrf\tPosition\n")
			for random in self.exonStats:
				for chunk in self.exonStats[random]:
					F.write("{0}\t{1}\t".format(options.Options().tag,random))
					F.write("{0}\t{1}\t".format(exon["switch"],exon["role"]))
					F.write("{0}\t{1}\t".format(exon["length"],exon["cdsLength"]))
					F.write("{0}\t{1}\t".format(exonInfo["relativeSize"],exon["position"]))
					F.write("{0}\t{1}\n".format(exon["keepORF"],exon["position"]))

	def printStructutalInfo(self):

		with open("{0}result_summary/structural_features{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tGene\tSymbol\tnTx\ttTx\tAnalysis\tRandom\tWhatsHappenning\t")
			F.write("Feature\tRelevant\tModel\tNoise\tDriver\tASDriver\tDriverType\n")
			
			for s in self.featuresTable:
				switchElements = s[0].split("_")
				gene = switchElements[0]
				symbol = switchElements[1]
				nTx = switchElements[2]
				tTx = switchElements[3]

				randomTag = s[1]
				featureDict = s[2]

				for pfamDom in featureDict["Pfam"]:
					F.write("{0}\t{1}\t{2}\t".format(options.Options().tag,gene,symbol))
					F.write("{0}\t{1}\tPfam\t{2}\t".format(nTx,tTx,randomTag))
					F.write("{0}\t{1}\t".format(pfamDom[1],pfamDom[0].replace(" ","_")))
					F.write("{0}\t{1}\t".format(featureDict["Relevant"],featureDict["Model"]))
					F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["Driver"]))
					F.write("{0}\t{1}\n".format(featureDict["ASDriver"],featureDict["DriverType"]))

				# if featureDict["IUPREDLong"]:
				# 	for iulong in featureDict["IUPREDLong"][0]:
				# 		F.write("{0}\t{1}\tIUPREDLong\t".format(options.Options().tag,gene))
				# 		F.write("{0}\t{1}\t{2}\t".format(symbol,nTx,tTx))

				# 		F.write("{0}\t{1}\t".format(featureDict["IUPREDLong"][1],iulong))
				# 		F.write("{0}\t{1}\t".format(featureDict["Relevant"],featureDict["Model"]))
				# 		F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["Driver"]))
				# 		F.write("{0}\t{1}\n".format(featureDict["ASDriver"],featureDict["DriverType"]))

				# if featureDict["IUPREDShort"]:
				# 	for iushort in featureDict["IUPREDShort"][0]:
				# 		F.write("{0}\t{1}\tIUPREDShort\t".format(options.Options().tag,gene))
				# 		F.write("{0}\t{1}\t{2}\t".format(symbol,nTx,tTx))
				# 		F.write("{0}\t{1}\t".format(featureDict["IUPREDShort"][1],iushort))
				# 		F.write("{0}\t{1}\t".format(featureDict["Relevant"],featureDict["Model"]))
				# 		F.write("{0}\t{1}\t".format(featureDict["Noise"],featureDict["Driver"]))
				# 		F.write("{0}\t{1}\n".format(featureDict["ASDriver"],featureDict["DriverType"]))

		with open("{0}result_summary/structural_loops{1}.tsv".format(options.Options().qout,options.Options().filetag), "w" ) as F:
			F.write("Cancer\tDifferent\t")
			F.write("Same\tOnly_nIso\tOnly_tIso\tNone\n")
			F.write("{0}\t".format(options.Options().tag ))
			F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["different"],self.loops["loopsChange"]["same"]) )
			F.write("{0}\t{1}\t".format(self.loops["loopsChange"]["onlyN"],self.loops["loopsChange"]["onlyT"]) )
			F.write("{0}\t{1}\n".format(self.loops["loopsChange"]["noLoops"]))

	def changedStructuralFeatures(self,random,gene,info,switchDict,switch):	

		tag = "{0}_{1}_{2}_{3}".format(gene,info["symbol"],switch.nTx,switch.tTx)

		if random:
			randomTag = "Random"
			filetag = "_random"
		else:
			randomTag = "NonRandom"
			filetag = ""

		switchFeatures = {}

		Pfam = []

		functionalChanges = ""
		for element in utils.readTable("{0}structural_analysis/interpro_analysis{1}.tsv".format(options.Options().qout,filetag)):
			if element[2]==switch.nTx and element[3]==switch.tTx:
				if element[5]=="Pfam":
					Pfam.append((element[7],element[4]))

		switchFeatures["Pfam"] = Pfam
		
		# if switch.disorderChange:
		# 	longDisorder = []
		# 	shortDisorder = []

		# 	for element in utils.readTable("{0}structural_analysis/iupred_analysis.tsv".format(options.Options().qout)):
		# 		if element[0]==gene and element[2] in [switch.nTx,switch.tTx] and element[4]:
		# 			if element[3]=='short':
		# 				if element[2] == switch.nTx:
		# 					shortDisorder = (element[4].split(","),'Lost in tumor')
		# 				else:
		# 					shortDisorder = (element[4].split(","),'Gained in tumor')
		# 			elif element[3]=='long':
		# 				if element[2] == switch.nTx:
		# 					longDisorder = (element[4].split(","),'Lost in tumor')
		# 				else:
		# 					longDisorder = (element[4].split(","),'Gained in tumor')
			
		# 	switchFeatures["IUPREDShort"] = shortDisorder
		# 	switchFeatures["IUPREDLong"] = longDisorder
		# else:
		# 	switchFeatures["IUPREDShort"] = []
		# 	switchFeatures["IUPREDLong"] = []
		
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
