#!/soft/devel/python-2.7/bin/python

from interface import interpro_analysis
from libs import options
from libs import utils
from methods import method

from collections import Counter
import fisher
import glob
import numpy as np
import os

class StructuralAnalysis(method.Method):
	def __init__(self,gn_network,tx_network,isRand=False):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self.isRandom = isRand
		if self.isRandom:
			#self.anchor_threshold = 1
			#self.iupred_threshold = 1
			#self.prosite_threshold = 1
			tag = "_random"
			if options.Options().parallelRange:
				tag += "_{}".format(options.Options().parallelRange)
		else:
			tag = ""
			if not options.Options().parallelRange:
				#self.joinFiles("_random")
				self.newJoinFiles("_random")
				
			else:
				tag = "_{}".format(options.Options().parallelRange)
			#self.anchor_threshold = self.getThreshold("anchor",80,tag="")
			#self.iupred_threshold = self.getThreshold("iupred",80,tag="")
			#self.prosite_threshold = self.getThreshold("prosite",80,tag="")

		self.anchor_threshold = 0.5
		self.iupred_threshold = 0.5

		self._interpro_file = "{0}structural_analysis/interpro_analysis{1}.tsv".format(options.Options().qout,tag)
		self.IP = open(self._interpro_file,"w")
		self.writeIPHeader(self.IP)
				
		self._iupred_file = "{0}structural_analysis/iupred_analysis{1}.tsv".format(options.Options().qout,tag)
		self.IU = open(self._iupred_file,"w")
		self.writeIUHeader(self.IU)
								
		self._anchor_file = "{0}structural_analysis/anchor_analysis{1}.tsv".format(options.Options().qout,tag)
		self.ANCHOR = open(self._anchor_file,"w")
		self.writeANCHORHeader(self.ANCHOR)
				
		self._prosite_file = "{0}structural_analysis/prosite_analysis{1}.tsv".format(options.Options().qout,tag)
		self.PROSITE = open(self._prosite_file,"w")
		self.newWritePrositeHeader(self.PROSITE)
						
		self._relevance_info = "{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,tag)
		self.REL = open(self._relevance_info,"w")
		self.writeSummaryHeader(self.REL)
				
	def run(self):		
		self.logger.info("Structural analysis.")

		# get loops families

		isoInfo = {}
		for line in open(options.Options().wd+"Data/TCGA/UnifiedFasta_iLoops13.fa"):
			if ">" in line:
				elements = line.strip().split("#")
				isoInfo[elements[0][1:]] = {}
				isoInfo[elements[0][1:]]["UniProt"] = elements[2]
				isoInfo[elements[0][1:]]["iLoopsFamily"] = elements[3]

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,partialCreation=True,removeNoise=True):
			
			if not thisSwitch.nIsoform or not thisSwitch.tIsoform: 
				continue

			thisSwitch._iloops_change 	  = self.archDBAnalysis(thisSwitch,gene,info,isoInfo)
			thisSwitch._functional_change = self.interproAnalysis(thisSwitch,gene,info)
			thisSwitch._disorder_change   = self.disorderAnalysis(thisSwitch,gene,info)
			thisSwitch._anchor_change     = self.anchorAnalysis(thisSwitch,gene,info)
			
			#thisSwitch._ptm_change 		  = self.prositeAnalysis(thisSwitch,gene,info)
			thisSwitch._ptm_change 		  = self.newPrositeAnalysis(thisSwitch,gene,info)

			self.REL.write("{0}\t{1}\t".format(gene,info["symbol"]))
			self.REL.write("{0}\t{1}\t".format(thisSwitch.nTx,thisSwitch.tTx))
			self.REL.write("{0}\t{1}\t".format(thisSwitch._iloops_change,thisSwitch._functional_change))
			self.REL.write("{0}\t{1}\t".format(thisSwitch._disorder_change,thisSwitch._anchor_change))
			self.REL.write("{0}\n".format(thisSwitch._ptm_change))

		self.IP.close()
		self.IU.close()
		self.ANCHOR.close()
		self.PROSITE.close()
		self.REL.close()

	def archDBAnalysis(self,thisSwitch,gene,info,isoInfo):

		self.logger.debug("iLoops: looking for loop changes for gene {0}.".format(gene) )
		
		if thisSwitch.nTx in isoInfo and thisSwitch.tTx in isoInfo:
			nLoops = isoInfo[thisSwitch.nTx]["iLoopsFamily"]
			tLoops = isoInfo[thisSwitch.tTx]["iLoopsFamily"]

			if nLoops != tLoops:
				self.logger.debug("iLoops: information found for gene {0}.".format(gene) )
				return True
			else:
				return False

		return False

	def interproAnalysis(self,thisSwitch,gene,info):

		self.logger.debug("InterPro: searching changes for gene {0}.".format(gene) )
		
		normalProtein = thisSwitch.nIsoform
		tumorProtein = thisSwitch.tIsoform

		for protein in [normalProtein,tumorProtein]:
			if not protein: continue

			interproFile = interpro_analysis.InterproAnalysis().launchAnalysis(protein.tx, protein.seq)
			if not interproFile: continue
			protein.getFeatures(interproFile)

		if not normalProtein or not tumorProtein: return False
		if not normalProtein._features and not tumorProtein._features: return False
		
		nIsoFeatures = Counter([ x["accession"] for x in normalProtein._features ])
		tIsoFeatures = Counter([ x["accession"] for x in tumorProtein._features ])

		nIso_uniqFeats = dict([ (x,nIsoFeatures[x]-tIsoFeatures.get(x,0)) for x in nIsoFeatures if nIsoFeatures[x]-tIsoFeatures.get(x,0) > 0 ])
		tIso_uniqFeats = dict([ (x,tIsoFeatures[x]-nIsoFeatures.get(x,0)) for x in tIsoFeatures if tIsoFeatures[x]-nIsoFeatures.get(x,0) > 0 ])

		iterationZip = zip([nIsoFeatures,tIsoFeatures],[nIso_uniqFeats,tIso_uniqFeats],[normalProtein,tumorProtein],["Lost_in_tumor","Gained_in_tumor"])

		for features,uniqFeatures,protein,whatShouldBeHappening in iterationZip:
			for f in features:
		
				whatsHappening = "Nothing"
				featInfo = [ x for x in protein._features if x["accession"]==f ][0]
				notInUniq = 0
				if f in uniqFeatures:
					whatsHappening = whatShouldBeHappening
					notInUniq = uniqFeatures[f]

				self.IP.write("{0}\t{1}\t{2}\t".format(gene,info["symbol"],thisSwitch.nTx))
				self.IP.write("{0}\t{1}\t".format(thisSwitch.tTx,whatsHappening))
				self.IP.write("{0}\t{1}\t".format(featInfo["analysis"],featInfo["accession"]))
				self.IP.write("{0}\t{1}\t".format(featInfo["description"].replace(" ","_"),features[f]))
				self.IP.write("{0}\n".format(notInUniq))

				toSave = [featInfo["analysis"],featInfo["accession"],featInfo["description"],notInUniq]
				if thisSwitch._domain_change is None:
					thisSwitch._domain_change = []
				thisSwitch._domain_change.append(toSave)

		if thisSwitch._domain_change: 
			self.logger.debug("InterPro: information found for gene {0}.".format(gene))
			return True
		else: return False

	def disorderAnalysis(self,thisSwitch,gene,info):

		self.logger.debug("IUPRED: Searching disorder for gene {0}.".format(gene))
		
		anyIUpredSeq = False
		normalProtein = thisSwitch.nIsoform
		tumorProtein = thisSwitch.tIsoform

		if not normalProtein or not tumorProtein: return False

		for protein,whatShouldBeHappening in zip([normalProtein,tumorProtein],["Lost_in_tumor","Gained_in_tumor"]):
			for mode in ["short","long"]:
				#Parse iupred output
				for line in self.getIupredOutput(protein,mode):
					resNum  = int(line[0])
					residue	= line[1]
					score 	= float(line[-1])

					thisRes = protein._structure[resNum-1]
					if residue != thisRes.res:
						self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(protein.tx))
						continue
					thisRes.set_iuPredScore(score)

				disordered = protein.getSegments("disordered",minLength=5,gap=2)
				isoform = protein.getSegments("isoform-specific")

				for disorderedRegion in disordered:
					whatsHappening = whatShouldBeHappening
					disorderedRegionSet = set(disorderedRegion)

					overlappingIsoSpecific = []

					for isoSpecific in isoform:
						isoSpecificSet = set(isoSpecific)

						if isoSpecificSet & disorderedRegionSet:
							overlappingIsoSpecific.extend(isoSpecific)
					
					if not overlappingIsoSpecific:
						jaccard = 0
						macroScore = 0
						microScore = 0
						whatsHappening = "Nothing"

					else:
						overlappingIsoSpecificSet = set(overlappingIsoSpecific)

						intersection = float(len(overlappingIsoSpecificSet & disorderedRegionSet))
						union = float(len(overlappingIsoSpecificSet | disorderedRegionSet))

						jaccard = intersection/union
						microScore = intersection/len(overlappingIsoSpecificSet)
						macroScore = intersection/len(disorderedRegionSet)

					motifSequence = ""
					start = float("inf")
					end = float("-inf")
					for thisRes in disorderedRegion:
						res = thisRes.res.upper() if thisRes.isoformSpecific else thisRes.res.lower()
						motifSequence += res
						if thisRes.num > end:
							end = thisRes.num
						if thisRes.num < start:
							start = thisRes.num

					significant = max(microScore,macroScore) > self.iupred_threshold

					self.IU.write("{0}\t{1}\t".format(gene,info["symbol"]))
					self.IU.write("{0}\t{1}\t".format(normalProtein.tx,tumorProtein.tx))
					self.IU.write("{0}\t{1}\t".format(whatsHappening,motifSequence))
					self.IU.write("{0}\t{1}\t".format(jaccard,microScore))
					self.IU.write("{0}\t{1}\n".format(macroScore,int(significant)))

					if significant:
						anyIUpredSeq = True

		return anyIUpredSeq

	def anchorAnalysis(self,thisSwitch,gene,info):

		self.logger.debug("ANCHOR: Searching anchoring regions for gene {0}.".format(gene))
		
		anyAnchorSeq = False
		normalProtein = thisSwitch.nIsoform
		tumorProtein = thisSwitch.tIsoform

		if not normalProtein or not tumorProtein: return False

		for protein,whatShouldBeHappening in zip([normalProtein,tumorProtein],["Lost_in_tumor","Gained_in_tumor"]):
			
			#Parse anchor output
			for line in self.getAnchorOutput(protein):
				resNum  = int(line[0])
				residue	= line[1]
				score 	= float(line[2].strip())

				thisRes = protein._structure[resNum-1]

				if residue != thisRes.res:
					self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(protein.tx))
					continue

				thisRes.set_anchorScore(score)

			anchor = protein.getSegments("anchor",minLength=5,gap=2)
			isoform = protein.getSegments("isoform-specific")

			for anchorRegion in anchor:
				whatsHappening = whatShouldBeHappening
				anchorRegionSet = set(anchorRegion)

				overlappingIsoSpecific = []

				for isoSpecific in isoform:
					isoSpecificSet = set(isoSpecific)

					if isoSpecificSet & anchorRegionSet:
						overlappingIsoSpecific.extend(isoSpecific)
				
				if not overlappingIsoSpecific:
					jaccard = 0
					macroScore = 0
					microScore = 0
					whatsHappening = "Nothing"
				else:
					overlappingIsoSpecificSet = set(overlappingIsoSpecific)

					intersection = float(len(overlappingIsoSpecificSet & anchorRegionSet))
					union = float(len(overlappingIsoSpecificSet | anchorRegionSet))

					jaccard = intersection/union
					microScore = intersection/len(overlappingIsoSpecificSet)
					macroScore = intersection/len(anchorRegionSet)

				motifSequence = ""
				start = float("inf")
				end = float("-inf")

				for thisRes in anchorRegion:
					res = thisRes.res.upper() if thisRes.isoformSpecific else thisRes.res.lower()
					motifSequence += res
					if thisRes.num > end:
						end = thisRes.num
					if thisRes.num < start:
						start = thisRes.num

				significant = max(microScore,macroScore) > self.anchor_threshold

				self.ANCHOR.write("{0}\t{1}\t".format(gene,info["symbol"]))
				self.ANCHOR.write("{0}\t{1}\t".format(normalProtein.tx,tumorProtein.tx))
				self.ANCHOR.write("{0}\t{1}\t".format(whatsHappening,motifSequence))
				self.ANCHOR.write("{0}\t{1}\t".format(jaccard,microScore))
				self.ANCHOR.write("{0}\t{1}\n".format(macroScore,int(significant)))

				if significant:
					anyAnchorSeq = True

		return anyAnchorSeq

	def prositeAnalysis(self,thisSwitch,gene,info):
		normalProtein = thisSwitch.nIsoform
		tumorProtein = thisSwitch.tIsoform
		happenings = {}

		anyPTM = False

		for protein,whatShouldBeHappening in zip([normalProtein,tumorProtein],["Lost_in_tumor","Gained_in_tumor"]):
			if not protein: continue

			txFile = "{0}Data/TCGA/ProSite/{1}.out".format(options.Options().wd,protein.tx)
			if os.stat(txFile).st_size == 0:
				continue
			
			for line in utils.readTable(txFile,header=False):
				if ' ' in line[0]:
					newLine = [ line[0].split()[0] ]
					newLine.append('')
					newLine.extend(line[1:4])
					line = newLine
				score = float(line[0])
				start = int(line[2][:-2])
				end = int(line[3])
				motif = line[4]

				for resNum in range(start,end):
					thisRes = protein._structure[resNum-1]
					thisRes._ptms.append(motif)

			ptms = []
			[ ptms.extend(x._ptms) for x in protein._structure ]
			ptms = set(ptms)

			ptmRegions = {}
			for ptm in ptms:
				ptmRegions[ptm] = protein.getSegments(ptm,minLength=1,gap=0)
			isoSpecificRegions = protein.getSegments("isoform-specific")

			for ptm in ptmRegions:
				for region in ptmRegions[ptm]:
					whatsHappening = whatShouldBeHappening

					regionSet = set(region)
					overlappingIsoSpecific = []
					[ overlappingIsoSpecific.extend(x) for x in isoSpecificRegions if set(x) & regionSet ]

					if not overlappingIsoSpecific:
						jaccard = 0
						macroScore = 0
						microScore = 0
						whatsHappening = "Nothing"
					else:
						overlappingIsoSpecificSet = set(overlappingIsoSpecific)

						intersection = float(len(overlappingIsoSpecificSet & regionSet))
						union = float(len(overlappingIsoSpecificSet | regionSet))

						jaccard = intersection/union
						microScore = intersection/len(overlappingIsoSpecificSet)
						macroScore = intersection/len(regionSet)

					self.PROSITE.write("{0}\t{1}\t{2}\t".format(gene,info["symbol"],thisSwitch.nTx))
					self.PROSITE.write("{0}\t{1}\t{2}\t".format(thisSwitch.tTx,whatsHappening,ptm.replace(" ","_")))
					self.PROSITE.write("{0}\t{1}\t".format(jaccard,self.prosite_threshold))
					self.PROSITE.write("{0}\t".format(int(jaccard > self.prosite_threshold)))
					self.PROSITE.write("{0}\t{1}\n".format(microScore,macroScore))

					if jaccard > self.prosite_threshold:
						anyPTM = True

		return anyPTM

	def getThreshold(self,sInputType,iQuantileThreshold,jaccardIndex=[],tag="_random"):

		sThresholdFile = "{0}structural_analysis/{1}_analysis{2}.tsv".format(options.Options().qout,sInputType,tag)
		
		# check if file exists
		if not os.path.isfile(sThresholdFile):
			return 1
		if not jaccardIndex:
			for line in utils.readTable(sThresholdFile):
				j = float(line[6])
				# remove zeros
				if j:
					jaccardIndex.append(j)
		if jaccardIndex:
			threshold = np.percentile(jaccardIndex,iQuantileThreshold)
		else:
			# list is empty, no threshold can be calculated
			return 1

		return threshold

	def recalculateScores(self,sInputType,iQuantileThreshold,valueList=[],tag=""):
		relevantSwitches = []
		outFile = "{0}structural_analysis/{1}_analysis{2}.tsv".format(options.Options().qout,sInputType,tag)
		files = glob.glob("{0}structural_analysis/{1}_analysis{2}_[0-9]*.tsv".format(options.Options().qout,sInputType,tag))

		# close writing files in case they are opened
		if (hasattr(self,'IP') and hasattr(self,'IU') and hasattr(self,'ANCHOR') and hasattr(self,'PROSITE')):
			self.IP.close()
			self.IU.close()
			self.ANCHOR.close()
			self.PROSITE.close()	
			
		if valueList:
			#threshold = self.getThreshold(sInputType,iQuantileThreshold,valueList)
			threshold = 0.5
			with open(outFile,"w") as OUT:
				if sInputType=='interpro': self.writeIPHeader(OUT)
				elif sInputType=='iupred': self.writeIUHeader(OUT)
				elif sInputType=='anchor': self.writeANCHORHeader(OUT)
				elif sInputType=='prosite': self.newWritePrositeHeader(OUT)
					
				for aFile in files:
					for line in utils.readTable(aFile):
						jaccard = float(line[6])

						OUT.write("{0}\t{1}\t".format(line[0],line[1]))
						OUT.write("{0}\t{1}\t".format(line[2],line[3]))
						OUT.write("{0}\t{1}\t".format(line[4],line[5]))
						OUT.write("{0}\t{1}\t".format(line[6],threshold))
						OUT.write("{0}\t".format(int(jaccard > threshold)))
						OUT.write("{0}\t{1}\n".format(line[9],line[10]))

						if jaccard > threshold:
							relevantSwitches.append("{0}_{1}_{2}_{3}".format(line[0],line[1],line[2],line[3]))
		else:
			with open(outFile,"w") as OUT:
				if sInputType=='interpro': self.writeIPHeader(OUT)
				elif sInputType=='iupred': self.writeIUHeader(OUT)
				elif sInputType=='anchor': self.writeANCHORHeader(OUT)
				elif sInputType=='prosite': self.newWritePrositeHeader(OUT)
					
				for aFile in files:
					with open(aFile) as IN:
						IN.readline()
						OUT.write(IN.read())
		
		# reopen files in case they can
		if (hasattr(self,'_interpro_file') and hasattr(self,'_iupred_file') and hasattr(self,'_anchor_file') and hasattr(self,'_prosite_file')):
			self.IP = open(self._interpro_file,"a")
			self.IU = open(self._iupred_file,"a")
			self.ANCHOR = open(self._anchor_file,"a")
			self.PROSITE = open(self._prosite_file,"a")

		return relevantSwitches

	def joinFiles(self,tag=""):

		roots = ['interpro','iupred','anchor','prosite']
		rootDict = {'iupred':[],'anchor':[],'prosite':[]}
		switches = {'iupred':[],'anchor':[],'prosite':[]}

		for root in roots:
			if root in rootDict:
				files = glob.glob("{0}structural_analysis/{1}_analysis{2}_[0-9]*.tsv".format(options.Options().qout,root,tag))
				for aFile in files:
					for line in utils.readTable(aFile):
						j = float(line[6])
						# remove zeroes
						if j:
							rootDict[root].append(j)

				switches[root] = self.recalculateScores(root,80,rootDict[root],tag=tag)
			else:
				self.recalculateScores(root,80,tag=tag)

		self.recalculateRelevance(switches,tag)

	def recalculateRelevance(self,switches,tag):
		files = glob.glob("{0}structural_analysis/structural_summary{1}_[0-9]*.tsv".format(options.Options().qout,tag))
		outFile = "{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,tag)

		with open(outFile,"w") as OUT:
			self.writeSummaryHeader(OUT)
			for aFile in files:
				for line in utils.readTable(aFile):
					switchId = "{0}_{1}_{2}_{3}".format(line[0],line[1],line[2],line[3])

					gene = line[0]
					symbol = line[1]
					ntx = line[2]
					ttx = line[3]
					rLoop = line[4]
					rDomain = line[5]

					if switchId in switches['iupred']:
						rDisorder = "True"
					else:
						rDisorder = "False"
					
					if switchId in switches['anchor']:
						rAnchor = "True"
					else:
						rAnchor = "False"

					if switchId in switches['prosite']:
						rProsite = "True"
					else:
						rProsite = "False"
					
					OUT.write("{0}\t{1}\t".format(gene,symbol))
					OUT.write("{0}\t{1}\t".format(ntx,ttx))
					OUT.write("{0}\t{1}\t".format(rLoop,rDomain))
					OUT.write("{0}\t{1}\t".format(rDisorder,rAnchor))
					OUT.write("{0}\n".format(rProsite))

	def getIupredOutput(self,protein,mode):

		outfile = "{0}Data/{1}/IUPred/{2}.{3}.txt".format(options.Options().wd,options.Options().inputType,protein.tx,mode)
		out = []

		if not os.path.isfile(outfile):
			outfile = "{0}IUPred/{1}.{2}.txt".format(options.Options().wd,protein.tx,mode)
			fFile = "{0}{1}{2}.fa".format(options.Options().qout,protein.tx,options.Options().filetag)
			with open(fFile,"w") as FASTA:
				FASTA.write(">{0}\n{1}\n".format(protein.tx,protein.seq))

			proc = utils.cmdOut(options.Options().wd+"Pipeline/libs/bin/iupred/iupred",fFile,mode)
			out = [ x.strip().split(" ") for x in proc.stdout if "#" not in x ]
			os.remove(fFile)

			with open(outfile,"w") as IUout:
				for line in out:
					IUout.write("{0}\t{1}\t{2}\n".format(line[0],line[1],line[-1]))
					yield line

		else:
			with open(outfile) as IUout:
				for line in IUout:
					yield line.strip().split()

	def getAnchorOutput(self,protein):
		outfile = "{0}Data/{1}/ANCHOR/{2}.txt".format(options.Options().wd,options.Options().inputType,protein.tx)
		out = []

		if not os.path.isfile(outfile):
			outfile = "{0}ANCHOR/{1}.txt".format(options.Options().wd,protein.tx)
			fFile = "{0}{1}{2}.fa".format(options.Options().qout,protein.tx,options.Options().filetag)
			with open(fFile,"w") as FASTA:
				FASTA.write(">{0}\n{1}\n".format(protein.tx,protein.seq))

			proc = utils.cmdOut(options.Options().wd+"Pipeline/libs/ANCHOR/anchor",fFile)
			out = [ x.strip().split() for x in proc.stdout if "#" not in x ]
			os.remove(fFile)

			with open(outfile,"w") as ANCHORout:
				for line in out:
					ANCHORout.write("{0}\t{1}\t{2}\n".format(line[0],line[1],line[2]))
					yield line

		else:
			with open(outfile) as ANCHORout:
				for line in ANCHORout:
					yield line.strip().split()

	def writeIPHeader(self,IP):
		IP.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
		IP.write("What\tAnalysis\tFeature_accesion\tFeature\tRepetitions\t")
		IP.write("Exclusive repetitions\n")

	def writeIUHeader(self,IU):
		IU.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
		IU.write("What\tSequence\tJaccard\tmicroScore\tmacroScore\t")
		IU.write("Significant\n")
		
	def writeANCHORHeader(self,ANCHOR):
		ANCHOR.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
		ANCHOR.write("What\tSequence\tJaccard\tmicroScore\tmacroScore\t")
		ANCHOR.write("Significant\n")

	def writePrositeHeader(self,PROSITE):
		PROSITE.write("Gene\tSymbol\tNormalTranscript\t")
		PROSITE.write("tTumorTranscript\tWhat\tMotif\t")
		PROSITE.write("Jaccard\tThreshold\tSignificant\tmicroScore\tmacroScore\n")

	def newWritePrositeHeader(self,PROSITE):
		PROSITE.write("Gene\tSymbol\tNormalTranscript\t")
		PROSITE.write("tTumorTranscript\tWhat\tFeature\t")
		PROSITE.write("normalReps\ttumorReps\tnMacroScore\t")
		PROSITE.write("nMicroScore\tnJaccard\ttMacroScore\t")
		PROSITE.write("tMicroScore\ttJaccard\n")
		
	def writeSummaryHeader(self,REL):
		REL.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
		REL.write("iLoops\tDomain\tDisorder\tAnchor\tPTM\n")

	'''
	def gpsAnalysis(self,switch,gene,info):
		normalProtein = switch.nIsoform
		tumorProtein = switch.tIsoform
		gpsFile = "{0}Data/TCGA/GPS.out".format(options.Options().wd)
		happenings = {}

		for protein,whatsHappening in zip([normalProtein,tumorProtein],["Lost_in_tumor","Gained_in_tumor"]):
			if not protein: continue
			reading = False
			for line in utils.readTable(gpsFile,header=True):
				if protein.tx not in line[0] and not reading:
					continue
				elif protein.tx in line[0] and not reading:
					reading = True
					continue
				elif ">" in line[0]:
					reading = False
					continue

				resNum 		= int(line[0])
				res 		= line[1]
				kinase 		= line[2]
				motif 		= line[3]
				score 		= float(line[4])
				threshold 	= float(line[5])

				if score < threshold:
					continue

				thisRes = protein._structure[resNum-1]
				thisRes._kinases.add(kinase)

			happenings[whatsHappening] = [ x for x in protein._structure if x._kinases and x.isoformSpecific ]

		anyPTM = False
		for whatsHappening in happenings:
			for thisRes in happenings[whatsHappening]:
				self.GPS.write("{0}\t{1}\t{2}\t".format(info["symbol"],gene,switch.nTx))
				self.GPS.write("{0}\t{1}\t{2}\t".format(switch.tTx,whatsHappening,thisRes.num))
				self.GPS.write("{0}\t{1}\n".format(thisRes.res,",".join(thisRes._kinases)))
				anyPTM = True

		return anyPTM
	'''
	'''
	def findBrokenSurfaces(self,switch,gene,info):

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
				self.I3D.write("{0}\t{1}\t{2}\t{3}\t".format(gene,info["symbol"],protein.tx,protein.uniprot))
				self.I3D.write("{0}\t{1}\t{2}\t{3}\n".format(percent,pval,isoInfo,isoSpec))
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
	'''
	def newPrositeAnalysis(self,thisSwitch,gene,info):
		anyPTM = False
		ptms = set()

		for tx,iso in zip([thisSwitch.nTranscript,thisSwitch.tTranscript],[thisSwitch.nIsoform,thisSwitch.tIsoform]):
			tx.readProsite()

			for ptm in tx._ptms:
				for start,end in tx._ptms[ptm]:
					for resNum in range(start,end):
						thisRes = iso._structure[resNum-1]
						thisRes._ptms.append(ptm)

				ptms.add(ptm)

		for ptm in ptms:
			nPtmRegions = thisSwitch.nIsoform.getSegments(ptm,minLength=1,gap=0)
			tPtmRegions = thisSwitch.tIsoform.getSegments(ptm,minLength=1,gap=0)
			nIsospRegions = thisSwitch.nIsoform.getSegments("isoform-specific")
			tIsospRegions = thisSwitch.tIsoform.getSegments("isoform-specific")

			nReps = 0 if not nPtmRegions else len(nPtmRegions)
			tReps = 0 if not tPtmRegions else len(tPtmRegions)
			maxReps = nReps if nReps > tReps else tReps

			for i in range(maxReps):

				nMacroScore = 0
				nMicroScore = 0
				nJaccard = 0
				tMacroScore = 0
				tMicroScore = 0
				tJaccard = 0
				
				if i+1 <= nReps:
					nThisIsosp = []
					[ nThisIsosp.extend(x) for x in nIsospRegions if set(x) & set(nPtmRegions[i]) ]
					intersection = float(len(set(nPtmRegions[i]) & set(nThisIsosp)))
					nPtmLength = float(len(set(nPtmRegions[i])))
					nIsoSpLength = float(len(set(nThisIsosp)))
					nMacroScore = intersection/nPtmLength
					nMicroScore =  "NA" if nIsoSpLength==0 else intersection/nIsoSpLength
					nJaccard = intersection/len(set(nPtmRegions[i]) | set(nThisIsosp))
					
				if i+1 <= tReps:
					tThisIsosp = []
					[ tThisIsosp.extend(x) for x in tIsospRegions if set(x) & set(tPtmRegions[i]) ]
					intersection = float(len(set(tPtmRegions[i]) & set(tThisIsosp)))
					tPtmLength = float(len(set(tPtmRegions[i])))
					tIsoSpLength = float(len(set(tThisIsosp)))
					tMacroScore = intersection/tPtmLength
					tMicroScore = "NA" if tIsoSpLength==0 else intersection/tIsoSpLength
					tJaccard = intersection/len(set(tPtmRegions[i]) | set(tThisIsosp))

				if i+1 <= nReps and i+1 <= tReps:
					whatsHappening = "Nothing"
				elif i+1 <= nReps:
					whatsHappening = "Lost_in_tumor"
					anyPTM = True
				elif i+1 <= tReps:
					whatsHappening = "Gained_in_tumor"
					anyPTM = True

				self.PROSITE.write("{0}\t{1}\t{2}\t".format(gene,info["symbol"],thisSwitch.nTx))
				self.PROSITE.write("{0}\t{1}\t{2}\t".format(thisSwitch.tTx,whatsHappening,ptm))
				self.PROSITE.write("{0}/{1}\t{2}/{3}\t".format(i+1,nReps,i+1,tReps))
				self.PROSITE.write("{0}\t{1}\t{2}\t".format(nMacroScore,nMicroScore,nJaccard))
				self.PROSITE.write("{0}\t{1}\t{2}\n".format(tMacroScore,tMicroScore,tJaccard))

		return anyPTM

	def newJoinFiles(self,tag=""):

		roots = ['interpro_analysis','iupred_analysis','anchor_analysis','prosite_analysis','structural_summary']
		
		for root in roots:
			outFile = "{0}structural_analysis/{1}{2}.tsv".format(options.Options().qout,root,tag)
			files = glob.glob("{0}structural_analysis/{1}{2}_[0-9]*.tsv".format(options.Options().qout,root,tag))

			# close writing files in case they are opened
			if (hasattr(self,'IP') and hasattr(self,'IU') and hasattr(self,'ANCHOR') and hasattr(self,'PROSITE') and hasattr(self,'REL')):
				self.IP.close()
				self.IU.close()
				self.ANCHOR.close()
				self.PROSITE.close()
				self.REL.close()
		
			with open(outFile,"w") as OUT:
				if root=='interpro_analysis': self.writeIPHeader(OUT)
				elif root=='iupred_analysis': self.writeIUHeader(OUT)
				elif root=='anchor_analysis': self.writeANCHORHeader(OUT)
				elif root=='prosite_analysis': self.newWritePrositeHeader(OUT)
				elif root=='structural_summary': self.writeSummaryHeader(OUT)
						
				for aFile in files:
					with open(aFile) as IN:
						IN.readline()
						OUT.write(IN.read())

	def newRecalculateRelevance(self,switches,tag):
		files = glob.glob("{0}structural_analysis/structural_summary{1}_[0-9]*.tsv".format(options.Options().qout,tag))
		outFile = "{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,tag)

		with open(outFile,"w") as OUT:
			self.writeSummaryHeader(OUT)
			for aFile in files:
				for line in utils.readTable(aFile):
					switchId = "{0}_{1}_{2}_{3}".format(line[0],line[1],line[2],line[3])

					gene = line[0]
					symbol = line[1]
					ntx = line[2]
					ttx = line[3]
					rLoop = line[4]
					rDomain = line[5]
					rProsite = line[8]

					if switchId in switches['iupred']:
						rDisorder = "True"
					else:
						rDisorder = "False"
					
					if switchId in switches['anchor']:
						rAnchor = "True"
					else:
						rAnchor = "False"

					OUT.write("{0}\t{1}\t".format(gene,symbol))
					OUT.write("{0}\t{1}\t".format(ntx,ttx))
					OUT.write("{0}\t{1}\t".format(rLoop,rDomain))
					OUT.write("{0}\t{1}\t".format(rDisorder,rAnchor))
					OUT.write("{0}\n".format(rProsite))
		
if __name__ == '__main__':
	pass
