from interface import interpro_analysis
from libs import options
from libs import utils
from methods import method

from collections import Counter
import glob
import operator

class StructuralAnalysis(method.Method):
	def __init__(self,gn_network,tx_network,isRand=False):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self.isRandom = isRand
		if self.isRandom:
			tag = "_random"
			
		else:
			tag = ""
			if not options.Options().parallelRange:
				self.joinFiles("_random")

		if options.Options().parallelRange:
			tag += "_{}".format(options.Options().parallelRange)

		self.anchor_threshold = 0.5
		self.iupred_threshold = 0.5

		self._interpro_file = "{0}structural_analysis/interpro_analysis{1}.tsv".format(options.Options().qout,tag)
		self.IP = open(self._interpro_file,"w")
		self.writeKnownFeatureHeader(self.IP)
				
		self._iupred_file = "{0}structural_analysis/iupred_analysis{1}.tsv".format(options.Options().qout,tag)
		self.IU = open(self._iupred_file,"w")
		self.writeDisorderedRegionHeader(self.IU)
								
		self._anchor_file = "{0}structural_analysis/anchor_analysis{1}.tsv".format(options.Options().qout,tag)
		self.ANCHOR = open(self._anchor_file,"w")
		self.writeDisorderedRegionHeader(self.ANCHOR)
				
		self._prosite_file = "{0}structural_analysis/prosite_analysis{1}.tsv".format(options.Options().qout,tag)
		self.PROSITE = open(self._prosite_file,"w")
		self.writeKnownFeatureHeader(self.PROSITE)
						
		self._relevance_info = "{0}structural_analysis/structural_summary{1}.tsv".format(options.Options().qout,tag)
		self.REL = open(self._relevance_info,"w")
		self.writeSummaryHeader(self.REL)
				
	def run(self):		
		self.logger.info("Structural analysis.")

		# get loops families
		isoInfo = {}
		for line in open("{}data/{}/sequences.uniprot.loops.fa".format(
			options.Options().wd,options.Options().annotation)):
			if ">" in line:
				elements = line.strip().split("#")
				isoInfo[elements[0][1:]] = {}
				isoInfo[elements[0][1:]]["UniProt"] = elements[2]
				isoInfo[elements[0][1:]]["iLoopsFamily"] = elements[3]

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,partialCreation=True,removeNoise=True):
			
			if not thisSwitch.nIsoform or not thisSwitch.tIsoform: 
				continue

			thisSwitch._iloops_change 	  							= self.archDBAnalysis(thisSwitch,gene,info,isoInfo)
			(thisSwitch._functional_change,thisSwitch._ptm_change) 	= self.knownFeaturesAnalysis(thisSwitch,gene,info)
			thisSwitch._disorder_change   							= self.disorderAnalysis(thisSwitch,gene,info)
			thisSwitch._anchor_change    							= self.anchorAnalysis(thisSwitch,gene,info)

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

	def disorderAnalysis(self,thisSwitch,gene,info):

		self.logger.debug("IUPRED: Searching disorder for gene {0}.".format(gene))
		
		anyIUpredSeq = False
		normalProtein = thisSwitch.nIsoform
		tumorProtein = thisSwitch.tIsoform

		if not normalProtein or not tumorProtein: return False

		for protein,whatShouldBeHappening in zip([normalProtein,tumorProtein],["Lost_in_tumor","Gained_in_tumor"]):
			for mode in ["short","long"]:
				protein.readIupred(mode)

				disordered = protein.getSegments("disordered",minLength=5,gap=2)
				isoform = protein.getSegments("isoform-specific")

				for disorderedRegion in disordered:

					whatsHappening = whatShouldBeHappening
					disorderedRegionSet = set(disorderedRegion)

					overlappingIsoSpecific = []
					[ overlappingIsoSpecific.extend(x) for x in isoform if set(x) & disorderedRegionSet]
					
					if not overlappingIsoSpecific:
						jaccard = "NA"
						macroScore = "NA"
						microScore = "NA"
						whatsHappening = "Nothing"
						significant = 0

					else:
						overlappingIsoSpecificSet = set(overlappingIsoSpecific)

						intersection = float(len(overlappingIsoSpecificSet & disorderedRegionSet))
						union = float(len(overlappingIsoSpecificSet | disorderedRegionSet))

						jaccard = intersection/union
						microScore = intersection/len(overlappingIsoSpecificSet)
						macroScore = intersection/len(disorderedRegionSet)
						significant = int(max(microScore,macroScore) > self.iupred_threshold)

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

					self.IU.write("{0}\t{1}\t".format(gene,info["symbol"]))
					self.IU.write("{0}\t{1}\t".format(normalProtein.tx,tumorProtein.tx))
					self.IU.write("{0}\t{1}\t".format(whatsHappening,motifSequence))
					self.IU.write("{0}\t{1}\t".format(start,end))
					self.IU.write("{0}\t{1}\t".format(jaccard,microScore))
					self.IU.write("{0}\t{1}\n".format(macroScore,significant))

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
			protein.readAnchor()
			
			anchor = protein.getSegments("anchor",minLength=5,gap=2)
			isoform = protein.getSegments("isoform-specific")

			for anchorRegion in anchor:
				whatsHappening = whatShouldBeHappening
				anchorRegionSet = set(anchorRegion)

				overlappingIsoSpecific = []
				[ overlappingIsoSpecific.extend(x) for x in isoform if set(x) & anchorRegionSet]
				
				if not overlappingIsoSpecific:
					jaccard = "NA"
					macroScore = "NA"
					microScore = "NA"
					whatsHappening = "Nothing"
					significant = 0
				else:
					overlappingIsoSpecificSet = set(overlappingIsoSpecific)

					intersection = float(len(overlappingIsoSpecificSet & anchorRegionSet))
					union = float(len(overlappingIsoSpecificSet | anchorRegionSet))

					jaccard = intersection/union
					microScore = intersection/len(overlappingIsoSpecificSet)
					macroScore = intersection/len(anchorRegionSet)
					significant = int(max(microScore,macroScore) > self.iupred_threshold)

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

				self.ANCHOR.write("{0}\t{1}\t".format(gene,info["symbol"]))
				self.ANCHOR.write("{0}\t{1}\t".format(normalProtein.tx,tumorProtein.tx))
				self.ANCHOR.write("{0}\t{1}\t".format(whatsHappening,motifSequence))
				self.ANCHOR.write("{0}\t{1}\t".format(start,end))
				self.ANCHOR.write("{0}\t{1}\t".format(jaccard,microScore))
				self.ANCHOR.write("{0}\t{1}\n".format(macroScore,significant))

				if significant:
					anyAnchorSeq = True

		return anyAnchorSeq

	def knownFeaturesAnalysis(self,thisSwitch,gene,info):
		anyFeature = {"Prosite":False,"InterPro":False}
		features = {"Prosite":set(),"InterPro":set()}
		prosites = set()
		interpros = set()

		for isoform in [thisSwitch.nIsoform,thisSwitch.tIsoform]:
			isoform.readProsite()
			isoform.readInterpro()

		[ features["Prosite"].add(x) for i in [thisSwitch.nIsoform,thisSwitch.tIsoform] for x in i._prosite ]
		[ features["InterPro"].add("{0}|{1}".format(x['accession'],x['description'])) for i in [thisSwitch.nIsoform,thisSwitch.tIsoform] for x in i._pfam ]

		for OUT,featType in zip([self.IP,self.PROSITE],["InterPro","Prosite"]):
			for feature in features[featType]:
				featInfo = {}
				for isoform in [thisSwitch.nIsoform,thisSwitch.tIsoform]:
					featInfo[isoform.tx] = []
					featRegions = isoform.getSegments(feature,minLength=1,gap=0)
					isospRegions = isoform.getSegments("isoform-specific")

					for region in featRegions:

						macroScore = float("-inf")
						microScore = float("-inf")
						jaccard = float("-inf")

						thisIsosp = []
						[ thisIsosp.extend(x) for x in isospRegions if set(x) & set(region) ]
						intersection = float(len(set(region) & set(thisIsosp)))
						featLength = float(len(set(region)))
						isoSpLength = float(len(set(thisIsosp)))
						macroScore = intersection/featLength
						microScore =  float("-inf") if isoSpLength==0 else intersection/isoSpLength
						jaccard = intersection/len(set(region) | set(thisIsosp))

						featInfo[isoform.tx].append({"macro": macroScore, 
										   "micro": microScore, "jaccard": jaccard})

					featInfo[isoform.tx] = sorted(featInfo[isoform.tx], key=operator.itemgetter("macro"))

				featInfoZipped = map(None, featInfo[thisSwitch.nIsoform.tx], featInfo[thisSwitch.tIsoform.tx])

				for i in range(len(featInfoZipped)):

					nDict = featInfoZipped[i][0]
					tDict = featInfoZipped[i][1]

					nMacroScore = "NA"
					nMicroScore = "NA"
					nJaccard = "NA"
					tMacroScore = "NA"
					tMicroScore = "NA"
					tJaccard = "NA"

					whatsHappening = "Nothing"

					if nDict:
						nMacroScore = "NA" if nDict["macro"] < 0 else nDict["macro"]
						nMicroScore = "NA" if nDict["micro"] < 0 else nDict["micro"]
						nJaccard = "NA" if nDict["jaccard"] < 0 else nDict["jaccard"]

					if tDict:
						tMacroScore = "NA" if tDict["macro"] < 0 else tDict["macro"]
						tMicroScore = "NA" if tDict["micro"] < 0 else tDict["micro"]
						tJaccard = "NA" if tDict["jaccard"] < 0 else tDict["jaccard"]

					if nDict and tDict is None:
						if nMacroScore != "NA" and nMacroScore > 0:
							whatsHappening = "Lost_in_tumor"
							anyFeature[featType] = True
					elif tDict and nDict is None:
						if tMacroScore != "NA" and tMacroScore > 0:
							whatsHappening = "Gained_in_tumor"
							anyFeature[featType] = True

					OUT.write("{0}\t{1}\t{2}\t".format(gene,info["symbol"],thisSwitch.nTx))
					OUT.write("{0}\t{1}\t{2}\t".format(thisSwitch.tTx,whatsHappening,feature))
					OUT.write("{0}/{1}\t".format(i+1,len(featInfo[thisSwitch.nIsoform.tx])))
					OUT.write("{0}/{1}\t".format(i+1,len(featInfo[thisSwitch.tIsoform.tx])))
					OUT.write("{0}\t{1}\t{2}\t".format(nMacroScore,nMicroScore,nJaccard))
					OUT.write("{0}\t{1}\t{2}\n".format(tMacroScore,tMicroScore,tJaccard))

		return anyFeature["InterPro"],anyFeature["Prosite"]

	def writeKnownFeatureHeader(self,OUT):
		OUT.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
		OUT.write("What\tFeature\tnormalReps\ttumorReps\tnMacroScore\t")
		OUT.write("nMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard\n")

	def writeDisorderedRegionHeader(self,OUT):
		OUT.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
		OUT.write("What\tSequence\tStartPos\tEndPos\tJaccard\t")
		OUT.write("microScore\tmacroScore\tSignificant\n")
		
	def writeSummaryHeader(self,REL):
		REL.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
		REL.write("iLoops\tDomain\tDisorder\tAnchor\tPTM\n")

	def joinFiles(self,tag=""):

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
				if root=='interpro_analysis': self.writeKnownFeatureHeader(OUT)
				elif root=='iupred_analysis': self.writeDisorderedRegionHeader(OUT)
				elif root=='anchor_analysis': self.writeDisorderedRegionHeader(OUT)
				elif root=='prosite_analysis': self.writeKnownFeatureHeader(OUT)
				elif root=='structural_summary': self.writeSummaryHeader(OUT)
						
				for aFile in files:
					with open(aFile) as IN:
						IN.readline()
						OUT.write(IN.read())
		
if __name__ == '__main__':
	pass
