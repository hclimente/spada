from interface import interpro_analysis
from libs import options
from libs import utils
from methods import method

from collections import Counter
import glob
import operator
import os

class StructuralAnalysis(method.Method):
	def __init__(self,gn_network,tx_network,isRand=False):

		if not os.path.exists("{}structural_analysis".format(options.Options().qout)):
			utils.cmd("mkdir","{}structural_analysis".format(options.Options().qout))

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

		self._interpro_file = "{}structural_analysis/interpro_analysis{}.tsv".format(options.Options().qout,tag)
		self.IP = open(self._interpro_file,"w")
		self.writeKnownFeatureHeader(self.IP)
				
		self._iupred_file = "{}structural_analysis/iupred_analysis{}.tsv".format(options.Options().qout,tag)
		self.IU = open(self._iupred_file,"w")
		self.writeDisorderedRegionHeader(self.IU)
								
		self._anchor_file = "{}structural_analysis/anchor_analysis{}.tsv".format(options.Options().qout,tag)
		self.ANCHOR = open(self._anchor_file,"w")
		self.writeDisorderedRegionHeader(self.ANCHOR)
				
		self._prosite_file = "{}structural_analysis/prosite_analysis{}.tsv".format(options.Options().qout,tag)
		self.PROSITE = open(self._prosite_file,"w")
		self.writeKnownFeatureHeader(self.PROSITE)
						
		self._relevance_info = "{}structural_analysis/structural_summary{}.tsv".format(options.Options().qout,tag)
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

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network,partialCreation=True,removeNoise=False):
			thisSwitch._iloops_change 	  							= self.archDBAnalysis(thisSwitch,gene,info,isoInfo)
			(thisSwitch._functional_change,thisSwitch._ptm_change) 	= self.knownFeaturesAnalysis(thisSwitch,gene,info)
			thisSwitch._disorder_change   							= self.disorderAnalysis(thisSwitch,gene,info)
			thisSwitch._anchor_change    							= self.anchorAnalysis(thisSwitch,gene,info)

			self.REL.write("{}\t{}\t".format(gene,info["symbol"]))
			self.REL.write("{}\t{}\t".format(thisSwitch.nTx,thisSwitch.tTx))
			self.REL.write("{}\t{}\t".format(thisSwitch._iloops_change,thisSwitch._functional_change))
			self.REL.write("{}\t{}\t".format(thisSwitch._disorder_change,thisSwitch._anchor_change))
			self.REL.write("{}\n".format(thisSwitch._ptm_change))

		self.IP.close()
		self.IU.close()
		self.ANCHOR.close()
		self.PROSITE.close()
		self.REL.close()

	def archDBAnalysis(self,thisSwitch,gene,info,isoInfo):

		self.logger.debug("iLoops: looking for loop changes for gene {}.".format(gene) )
		
		if thisSwitch.nTx in isoInfo and thisSwitch.tTx in isoInfo:
			nLoops = isoInfo[thisSwitch.nTx]["iLoopsFamily"]
			tLoops = isoInfo[thisSwitch.tTx]["iLoopsFamily"]

			if nLoops != tLoops:
				self.logger.debug("iLoops: information found for gene {}.".format(gene) )
				return True
			else:
				return False

		elif thisSwitch.nTx in isoInfo and thisSwitch.tIsoform is None:
			if isoInfo[thisSwitch.nTx]["iLoopsFamily"]:
				return True

		elif thisSwitch.tTx in isoInfo and thisSwitch.nIsoform is None:
			if isoInfo[thisSwitch.tTx]["iLoopsFamily"]:
				return True

		return False

	def disorderAnalysis(self,thisSwitch,gene,info):

		self.logger.debug("IUPRED: Searching disorder for gene {}.".format(gene))
		
		anyIUpredSeq = False

		for protein,whatShouldBeHappening in zip([thisSwitch.nIsoform,thisSwitch.tIsoform],["Lost_in_tumor","Gained_in_tumor"]):
			if not protein:
				continue

			for mode in ["short","long"]:
				protein.readIupred(mode)

				disordered = protein.getSegments("disordered",minLength=5,gap=2)
				isoform = protein.getSegments("isoform-specific")

				for disorderedRegion in disordered:

					whatsHappening = whatShouldBeHappening
					disorderedRegionSet = set(disorderedRegion)

					overlappingIsoSpecific = []
					if None not in [thisSwitch.nIsoform,thisSwitch.tIsoform]:
						[ overlappingIsoSpecific.extend(x) for x in isoform if set(x) & disorderedRegionSet]
					else:
						overlappingIsoSpecific = disorderedRegionSet
					
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

					self.IU.write("{}\t{}\t".format(gene,info["symbol"]))
					self.IU.write("{}\t{}\t".format(thisSwitch.nTx,thisSwitch.tTx))
					self.IU.write("{}\t{}\t".format(whatsHappening,motifSequence))
					self.IU.write("{}\t{}\t".format(start,end))
					self.IU.write("{}\t{}\t".format(jaccard,microScore))
					self.IU.write("{}\t{}\n".format(macroScore,significant))

					if significant:
						anyIUpredSeq = True

		return anyIUpredSeq

	def anchorAnalysis(self,thisSwitch,gene,info):

		self.logger.debug("ANCHOR: Searching anchoring regions for gene {}.".format(gene))
		
		anyAnchorSeq = False

		for protein,whatShouldBeHappening in zip([thisSwitch.nIsoform,thisSwitch.tIsoform],["Lost_in_tumor","Gained_in_tumor"]):
			if not protein:
				continue
			#Parse anchor output
			protein.readAnchor()
			
			anchor = protein.getSegments("anchor",minLength=5,gap=2)
			isoform = protein.getSegments("isoform-specific")

			for anchorRegion in anchor:
				whatsHappening = whatShouldBeHappening
				anchorRegionSet = set(anchorRegion)

				overlappingIsoSpecific = []
				if None not in [thisSwitch.nIsoform,thisSwitch.tIsoform]:
					[ overlappingIsoSpecific.extend(x) for x in isoform if set(x) & anchorRegionSet]
				else:
					overlappingIsoSpecific = anchorRegionSet
				
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

				self.ANCHOR.write("{}\t{}\t".format(gene,info["symbol"]))
				self.ANCHOR.write("{}\t{}\t".format(thisSwitch.nTx,thisSwitch.tTx))
				self.ANCHOR.write("{}\t{}\t".format(whatsHappening,motifSequence))
				self.ANCHOR.write("{}\t{}\t".format(start,end))
				self.ANCHOR.write("{}\t{}\t".format(jaccard,microScore))
				self.ANCHOR.write("{}\t{}\n".format(macroScore,significant))

				if significant:
					anyAnchorSeq = True

		return anyAnchorSeq

	def knownFeaturesAnalysis(self,thisSwitch,gene,info):
		anyFeature = {"prosite":False,"pfam":False}
		features = {"prosite":set(),"pfam":set()}
		prosites = set()
		interpros = set()

		for isoform in [thisSwitch.nIsoform,thisSwitch.tIsoform]:
			if not isoform:
				continue
			isoform.readProsite()
			isoform.readInterpro()
			[ features["prosite"].add(x) for x in isoform._prosite ]
			[ features["pfam"].add("{}|{}".format(x['accession'],x['description'])) for x in isoform._pfam ]

		for OUT,featType in zip([self.IP,self.PROSITE],["pfam","prosite"]):
			for feature in features[featType]:
				featInfo = { thisSwitch.nTx: [], thisSwitch.tTx: []}
				for isoform in [thisSwitch.nIsoform,thisSwitch.tIsoform]:
					if not isoform:
						continue

					featRegions = isoform.getSegments(feature,minLength=1,gap=0)
					isospRegions = isoform.getSegments("isoform-specific")

					for region in featRegions:

						macroScore = float("-inf")
						microScore = float("-inf")
						jaccard = float("-inf")

						thisIsosp = []
						if None not in [thisSwitch.nIsoform,thisSwitch.tIsoform]:
							[ thisIsosp.extend(x) for x in isospRegions if set(x) & set(region) ]
						else:
							thisIsosp = set(region)
						intersection = float(len(set(region) & set(thisIsosp)))
						featLength = float(len(set(region)))
						isoSpLength = float(len(set(thisIsosp)))
						macroScore = intersection/featLength
						microScore =  float("-inf") if isoSpLength==0 else intersection/isoSpLength
						jaccard = intersection/len(set(region) | set(thisIsosp))

						featInfo[isoform.tx].append({"macro": macroScore,"micro": microScore, "jaccard": jaccard})

					featInfo[isoform.tx] = sorted(featInfo[isoform.tx], key=operator.itemgetter("macro"))

				i = 1
				for nDict,tDict in zip(featInfo[thisSwitch.nTx], featInfo[thisSwitch.tTx]):

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
						if (nMacroScore != "NA" and nMacroScore > 0):
							whatsHappening = "Lost_in_tumor"
							anyFeature[featType] = True
					elif tDict and nDict is None:
						if (tMacroScore != "NA" and tMacroScore > 0):
							whatsHappening = "Gained_in_tumor"
							anyFeature[featType] = True

					OUT.write("{}\t{}\t{}\t".format(gene,info["symbol"],thisSwitch.nTx))
					OUT.write("{}\t{}\t{}\t".format(thisSwitch.tTx,whatsHappening,feature))
					OUT.write("{}/{}\t".format(i,len(featInfo[thisSwitch.nTx])))
					OUT.write("{}/{}\t".format(i,len(featInfo[thisSwitch.tTx])))
					OUT.write("{}\t{}\t{}\t".format(nMacroScore,nMicroScore,nJaccard))
					OUT.write("{}\t{}\t{}\n".format(tMacroScore,tMicroScore,tJaccard))
					i += 1

		return anyFeature["pfam"],anyFeature["prosite"]

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
			outFile = "{}structural_analysis/{}{}.tsv".format(options.Options().qout,root,tag)
			files = glob.glob("{}structural_analysis/{}{}_[0-9]*.tsv".format(options.Options().qout,root,tag))

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

		if tag=="":
			saveName = "geneNetwork.pkl"
			g = self._gene_network
		elif tag=="_random":
			saveName = "randomGeneNetwork_fixNormal.pkl"
			import pickle
			g = pickle.load(open("{}{}".format(options.Options().qout,saveName),"rb"))
			g.createLogger()

		for elements in utils.readTable("{}structural_analysis/structural_summary{}.tsv".format(options.Options().qout,tag)):
			gene = elements[0]
			nTx = elements[2]
			tTx = elements[3]

			for switchDict in g._net.node[gene]["isoformSwitches"]:
				if switchDict["nIso"]==nTx and switchDict["tIso"]==tTx:
					thisSwitch = g.createSwitch(switchDict,self._transcript_network,True)
					switchDict["functional"] = thisSwitch.is_functional

		g.saveNetwork(saveName)
		
if __name__ == '__main__':
	pass
