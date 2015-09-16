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
			tag = "_random"
			if options.Options().parallelRange:
				tag += "_{}".format(options.Options().parallelRange)
		else:
			tag = ""
			if not options.Options().parallelRange:
				self.joinFiles("_random")
				
			else:
				tag = "_{}".format(options.Options().parallelRange)

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
		self.writePrositeHeader(self.PROSITE)
						
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
			thisSwitch._ptm_change 		  = self.prositeAnalysis(thisSwitch,gene,info)

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
				protein.readIupred(mode)

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
			protein.readAnchor()
			
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

				nMacroScore = "NA"
				nMicroScore = "NA"
				nJaccard = "NA"
				tMacroScore = "NA"
				tMicroScore = "NA"
				tJaccard = "NA"

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
		PROSITE.write("tTumorTranscript\tWhat\tFeature\t")
		PROSITE.write("normalReps\ttumorReps\tnMacroScore\t")
		PROSITE.write("nMicroScore\tnJaccard\ttMacroScore\t")
		PROSITE.write("tMicroScore\ttJaccard\n")
		
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
				if root=='interpro_analysis': self.writeIPHeader(OUT)
				elif root=='iupred_analysis': self.writeIUHeader(OUT)
				elif root=='anchor_analysis': self.writeANCHORHeader(OUT)
				elif root=='prosite_analysis': self.writePrositeHeader(OUT)
				elif root=='structural_summary': self.writeSummaryHeader(OUT)
						
				for aFile in files:
					with open(aFile) as IN:
						IN.readline()
						OUT.write(IN.read())
		
if __name__ == '__main__':
	pass
