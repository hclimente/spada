#!/soft/devel/python-2.7/bin/python

from interface import interpro_analysis
from libs import options
from libs import utils
from methods import method

from collections import Counter
from Bio import pairwise2
import fisher
import glob
import numpy as np
import os

import pdb

class StructuralAnalysis(method.Method):
	def __init__(self,gn_network,tx_network,isRand=False):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self.isRandom = isRand
		if self.isRandom:
			self.anchor_threshold = 0
			self.iupred_threshold = 0
			tag = "_random"
		else:
			tag = ""
			self.joinRandomFiles()
			#self.anchor_threshold = self.getThresholdFromRandom("anchor",95)
			#self.iupred_threshold = self.getThresholdFromRandom("iupred",95)

		self._interpro_file = "{0}structural_analysis/InterPro_report{1}{2}.tsv".format(options.Options().qout,tag,options.Options().filetag)
		self.IP = open(self._interpro_file,"w")
		self.IP.write("Symbol\tGene\tNormalTranscript\ttTumorTranscript\tAnalysis\tFeature_accesion\t")
		self.IP.write("Feature\t(Additional) repetitions\n")

		self._iupred_file = "{0}structural_analysis/iupred_analysis{1}{2}.tsv".format(options.Options().qout,tag,options.Options().filetag)
		self.IU = open(self._iupred_file,"w")
		self.IU.write("Symbol\tGene\tNormalTranscript\ttTumorTranscript\t")
		self.IU.write("What\tSequence\tJaccard\tThreshold\tSignificant\n")
				
		self._anchor_file = "{0}structural_analysis/anchor_analysis{1}{2}.tsv".format(options.Options().qout,tag,options.Options().filetag)
		self.ANCHOR = open(self._anchor_file,"w")
		self.ANCHOR.write("Symbol\tGene\tNormalTranscript\ttTumorTranscript\t")
		self.ANCHOR.write("What\tSequence\tJaccard\tThreshold\tSignificant\n")

		self._prosite_file = "{0}structural_analysis/prosite_analysis{1}{2}.tsv".format(options.Options().qout,tag,options.Options().filetag)
		self.PROSITE = open(self._prosite_file,"w")
		self.PROSITE.write("Symbol\tGene\tNormalTranscript\t")
		self.PROSITE.write("tTumorTranscript\tWhat\tMotif\n")

		self._relevance_info = "{0}structural_analysis/structural_summary{1}{2}.tsv".format(options.Options().qout,tag,options.Options().filetag)
		self.REL = open(self._relevance_info,"w")
		self.REL.write("Gene\tNormalTranscript\tTumorTranscript\tiLoops\t")
		self.REL.write("Domain\tDisorder\tAnchor\tPTM\n")

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


		for gene,info,switchDict,switch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,partialCreation=True,removeNoise=True):
			switch._iloops_change 	  = self.findDiffLoops(switch,gene,info,isoInfo)
			switch._functional_change = self.interProAnalysis(switch,gene,info)
			switch._disorder_change   = self.disorderAnalysis(switch,gene,info)
			switch._anchor_change     = self.anchorAnalysis(switch,gene,info)
			switch._ptm_change 		  = self.prositeAnalysis(switch,gene,info)

			self.REL.write("{0}\t{1}\t{2}\t".format(gene,switch.nTx,switch.tTx))
			self.REL.write("{0}\t{1}\t".format(switch._iloops_change,switch._functional_change))
			self.REL.write("{0}\t{1}\t".format(switch._disorder_change,switch._anchor_change))
			self.REL.write("{0}\n".format(switch._ptm_change))

		self.IP.close()
		self.IU.close()
		self.ANCHOR.close()
		self.REL.close()

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

		return False

	def interProAnalysis(self,switch,gene,info):

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

				self.IP.write("{0}\t{1}\t{2}\t".format(info["symbol"],gene,switch.nTx))
				self.IP.write("{0}\t{1}\t".format(switch.tTx,whatsHappening))
				self.IP.write("{0}\t{1}\t".format(featInfo["analysis"],featInfo["accession"]))
				self.IP.write("{0}\t{1}\n".format(featInfo["description"],reps))
				# self.IP.write("{0}\n".format(featInfo["percentAffected"]))

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

	def disorderAnalysis(self,switch,gene,info):

		self.logger.debug("IUPRED: Searching disorder for gene {0}.".format(gene))
		
		anyIUpredSeq = False
		normalProtein = switch.nIsoform
		tumorProtein = switch.tIsoform

		if not normalProtein or not tumorProtein: return False

		for protein,whatsHappening in zip([normalProtein,tumorProtein],["Lost in tumor","Gained in tumor"]):
			fFile = "{0}{1}.fa".format(options.Options().qout,protein.tx)
			with open(fFile,"w") as FASTA:
				FASTA.write(">{0}\n{1}\n".format(protein.tx,protein.seq))

			for mode in ["short","long"]:
				proc = utils.cmdOut(options.Options().wd+"Pipeline/libs/bin/iupred/iupred",fFile,mode)

				#Parse iupred output
				for line in [ x.strip().split(" ") for x in proc.stdout if "#" not in x ]:
					resNum  = int(line[0])
					residue	= line[1]
					score 	= float(line[-1])

					thisRes = protein._structure[resNum-1]
					if residue != thisRes.res:
						self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(protein.tx))
						os.remove(fFile)
						continue
					thisRes.set_iuPredScore(score)

				disordered = protein.getSegments("disordered",minLength=5,gap=2)
				isoform = protein.getSegments("isoform-specific")

				for disorderedRegion in disordered:
					disorderedRegionSet = set(disorderedRegion)

					overlappingIsoSpecific = []

					for isoSpecific in isoform:
						isoSpecificSet = set(isoSpecific)

						if isoSpecificSet & disorderedRegionSet:
							overlappingIsoSpecific.extend(isoSpecific)
					
					if not overlappingIsoSpecific:
						continue

					overlappingIsoSpecificSet = set(overlappingIsoSpecific)

					intersection = float(len(overlappingIsoSpecificSet & disorderedRegionSet))
					union = float(len(overlappingIsoSpecificSet | disorderedRegionSet))

					jaccard = intersection/union

					motifSequence = ""
					for thisRes in disorderedRegion:
						res = thisRes.res.upper() if thisRes.isoformSpecific else thisRes.res.lower()
						motifSequence += res

					self.IU.write("{0}\t{1}\t".format(gene,info["symbol"]))
					self.IU.write("{0}\t{1}\t".format(normalProtein.tx,tumorProtein.tx))
					self.IU.write("{0}\t{1}\t".format(whatsHappening,motifSequence))
					self.IU.write("{0}\t{1}\t".format(jaccard,self.iupred_threshold))
					self.IU.write("{0}\n".format(int(jaccard > self.iupred_threshold)))

					if jaccard > self.iupred_threshold:
						anyIUpredSeq = True
					
			os.remove(fFile)

		return anyIUpredSeq

	def anchorAnalysis(self,switch,gene,info):

		self.logger.debug("ANCHOR: Searching anchoring regions for gene {0}.".format(gene))
		
		anyAnchorSeq = False
		normalProtein = switch.nIsoform
		tumorProtein = switch.tIsoform

		if not normalProtein or not tumorProtein: return False

		for protein,whatsHappening in zip([normalProtein,tumorProtein],["Lost in tumor","Gained in tumor"]):
			fFile = "{0}{1}.fa".format(options.Options().qout,protein.tx)
			with open(fFile,"w") as FASTA:
				FASTA.write(">{0}\n{1}\n".format(protein.tx,protein.seq))

			proc = utils.cmdOut(options.Options().wd+"Pipeline/libs/ANCHOR/anchor",fFile)

			#Parse anchor output
			for line in [ x.strip().split("\t") for x in proc.stdout if "#" not in x ]:
				resNum  = int(line[0])
				residue	= line[1]
				score 	= float(line[2].strip())

				thisRes = protein._structure[resNum-1]

				if residue != thisRes.res:
					self.logger.error("Not matching residue in ANCHOR analysis, transcript {0}.".format(protein.tx))
					os.remove(fFile)
					continue

				thisRes.set_anchorScore(score)

			anchor = protein.getSegments("anchor",minLength=5,gap=2)
			isoform = protein.getSegments("isoform-specific")

			for anchorRegion in anchor:
				anchorRegionSet = set(anchorRegion)

				overlappingIsoSpecific = []

				for isoSpecific in isoform:
					isoSpecificSet = set(isoSpecific)

					if isoSpecificSet & anchorRegionSet:
						overlappingIsoSpecific.extend(isoSpecific)
				
				if not overlappingIsoSpecific:
					continue

				overlappingIsoSpecificSet = set(overlappingIsoSpecific)

				intersection = float(len(overlappingIsoSpecificSet & anchorRegionSet))
				union = float(len(overlappingIsoSpecificSet | anchorRegionSet))

				jaccard = intersection/union

				motifSequence = ""
				for thisRes in anchorRegion:
					res = thisRes.res.upper() if thisRes.isoformSpecific else thisRes.res.lower()
					motifSequence += res

				self.ANCHOR.write("{0}\t{1}\t".format(gene,info["symbol"]))
				self.ANCHOR.write("{0}\t{1}\t".format(normalProtein.tx,tumorProtein.tx))
				self.ANCHOR.write("{0}\t{1}\t".format(whatsHappening,motifSequence))
				self.ANCHOR.write("{0}\t{1}\t".format(jaccard,self.anchor_threshold))
				self.ANCHOR.write("{0}\n".format(int(jaccard > self.anchor_threshold)))

				if jaccard > self.anchor_threshold:
					anyAnchorSeq = True

			os.remove(fFile)

		return anyAnchorSeq

	def prositeAnalysis(self,switch,gene,info):

		normalProtein = switch.nIsoform
		tumorProtein = switch.tIsoform
		happenings = {}

		for protein,whatsHappening in zip([normalProtein,tumorProtein],["Lost in tumor","Gained in tumor"]):
			if not protein: continue

			txFile = "{0}Data/TCGA/ProSite/{1}.out".format(options.Options().qout,protein.tx)
			if not os.path.isfile(txFile) or os.stat(txFile).st_size == 0:
				continue
			
			for line in utils.readTable(txFile,header=False):
				score = float(line[1])
				start = int(line[3][:-2])
				end = int(line[4])
				motif = line[5]

				for resNum in range(start,end):
					thisRes = protein._structure[resNum-1]
					thisRes._ptm.append(motif)

			happenings[whatsHappening] = [ x for x in protein._structure if x._ptm and x.isoformSpecific ]

		anyPTM = False
		for whatsHappening in happenings:
			for thisRes in happenings[whatsHappening]:
				self.PROSITE.write("{0}\t{1}\t{2}\t".format(info["symbol"],gene,switch.nTx))
				self.PROSITE.write("{0}\t{1}\t{2}\n".format(switch.tTx,whatsHappening,",".join(thisRes._ptm)))
				anyPTM = True

		return anyPTM

	def getThresholdFromRandom(self,sInputType,iQuantileThreshold):
		sRandomFile = "{0}structural_analysis/{1}_analysis_random.tsv".format(options.Options().qout,sInputType)
		jaccardIndex = []
		for line in utils.readTable(sRandomFile):
			jaccardIndex.append(int(line[6]))

		threshold = np.percentile(jaccardIndex,iQuantileThreshold)

		return threshold

	def recalculateScores(self,sInputType,iQuantileThreshold):
		sOldRandomFile = "{0}structural_analysis/{1}_analysis_random.tsv".format(options.Options().qout,sInputType)
		sNewRandomFile = "{0}structural_analysis/{1}_analysis_newRandom.tsv".format(options.Options().qout,sInputType)
		
		threshold = self.getThresholdFromRandom(sInputType,iQuantileThreshold)

		with open(sNewRandomFile,"w") as newRandom:
			
			newRandom.write("Symbol\tGene\tNormalTranscript\ttTumorTranscript\t")
			newRandom.write("What\tSequence\tJaccard\tThreshold\tSignificant\n")

			for line in utils.readTable(sOldRandomFile):
				jaccard = int(line[6])

				newRandom.write("{0}\t{1}\t".format(line[0],line[1]))
				newRandom.write("{0}\t{1}\t".format(line[2],line[3]))
				newRandom.write("{0}\t{1}\t".format(line[4],line[5]))
				newRandom.write("{0}\t{1}\t".format(line[6],threshold))
				newRandom.write("{0}\n".format(int(jaccard > threshold)))

		# ASEGURARSE DE QUE EL ARCHIVO ES CORRECTO
		# os.rename(sNewRandomFile,sOldRandomFile)

	def joinRandomFiles(self):

		roots = ['InterPro_report','iupred_analysis','anchor_analysis','prosite_analysis']

		for root in roots:
			files = glob.glob("{0}structural_analysis/{1}_random_*.tsv".format(options.Options().qout,root))

			outFile = "{0}structural_analysis/{1}_random.tsv".format(options.Options().qout,root)
			if os.path.exists(outFile):
				continue
			with open(outFile,"w") as OUT:
				for aFile in files:
					with open(aFile) as IN:
						IN.readline()
						OUT.write(IN.read())

		# QUITAR get thresholds from random switches
		#self.recalculateScores("anchor",95)
		#self.recalculateScores("iupred",95)

	'''
	def gpsAnalysis(self,switch,gene,info):
		normalProtein = switch.nIsoform
		tumorProtein = switch.tIsoform
		gpsFile = "{0}Data/TCGA/GPS.out".format(options.Options().wd)
		happenings = {}

		for protein,whatsHappening in zip([normalProtein,tumorProtein],["Lost in tumor","Gained in tumor"]):
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
		
if __name__ == '__main__':
	pass
