#!/soft/devel/python-2.7/bin/python

from interface import iLoops_parser as parser
from interface import iLoops_outputPruner as pruner
from libs import options
from libs import utils
from network import gene_network, isoform_network
from methods import method

import copy
import fnmatch
import glob
import networkx as nx
import os

class CalculateInteractions(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		# remove sticky
		self._gene_network.removeLogger()
		self._modGeneNetwork = copy.deepcopy(self._gene_network)
		self._gene_network.createLogger()
		self._modGeneNetwork.createLogger()
		
		for gene,info in self._modGeneNetwork.iterate_genes_ScoreWise():
			if self._modGeneNetwork._net.degree(gene) > 20:
				self._modGeneNetwork._net.remove_node(gene)

		# get all genes at d<=2 of a driver
		d0 = set()
		d1 = set()
		d2 = set()
		for gene,info in self._modGeneNetwork.iterate_genes_ScoreWise():
			if info["Driver"]:
				d0.add(gene)
				dist1 = set(self._modGeneNetwork._net.neighbors(gene))
				d1 = d1 and dist1
				dist2 = set(nx.single_source_shortest_path_length(self._modGeneNetwork._net,gene,cutoff=2))
				d2 = d2 and dist2

		self.interestingNeighborhood = []
		self.interestingNeighborhood.extend(list(d0))
		self.interestingNeighborhood.extend(list(d1))
		# dont analyze d2 yet
		# self.interestingNeighborhood.extend(list(d2))

		# remove files from Input and logs dirs
		[ os.remove(x) for x in glob.glob("{0}Input/*".format(options.Options().qout)) ]
		[ utils.cmd("rm","-r",x) for x in glob.glob("{0}Output/*".format(options.Options().qout)) ]
		[ os.remove(x) for x in glob.glob("{0}logs/*".format(options.Options().qout)) ]
		
	def run(self):
		self.logger.info("Examining iLoops results.")

		#filetag = "_"
		filetag = "_gainedDrivers_D0_and_D1_lost_"
		self.createPartnersFastq(filetag)
		self.selectIloopsSwitches(filetag)

	def selectIloopsSwitches(self,filetag=""):

		loopFamilies = {}

		for tx,family in [ (x[0],x[1]["iLoopsFamily"]) for x in self._transcript_network.nodes(data=True) if x[1]["iLoopsFamily"] is not None ]:
			if not family in loopFamilies: 	
				loopFamilies[family] = set()
			loopFamilies[family].add(tx)

		analyzedLoops = {}

		for gene in self.interestingNeighborhood:
			gInfo = self._gene_network._net.node[gene]

			# write gene as potential interactor
			for switchDict in [ x for x in gInfo["isoformSwitches"] if not x["noise"] and x["model"]]:
				thisSwitch = self._gene_network.createSwitch(switchDict,self._transcript_network,partialCreation=True)
				nIso = thisSwitch.nTx
				tIso = thisSwitch.tTx
				nInfo = self._transcript_network._net.node[nIso]
				tInfo = self._transcript_network._net.node[tIso]

				analyze = 0
				comment = "To analyze."

				if not thisSwitch.cds_diff: 
					analyze = -1
					comment = "No CDS change."
				elif not thisSwitch.cds_overlap: 
					analyze = -1
					comment = "No overlap between CDS."
				elif not nInfo["iLoopsFamily"] or not tInfo["iLoopsFamily"]:
					analyze = -1
					comment = "No loops mapped by {0}.".format(options.Options().iLoopsVersion)
				elif nInfo["iLoopsFamily"] == tInfo["iLoopsFamily"]:
					analyze = -1
					comment = "No different loops mapped by {0}.".format(options.Options().iLoopsVersion)
				elif thisSwitch.nIsoform is None or thisSwitch.tIsoform is None:
					analyze = -1
					comment = "No UniProt mapped to at least one isoform."
				if analyze < 0:
					continue

				for isoform,thisLoopPattern in zip([nIso,tIso],[nInfo["iLoopsFamily"],tInfo["iLoopsFamily"]]):
					allProteomeOutput = "{0}iLoops/{1}/{2}/{3}.tar.gz".format(options.Options().qout,
																	   options.Options().inputType,
																	   options.Options().iLoopsVersion,
																	   isoform)
					expectedOutput = "{0}iLoops/{1}/{2}/{3}{4}.tar.gz".format(options.Options().qout,
																	   options.Options().inputType,
																	   options.Options().iLoopsVersion,
																	   isoform,filetag)

					if os.path.isfile(allProteomeOutput) or os.path.isfile(expectedOutput):
						analyze = 1
						comment = "Already analyzed."
					elif thisLoopPattern in analyzedLoops:
						analyze = 2
						if isoform in analyzedLoops[thisLoopPattern]:
							comment = "Already being analyzed in this batch."
						else:
							comment = "Analyzing relative {0}.".format(analyzedLoops[thisLoopPattern])
					else:
						for iso in loopFamilies[thisLoopPattern]:
							if os.path.isfile("{0}iLoops/TCGA/{1}/{2}.tar.gz".format(options.Options().wd, options.Options().iLoopsVersion,iso) ):
								analyze = 2
								comment = "Analyzed relative {0}.".format(iso)
								break

					self.launchIloops(gene,thisSwitch,filetag)

					if analyze == 0: 
						analyzedLoops[thisLoopPattern] = isoform

	def createPartnersFastq(self,filetag="",gene=""):
				
		familiesAdded = set()
		wannaWrite = False
		fileCounter = 1
		transcriptCounter = 0

		candidatePartners = []
		if gene:
			self.logger.info("Preparing FASTA files for all interactors of {0}.".format(gene))
			basename = "{0}candidateLostPartners_{1}".format(options.Options().qout,gene)
			candidatePartners = [ x for x,y in self._transcript_network.nodes(data=True) if y["gene_id"] in self._gene_network._net.neighbors(gene) ]
			if not candidatePartners:
				return False
		elif filetag:
			self.logger.info("Preparing FASTA files for all interesting partners.")
			basename = "{0}candidateGainPartners_".format(options.Options().qout)
			candidatePartners = [ x for x,y in self._transcript_network.nodes(data=True) if y["gene_id"] in self.interestingNeighborhood ]
			if not candidatePartners:
				return False
		else:
			self.logger.info("Preparing FASTA files for all genes.")
			basename = "{0}allProteome_".format(options.Options().qout)

		[ os.remove(x) for x in glob.glob("{0}*".format(basename)) ]

		MULTIFASTA = open("{0}{1}.fasta".format(basename,fileCounter), "w")

		with open("Data/{0}/UnifiedFasta_{1}.fa".format(options.Options().inputType,options.Options().iLoopsVersion),"r") as gcMULTIFASTA:
			for line in gcMULTIFASTA:
				if ">" in line:
					wannaWrite = False
					loopFamily = line.strip().split("#")[3]
					identifier = line[1:].strip().split("#")[0]

					if transcriptCounter == 1500:
						fileCounter += 1
						MULTIFASTA.close()
						MULTIFASTA = open("{0}{1}.fasta".format(basename,fileCounter), "w")
						transcriptCounter = 0

					if loopFamily:
						if candidatePartners:
							if identifier in candidatePartners:
								wannaWrite = True
								transcriptCounter += 1
								MULTIFASTA.write(">" + identifier + "\n")
						elif loopFamily not in familiesAdded:
							familiesAdded.add(loopFamily)
							wannaWrite = True
							transcriptCounter += 1
							MULTIFASTA.write(">" + identifier + "\n")					

				elif wannaWrite:
					MULTIFASTA.write(line)

	def launchIloops(self,gene,thisSwitch,filetag):
		
		self.createPartnersFastq(gene=gene)

		gainedInteractions = glob.glob("{0}candidateGainPartners_*.fasta".format(options.Options().qout))
		lostInteractions = glob.glob("{0}candidateLostPartners_{1}_*.fasta".format(options.Options().qout,gene))
		allProteome = glob.glob("{0}allProteome_*.fasta".format(options.Options().qout))

		partnersToTest = []
		partnersToTest.extend(gainedInteractions)
		partnersToTest.extend(lostInteractions)
		partnersToTest.extend(allProteome)

		batch = 0
		for tx,seq in zip([thisSwitch.nTx,thisSwitch.tTx],[thisSwitch.nIsoform.seq,thisSwitch.tIsoform.seq]):
		 	for fastaFile in partnersToTest:
		 		batch += 1
		 		tag = "{0}_{1}".format(tx,batch)
		 		self.getFinalFASTAandPairs(fastaFile,tx,seq,batch)
		 		self.logger.info("Launching iLoops {0}/{1} for switch {2},{3}".format(batch,len(partnersToTest),thisSwitch.nTx,thisSwitch.tTx))
				utils.cmd("/soft/devel/python-2.7/bin/python",
						  "/sbi/programs/{0}/iLoops.py".format(options.Options().iLoopsVersion),
		 				  "-f {0}Input/{1}.fasta".format(options.Options().qout,tag),
		 				  "-q {0}Input/{1}.net".format(options.Options().qout,tag),
		 				  "-j {0}Output/{1}".format(options.Options().qout,tag),
		 				  "-x {0}.xml".format(tag),
		 				  "-v",
		 				  "-g all",
		 				  "-n 25",
		 				  "-Q sbi",
		 				  "-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50",
		 				  "2>&1 >{0}logs/{1}.log".format(options.Options().qout,tag) )

			task = "/sbi/users/hectorc/SmartAS_experimental/Pipeline/interface/iLoops_outputPruner.py "
			task += "{0} {1}Output/ {2} ".format(tx,options.Options().qout,filetag)
			task += "{0} {1}".format(options.Options().inputType,options.Options().iLoopsVersion)
			utils.launchSingleJob(task,"iLoops_{0}{1}".format(tx,filetag))

	def getFinalFASTAandPairs(self,fastaFile,tx,seq,batch):
		tag = "{0}_{1}".format(tx,batch)
		utils.cmd("cp",fastaFile,"{0}Input/{1}.fasta".format(options.Options().qout,tag))

	 	with open("{0}Input/{1}.fasta".format(options.Options().qout,tag),"r") as SEQS, \
	 		 open("{0}Input/{1}.net".format(options.Options().qout,tag),"w") as PAIRS:
	 		for line in SEQS:
	 			if ">" in line:
	 				partner = line[1:].strip()
	 				PAIRS.write("{0}\t{1}\n".format(tx,partner))

	 	with open("{0}Input/{1}.fasta".format(options.Options().qout,tag),"a") as SEQS:
	 		SEQS.write(">{0}\n{1}\n".format(tx,seq))

if __name__ == '__main__':
	pass
