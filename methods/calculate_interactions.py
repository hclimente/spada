#!/soft/devel/python-2.7/bin/python

from interface import iLoops_parser as parser
from interface import iLoops_outputPruner as pruner
from libs import options
from libs import utils
from network import gene_network, isoform_network
from methods import method

import copy
import fnmatch
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

	def run(self):
		self.logger.info("Examining iLoops results.")
		self.selectIloopsSwitches("Drivers")

	def selectIloopsSwitches(self,filetag=""):

		loopFamilies = {}

		for tx,family in [ (x[0],x[1]["iLoopsFamily"]) for x in self._transcript_network.nodes(data=True) if x[1]["iLoopsFamily"] is not None ]:
			if not family in loopFamilies: 	
				loopFamilies[family] = set()
			loopFamilies[family].add(tx)

		analyzedLoops = {}

		for gene in self.interestingNeighborhood:
			gInfo = self._modGeneNetwork._net.node[gene]

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

				if analyze < 0:
					continue

				allProteomeOutput = "{0}iLoops/{1}/{2}/{3}.tar.gz".format(options.Options().qout,
																	   options.Options().inputType,
																	   options.Options().iLoopsVersion,
																	   isoform)
				expectedOutput = "{0}iLoops/{1}/{2}/{3}{4}.tar.gz".format(options.Options().qout,
																	   options.Options().inputType,
																	   options.Options().iLoopsVersion,
																	   isoform,filetag)

				for isoform,thisLoopPattern in zip([nIso,tIso],[nInfo["iLoopsFamily"],tInfo["iLoopsFamily"]]):
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

					self.analyzeSwitch(gene,thisSwitch,filetag)

					if analyze == 0: 
						analyzedLoops[thisLoopPattern] = isoform

		def analyzeSwitch(self,gene,thisSwitch,filetag):
			self.createPartnersFastq(thisSwitch,filetag)
			self.launchIloops(thisSwitch)

		def createPartnersFastq(self,gene,thisSwitch,filetag=""):
			self.logger.info("Preparing FASTA files for all transcripts.")
			basename = "{0}testedPartners_".format(options.Options().qout)

			wannaWrite = False
			fileCounter = 1
			transcriptCounter = 0

			familiesAdded = set()

			candidatePartners = []
			if filetag:
				candidatePartners.extend([ x for x,y in self._transcript_network.nodes(data=True) if y["gene_id"] in self._modGeneNetwork._net.neighbors(gene) ])

			MULTIFASTA = open("{0}{1}{2}.fasta".format(options.Options().qout,basename,fileCounter), "w")

			with open("Data/{0}/UnifiedFasta_{1}.fa".format(inputType,iLoopsVersion),"r") as gcMULTIFASTA:
				for line in gcMULTIFASTA:
					if ">" in line:
						wannaWrite = False
						loopFamily = line.strip().split("#")[3]
						identifier = line[1:].strip().split("#")[0]

						if transcriptCounter >= 1499:
							fileCounter += 1
							MULTIFASTA.close()
							MULTIFASTA = open("{0}{1}{2}.fasta".format(options.Options().qout,basename,fileCounter), "w")
							transcriptCounter = 0
						
						if not candidatePartners:
							if loopFamily and loopFamily not in familiesAdded:
								familiesAdded.add(loopFamily)
								wannaWrite = True
								transcriptCounter += 1
								MULTIFASTA.write(">" + identifier + "\n")
						elif identifier in candidatePartners:
							familiesAdded.add(loopFamily)
							wannaWrite = True
							transcriptCounter += 1
							MULTIFASTA.write(">" + identifier + "\n")
							
					elif wannaWrite:
						MULTIFASTA.write(line)

	def launchIloops(self,thisSwitch):
		
		self.logger.info("Launching iLoops.")

		for tx,seq in zip([thisSwitch.nTx,thisSwitch.tTx],[thisSwitch.nIsoform.seq,thisSwitch.tIsoform.seq]):
		 	for fastaFile in fnmatch.filter(os.listdir(options.Options().qout), "testedPartners_*.fasta"):
		 		batch = (fastaFile.split(".")[0]).split("_")[1]
		 		tag = tx + "_" + batch
		 		self.getFinalFASTAandPairs(tx,seq,batch)
		 			
				utils.cmd(	"/soft/devel/python-2.7/bin/python",
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

			p = pruner.iLoopsOutput_pruner(tx, options.Options().qout + "Output/")
			p.joinFiles()
			if p.makeLiteVersion():
				utils.cmd("scp","-r", 
					"{0}Output/{1}.tar.gz".format(options.Options().qout,tx), 
					"hector@gencluster:~/iLoops/{0}/{1}".format(
						options.Options().inputType,
						options.Options().iLoopsVersion)
					)
			else:
				self.logger.error("Error in generation of file.")

	def getFinalFASTAandPairs(self,tx,seq,batch):
		tag = tx + "_" + batch
		utils.cmd("cp", 
			"{0}testedPartners_{1}.fasta".format(options.Options().qout,batch), 
			"{0}Input/{1}.fasta".format(options.Options().qout,tag))

	 	with open("{0}Input/{1}.fasta".format(options.Options().qout,tag), "r") as SEQS, \
	 		 open("{0}Input/{1}.net".format(options.Options().qout,tag), "w") as PAIRS:
	 		for line in SEQS:
	 			if ">" in line:
	 				partner = line[1:].strip()
	 				PAIRS.write("{0}\t{1}\n".format(tx,partner))

	 	with open("{0}Input/{1}.fasta".format(options.Options().qout,tag),"a") as SEQS:
	 		SEQS.write(">{0}\n{1}\n".format(tx,seq))

if __name__ == '__main__':
	pass