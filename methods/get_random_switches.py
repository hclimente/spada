from biological_entities import switch
from libs import options
from libs import utils
from methods import method
from methods import structural_analysis

import copy
import cPickle as pickle
import itertools
import operator
import random
import os

class GetRandomSwitches(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)

		self.MAX_SWITCHES = 5

	def run(self):
		if not os.path.exists("{}randomGeneNetwork.pkl".format(options.Options().qout)):
			self.logger.info("Generating random switches.")

			# copy original gene network and remove real switches
			self._gene_network.removeLogger()
			gcopy = copy.deepcopy(self._gene_network)
			self._gene_network.createLogger()
			gcopy.createLogger()
			gcopy.cleanNetwork()
			gcopy.removeLogger()
			
			# random: two random transcripts among those expressed in tumor or normal
			gcopy_random = copy.deepcopy(gcopy)
			gcopy_random.createLogger()
			self.sampleTranscripts_random(gcopy_random)
			gcopy_random.saveNetwork("randomGeneNetwork_random.pkl")			

			# random: most expressed isoform is normal
			gcopy_fixNormal = copy.deepcopy(gcopy)
			gcopy_fixNormal.createLogger()
			self.sampleTranscripts_fixNormal(gcopy_fixNormal)
			gcopy_fixNormal.saveNetwork("randomGeneNetwork_fixNormal.pkl")

			utils.launchJobs(self._gene_network,"random")

		else:
			# calculate functional switches
			self.logger.info("Calculating features for random switches.")
			gcopy = pickle.load(open("{}randomGeneNetwork.pkl".format(options.Options().qout)))
			S = structural_analysis.StructuralAnalysis(gcopy,self._transcript_network,isRand=True)
			S.run()	

	def sampleTranscripts_fixNormal(self,gcopy):
		for gene,info in [ (x,y) for x,y in gcopy.nodes(data=True) if len(set(y["expressedTxsNormal"]) & set(y["expressedTxsTumor"]))>1]:
			# set normal isoform as the most expressed in normal, shuffle the rest for tumor
			txExpression = [ (x,self._transcript_network._net.node[x]["median_TPM_N"]) for x in info["expressedTxsNormal"] ]
			nTx = max(txExpression, key=operator.itemgetter(1))[0]
			nTxExpression = max(txExpression, key=operator.itemgetter(1))[1]

			if nTx not in info["expressedTxsNormal"]:
				self.logger.warning("Median most expressed transcript from gene {} is not considered expressed. \
					Probably due to the threshold applied. TPM={}. Will be skipped. ".format(gene,txExpression))
				continue

			txs = list(info["expressedTxsTumor"])
			if nTx in txs:
				txs.remove(nTx)

			numSwitches = self.MAX_SWITCHES
			if len(txs) < self.MAX_SWITCHES:
				numSwitches = len(txs)

			random.shuffle(txs)

			for i in range(numSwitches):
				switchDict = {}
				switchDict["nIso"] = nTx
				switchDict["tIso"] = txs[i]

				switchDict["patients"] = []
				switchDict["noise"] = False
				switchDict["functional"] = None

				info["isoformSwitches"].append(switchDict)

	def sampleTranscripts_random(self, gcopy):
		
		for gene,info in [ (x,y) for x,y in gcopy.nodes(data=True) if len(set(y["expressedTxsNormal"]) & set(y["expressedTxsTumor"]))>1]:
			allExpressedTxs = set(y["expressedTxsNormal"]) & set(y["expressedTxsTumor"])

			allSwitches = [ x for x in itertools.combinations(allExpressedTxs,2) ]
			random.shuffle(allSwitches)
			
			allSwitches = allSwitches[0:self.MAX_SWITCHES]
			
			for oneSwitch in allSwitches:
				switchDict = {}
				switchDict["nIso"] = oneSwitch[0]
				switchDict["tIso"] = oneSwitch[1]
				switchDict["patients"] = []
				switchDict["noise"] = False
				switchDict["model"] = True
				switchDict["functional"] = None
				info["isoformSwitches"].append(switchDict)

if __name__ == '__main__':
	pass
