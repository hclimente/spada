from biological_entities import switch
from libs import options
from libs import utils
from methods import method
from methods import structural_analysis

import copy
import itertools
import random

class GetRandomSwitches(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)

		self.MAX_SWITCHES = 10

	def run(self):
		self.logger.info("Generating random switches.")

		# copy original gene network
		self._gene_network.removeLogger()
		gnNetCopy = copy.deepcopy(self._gene_network)
		self._gene_network.createLogger()
		gnNetCopy.createLogger()

		# remove real switches and calculate new ones
		self.sampleSwitches(gnNetCopy)

		# calculate relevant switches
		S = structural_analysis.StructuralAnalysis(gnNetCopy,self._transcript_network,isRand=True)
		S.run()

		gnNetCopy.saveNetwork("randomGeneNetwork.pkl")

	def sampleSwitches(self,gnNetCopy):
		# remove real switches
		for gene,info in gnNetCopy.nodes(data=True):
			info["isoformSwitches"] = []

		# create new ones
		for gene,info in [ (x,y) for x,y in gnNetCopy.nodes(data=True) if len(y["ExpressedTranscripts"])>1 ]:
			allSwitches = [ x for x in itertools.combinations(info["ExpressedTranscripts"],2) ]
			random.shuffle(allSwitches)
			allSwitches = allSwitches[0:self.MAX_SWITCHES]

			for oneSwitch in allSwitches:
				switchDict = {}
				switchDict["nIso"] = oneSwitch[0]
				switchDict["tIso"] = oneSwitch[1]

				switchDict["score"] = 0.0
				switchDict["patients"] = 0.0
				switchDict["precision"] = 0.0
				switchDict["sensitivity"] = 0.0

				info["isoformSwitches"].append(switchDict)

if __name__ == '__main__':
	pass