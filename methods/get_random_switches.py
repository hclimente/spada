from biological_entities import switch
from libs import options
from libs import utils
from methods import method
from methods import structural_analysis

import copy
import cPickle
import itertools
import random
import os

class GetRandomSwitches(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)

		self.MAX_SWITCHES = 5

	def run(self):
		if not os.path.exists("{0}randomGeneNetwork.pkl".format(options.Options().qout)):
			self.logger.info("Generating random switches.")

			# copy original gene network
			self._gene_network.removeLogger()
			gnNetCopy = copy.deepcopy(self._gene_network)
			self._gene_network.createLogger()
			gnNetCopy.createLogger()

			# remove real switches and calculate new ones
			self.sampleSwitches(gnNetCopy)
			gnNetCopy.saveNetwork("randomGeneNetwork.pkl")

			utils.launchJobs(self._gene_network,"random")

		else:
			# calculate relevant switches
			self.logger.info("Calculating features for random switches.")
			gnNetCopy = cPickle.load(open("{0}randomGeneNetwork.pkl".format(options.Options().qout)))
			S = structural_analysis.StructuralAnalysis(gnNetCopy,self._transcript_network,isRand=True)
			S.run()	

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
				switchDict["noise"] = False
				switchDict["model"] = True

				info["isoformSwitches"].append(switchDict)

	def launchJobs(self):
		import drmaa

		s = drmaa.Session()
		s.initialize()

		natSpec = ""
		natSpec += "-q short,normal -l 'qname=short|normal' "
		natSpec += "-cwd "
		natSpec += "-V "

		randoms = []

		for startingNode in range(1,len(self._gene_network.nodes()),20):
			cfg = options.Options().printToFile(filename="random_{0}_node{1}".format(options.Options().tag,startingNode),parallelRange=startingNode,onlyModels=False)
			jt = s.createJobTemplate()
			
			jt.remoteCommand = 'Pipeline/gSmartAS.py'
			jt.args = ['-f',cfg]
			jt.joinFiles=True
			jt.nativeSpecification = natSpec
			jt.nativeSpecification += "-N {0}_random_{1} ".format(options.Options().tag,startingNode)
			jt.nativeSpecification += "-e {0}logs/random_{1}.txt ".format(options.Options().qout,startingNode)
			jt.nativeSpecification += "-o {0}logs/random_{1}.txt".format(options.Options().qout,startingNode)

			if not os.path.exists("{0}logs".format(options.Options().qout)):
			    os.makedirs("{0}logs".format(options.Options().qout))

			jobid = s.runJob(jt)
			randoms.append(jobid)
			s.deleteJobTemplate(jt)


if __name__ == '__main__':
	pass
