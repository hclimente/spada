from libs import options
from libs import utils

class Export2MSAnalysis:
	def __init__(self):
		pass

	def generateFile(self,gn_network,tx_network):
		with open(options.Options().qout+"msInput.txt","w") as MS_INPUT:
			for gene,info,switchDict,switch in gn_network.iterate_switches_ScoreWise(tx_network,removeNoise=False,partialCreation=False):
				if switch.nIsoform and switch.tIsoform:

					nIso = [ len(x) for x in switch.nIsoform.getSegments('isoform-specific') ]
					nIso = sum(nIso)
					tIso = [ len(x) for x in switch.tIsoform.getSegments('isoform-specific') ]
					tIso = sum(tIso)

					MS_INPUT.write("{0}_{1}_{2}_{3}\t".format(gene,info["symbol"],
															switch.nTx,switch.tTx))
					MS_INPUT.write("{0}\t{1}\t".format(switch.nIsoform.seq,switch.tIsoform.seq))
					MS_INPUT.write("{0}\t".format(",".join(switch.patients)))
					MS_INPUT.write("{0}\t{1}\n".format(nIso,tIso))

if __name__ == '__main__':
	pass