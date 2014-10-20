from libs import options
from libs import utils

class Export2MSAnalysis:
	def __init__(self):
		pass

	def generateFile(self,gn_network):
		with open(options.Options().qout+"msInput.txt","w") as MS_INPUT:
			for gene,info,switch in utils.iterate_switches_ScoreWise(gn_network):
				if switch.nIsoform and switch.tIsoform:
					MS_INPUT.write("{0}_{1}_{2}_{3}\t".format(gene,info["symbol"],
															switch.nTx,switch.tTx))
					MS_INPUT.write("{0}\t{1}\t".format(switch.nIsoform.seq,switch.tIsoform.seq))
					MS_INPUT.write("{0}\n".format(",".join(switch.patients)))

if __name__ == '__main__':
	pass