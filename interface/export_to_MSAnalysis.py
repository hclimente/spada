from libs import options
from libs import utils

class Export2MSAnalysis:
	def __init__(self):
		pass

	def generateFile(self,gn_network,tx_network):
		with open(options.Options().qout+"msInput.txt","w") as MS_INPUT:
			for gene,info,switchDict,switch in gn_network.iterate_switches_ScoreWise(tx_network,removeNoise=False,partialCreation=False):
				if switch.nIsoform and switch.tIsoform:

					nIso = switch.nIsoform.getSegments('isoform-specific')
					tIso = switch.tIsoform.getSegments('isoform-specific')

					nSp = []
					tSp = []

					[nSp.extend(x) for x in nIso]
					[tSp.extend(x) for x in tIso]

					skipping = sum([ (switch.nIsoform._structure[0]==x or switch.nIsoform._structure[-1]==x ) for x in nSp ])
					skippingTag = "S" if skipping else "NS"

					MS_INPUT.write("{0}_{1}_{2}_{3}\t".format(gene,info["symbol"],
															switch.nTx,switch.tTx))
					MS_INPUT.write("{0}\t{1}\t".format(switch.nIsoform.seq,switch.tIsoform.seq))
					MS_INPUT.write("{0}\t".format(",".join(switch.patients)))
					MS_INPUT.write("{0}\t{1}\t".format(len(nSp),len(tSp)))
					MS_INPUT.write("{0}\n".format(skippingTag))

if __name__ == '__main__':
	pass