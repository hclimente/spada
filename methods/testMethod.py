from libs import options
from libs import utils
from methods import method

import cPickle
import pdb

class Test(method.Method):
	def __init__(self,gn_network,tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		import glob
		files = glob.glob("{0}structural_analysis/structural_summary_random_*.tsv".format(options.Options().qout))

		for afile in files:
			name=afile.split("/")[-1]
			idx=name.split("_")[3]

			newfile = "{0}structural_analysis/structural_summary_newrandom_{1}".format(options.Options().qout,idx)
			with open(newfile,"w") as OUT:
				OUT.write("Gene\tSymbol\tNormalTranscript\tTumorTranscript\t")
				OUT.write("iLoops\tDomain\tDisorder\tAnchor\tPTM\n")
				for line in utils.readTable(afile):
					symbol = self._gene_network._net.node[line[0]]["symbol"]
						
					OUT.write("{0}\t{1}\t".format(line[0],symbol))
					OUT.write("{0}\t{1}\t".format(line[1],line[2]))
					OUT.write("{0}\t{1}\t".format(line[3],line[4]))
					OUT.write("{0}\t{1}\t".format(line[5],line[5]))
					OUT.write("{0}\n".format(line[7]))

			os.rename(afile,afile+'_old')
			os.rename(newfile,afile)
		