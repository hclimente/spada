from libs import options
from libs import utils
from methods import method

import glob
import time

class MEPPI(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)

	def clean(self):
		utils.cmd("rm","-r","{}ppi".format(options.Options().qout))
		utils.cmd("mkdir","-p","{}ppi".format(options.Options().qout))

	def run(self):

		self.logger.info("Preparing the environment variables.")

		utils.cmd('Rscript', 'pipeline/methods/me_ppi.environment.R',
				  options.Options().wd, options.Options().qout)

		self.logger.info("Launching jobs.")
		i = 0
		for symbol in utils.readTable("{}ppi/genes".format(options.Options().qout)):
			symbol = symbol[0]
			utils.geneclusterLaunch("{}-ppi".format(symbol), "~/anaconda/bin/Rscript",
				'pipeline/methods/me_ppi.R', options.Options().wd, options.Options().qout, symbol)

			i += 1

			if i%1000 == 0:
				time.sleep(200)

		self.logger.info("Waiting for the all the threads to finish.")
		time.sleep(1800)

		files = glob.glob("{}ppi/*.txt".format(options.Options().qout))
		outFile = "{}ppi/me_ppi.tsv".format(options.Options().qout)

		with open(outFile,"w") as OUT:
			OUT.write("Gene\tfisher.p\tsampling.p\tn.patients\tn.genes\tn.drivers\n")
			for aFile in files:
				with open(aFile) as IN:
					IN.readline()
					OUT.write(IN.read())

if __name__ == '__main__':
	pass