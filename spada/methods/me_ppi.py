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

		self.logger.info("Calculating proportion of pannegative.")

		utils.cmd('Rscript', 'pipeline/methods/me_ppi.R',
				  options.Options().wd, options.Options().qout)

if __name__ == '__main__':
	pass