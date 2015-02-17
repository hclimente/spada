#!/soft/devel/python-2.7/bin/python

from libs import options

import cPickle
import gridmap
import logging

class Method:
	def __init__(self, name, gn_network,tx_network,gn_subnetwork=False):
		self.logger = logging.getLogger(name)

		if gn_network:
			if isinstance(gn_network,bool):
				self._gene_network = cPickle.load(open(options.Options().qout + "geneNetwork.pkl","r"))
			elif isinstance(gn_network,str):
				self._gene_network = cPickle.load(open(options.Options().qout+gn_network,"r"))
			else:
				self._gene_network = gn_network
			self._gene_network.createLogger()
		
		if tx_network:
			if isinstance(tx_network,bool):
				self._transcript_network = cPickle.load(open(options.Options().qout + "txNetwork.pkl","r"))
			elif isinstance(tx_network,str):
				self._transcript_network = cPickle.load(open(options.Options().qout+tx_network,"r"))
			else:
				self._transcript_network = tx_network
			self._transcript_network.createLogger()
			
		if gn_subnetwork:
			if isinstance(gn_subnetwork,bool):
				self._gene_subnetwork = cPickle.load(open(options.Options().qout + "geneSubnetwork.pkl","r"))
			elif isinstance(gn_subnetwork,str):
				self._gene_subnetwork = cPickle.load(open(options.Options().qout+gn_subnetwork,"r"))
			else:
				self._gene_subnetwork = gn_subnetwork
			self._gene_subnetwork.createLogger()

	def grid(self,f,*args):
		arguments = [ x for x in args ]
		funct = gridmap.Job(f,arguments,queue="normal")
		out = gridmap.process_jobs([funct])
		return out