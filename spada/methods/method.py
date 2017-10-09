from libs import options

import pickle
import logging

class Method:
	def __init__(self, name, gn_network,tx_network,gn_subnetwork=False):
		self.logger = logging.getLogger(name)

		if gn_network:
			if isinstance(gn_network,bool):
				self._gene_network = pickle.load(open(options.Options().qout + "geneNetwork.pkl","rb"))
			elif isinstance(gn_network,str):
				self._gene_network = pickle.load(open(options.Options().qout+gn_network,"rb"))
			else:
				self._gene_network = gn_network
			self._gene_network.createLogger()

		if tx_network:
			if isinstance(tx_network,bool):
				self._transcript_network = pickle.load(open(options.Options().qout + "txNetwork.pkl","rb"))
			elif isinstance(tx_network,str):
				self._transcript_network = pickle.load(open(options.Options().qout+tx_network,"rb"))
			else:
				self._transcript_network = tx_network
			self._transcript_network.createLogger()

		if gn_subnetwork:
			if isinstance(gn_subnetwork,bool):
				self._gene_subnetwork = pickle.load(open(options.Options().qout + "geneSubnetwork.pkl","rb"))
			elif isinstance(gn_subnetwork,str):
				self._gene_subnetwork = pickle.load(open(options.Options().qout+gn_subnetwork,"rb"))
			else:
				self._gene_subnetwork = gn_subnetwork
			self._gene_subnetwork.createLogger()
