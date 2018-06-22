import gzip
import pickle
import logging
import os.path

class Method:
	def __init__(self, name, annotation):

		self.logger = logging.getLogger(name)

		if isinstance(annotation, str):
			self._genes,self._txs = self.loadNetworks(annotation)
			self._genes.createLogger()
			self._txs.createLogger()
		elif isinstance(annotation, tuple):
			self._genes,self._txs = annotation
		else:
			self._genes,self._txs = None,None

	def saveNetworks(self, filename = 'annotation.pklz'):

		self.logger.debug("Saving annotation at {}.".format(filename))
		# unattach logger to save without thread problems
		self._genes.removeLogger()
		self._txs.removeLogger()

		pickle.dump((self._genes, self._txs), gzip.open(filename, "wb"), -1)

		self._genes.createLogger()
		self._txs.createLogger()

	def loadNetworks(self, annotation):
		return pickle.load(gzip.open(annotation, "rb"))
