import pickle
import logging

class Method:
	def __init__(self, name, genes, transcripts):
		self.logger = logging.getLogger(name)

		if genes:
			if isinstance(genes, bool):
				self._genes = pickle.load(open("genes.pkl", "rb"))
			elif isinstance(genes, str):
				self._genes = pickle.load(open(genes, "rb"))

			self._genes.createLogger()

		if transcripts:
			if isinstance(transcripts, bool):
				self._txs = pickle.load(open("transcripts.pkl", "rb"))
			elif isinstance(transcripts, str):
				self._txs = pickle.load(open(transcripts, "rb"))

			self._txs.createLogger()
