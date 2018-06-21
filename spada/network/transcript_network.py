from spada.io import io
from spada.network.network import Network

import numpy as np
import abc

class TranscriptNetwork(Network):
	"""docstring for TranscriptNetwork
	TranscriptNetwork contains a network of isoforms.

	Node information:
		id(str)						Transcript Id
		gene_id(str)				Gene Id of the parent gene.
		exons(list,None)			List of lists, each of them containing the limits of an exon.
		txCoords(list,None)			List with the starting and the ending genome positions
									of the trancript.
		CDS(list,None)				Starting and the ending genome positions of the CDS.
		strand(str,None)			Strand.
		chr(str,None)				Chromosome.
		main(bool,False)			Is it the main transcript of the gene?
		start_codon(int,None)		First residue of the start codon.
		stop_codon(int,None)		First residue of the stop codon.
		proteinSequence(str,None)	Protein sequence.

	Edge information:
		Id1(str)					Transcript id of interactor 1.
		Id2(str)					Transcript id of interactor 2.
	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):
		Network.__init__(self, name)

	@abc.abstractmethod
	def acceptCDS(self, **kwds):
		raise NotImplementedError()

	@abc.abstractmethod
	def genenameFilter(self, **kwds):
		raise NotImplementedError()

	@abc.abstractmethod
	def isMain(self, **kwds):
		raise NotImplementedError()

	def add_node(self, tx_name, gene_full_name):

		if tx_name in self.nodes():
			self.logger.debug("Transcript {}, gene {} already in the network.".format(tx_name,  gene_full_name))
			return True

		gene = self.genenameFilter(full_name = gene_full_name)[0]

		if not gene:
			self.logger.debug("No gene could be extracted for transcript {}, gene {}.".format(tx_name,  gene_full_name))
			return False

		self._net.add_node( tx_name,
							gene_id			= gene,
							exons			= [],
							txCoords		= None,
							CDS				= None,
							strand			= None,
							chr				= None,
							main 			= False,
							start_codon		= None,
							stop_codon		= None,
							proteinSequence	= None,
							Pfam			= {},
							Prosite			= {},
							IDR				= {})

		return True

	def update_node(self, tx_name, key, value, secondKey = ""):

		override = False

		# CDS
		if key == 'CDS' and self.nodes(data=True)[tx_name][key]:
			override = True
			i = self._net.node[tx_name]['strand'] == '-'
			value[i] = self._net.node[tx_name][key][i]
		elif key in ['start_codon','stop_codon'] and self.nodes(data=True)[tx_name][key]:
			override = True
			old = self._net.node[tx_name][key]
			value = min(old, value) if self._net.node[tx_name]['strand'] == '+' else max(old, value)

		return self._update_node(tx_name, key, value, secondKey, override)

	def add_edge(self, tx1, tx2, **kwargs):
		self._net.add_edge(tx1, tx2, **kwargs)

	def transcripts(self, onlyMain = False):
		'''
		Iterate transcripts.
		'''

		for tx,info in self.nodes(data=True):
			if onlyMain and not info['main']:
				continue
			yield tx, info
