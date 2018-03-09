from spada import utils
from spada.network import network

import pandas as pd
import numpy as np
import abc

class TranscriptNetwork(network.Network):
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
		median_TPM_N(float,None)	Median TPM of the isoform in the normal patients.
		median_PSI_N(float,None)	Median PSI of the isoform in the normal patients.
		median_TPM_T(float,None)	Median TPM of the isoform in the tumor patients.
		median_PSI_T(float,None)	Median PSI of the isoform in the tumor patients.
		proteinSequence(str,None)	Protein sequence.

	Edge information:
		Id1(str)					Transcript id of interactor 1.
		Id2(str)					Transcript id of interactor 2.
	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self, name):
		network.Network.__init__(self, name)

	@abc.abstractmethod
	def acceptCDS(self, **kwds):
		raise NotImplementedError()

	@abc.abstractmethod
	def genenameFilter(self, **kwds):
		raise NotImplementedError()

	def add_node(self,tx,geneFullname):

		if tx in self.nodes():
			self.logger.debug("Transcript {}, gene {} already in the network.".format(tx, geneFullname))
			return True

		gene = self.genenameFilter(full_name=geneFullname)[0]

		if not gene:
			self.logger.debug("No gene could be extracted for transcript {}, gene {}.".format(tx, geneFullname))
			return False

		self._net.add_node( tx,
							gene_id			= gene,
							exons	= [],
							txCoords		= None,
							CDS		= None,
							strand			= None,
							chr				= None,
							median_TPM_N	= None,
							median_PSI_N	= None,
							median_TPM_T	= None,
							median_PSI_T	= None,
							proteinSequence	= None,
							Pfam			= {},
							Prosite			= {},
							IDR				= {})

		return True

	def update_node(self, tx, key, value, secondKey = ""):

		override = False

		# CDS
		if key == 'CDS' and self.nodes(data=True)[tx][key]:
			override = True
			i = self._net.node[tx]['strand'] == '-'
			value[i] = self._net.node[tx][key][i]

		return self._update_node(tx, key, value, secondKey, override)

	def update_nodes(self, key, values):
		for tx, value in values.items():
			if isinstance(value, set):
				for v in value:
					self.update_node(tx, key, v)
			else:
				self.update_node(tx, key, value)

	def add_edge(self, tx1, tx2, **kwargs):
		self._net.add_edge(tx1, tx2, **kwargs)

	def update_edge(self, tx1, tx2, key, value):
		return self._update_edge(tx1, tx2, key, value)

	def transcripts(self):
		'''
		Iterate transcripts.
		'''

		for tx,info in self.nodes(data=True):
			yield tx, info