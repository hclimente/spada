import pandas as pd
import numpy as np
from libs import utils as ut
import network

class IsoformNetwork(network.Network):
	def __init__(self):
		network.Network.__init__(self)

	def add_node(self, isoform_id, mTPM_N=0.0, mPSI_N=0.0, mTPM_T=0.0, mPSI_T=0.0, gene_full_name=""):
		debug = False
		
		nameComponents 	= gene_full_name.split("|")
		if len(nameComponents) == 1:
			if debug: print(gene_full_name + " not added. Unable to get gene id and symbol.")
			return
		
		geneSymbol 	= nameComponents[0]
		geneID 		= nameComponents[1]

		if "locus" in gene_full_name:
			if debug: print( "Unknown gene " + gene_full_name + " not added.")
			return
		elif "locus" in gene_full_name:
			if debug: print( "Unknown gene " + gene_full_name + " not added.")
			return

		self._net.add_node( isoform_id, median_TPM_N=mTPM_N, median_PSI_N=mPSI_N, 
							median_TPM_T=mTPM_T, median_PSI_T=mPSI_T, 
							gene_id=geneID, gene_symbol=geneSymbol
						  )

	def update_node(self, id, properties):
		if id in self._net.nodes():
			for p in properties.keys(): self._net.node[id][p] = properties[p]

	def add_edge(self, isoform_id1, isoform_id2, Weight=0.0, iLoops=False, exp=False):
		self._net.add_node( isoform_id1, isoform_id2, weight=Weight, iLoops_prediction=iLoops, experiment=exp )

	def importTranscriptome(self, path):

		for line in ut.readTable(path + "expression_normal.tsv"):
			txName 			= line[1]
			gene_full_name 	= line[0]
			if line[3] == "NA":	median_PSI 	= 0.0
			else: 				median_PSI 	= float(line[3])

			if line[5] == "NA": median_TPM 	= 0.0
			else: 				median_TPM 	= float(line[5])

			self.add_node(txName, mTPM_N=median_TPM, mPSI_N=median_PSI, gene_full_name=gene_full_name)

		for line in ut.readTable(path + "expression_tumor.tsv"):
			txName 			= line[1]
			if line[3] ==  "NA":median_PSI 	= 0.0
			else: 				median_PSI 	= float(line[3])

			if line[5] == "NA": median_TPM 	= 0.0
			else: 				median_TPM 	= float(line[5])
				
			self.update_node(txName, {"median_TPM_T": median_TPM, "median_PSI_T" : median_PSI } )