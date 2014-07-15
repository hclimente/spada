import pandas as pd
import numpy as np
from libs import utils as ut
from libs import options
import network

class IsoformNetwork(network.Network):
	"""docstring for IsoformNetwork
	IsoformNetwork contains a network of isoforms.

	Node information:
		Id - str Transcript Id
		median_TPM_N - float Median TPM of the isoform in the normal patients.
		median_PSI_N - float Median PSI of the isoform in the normal patients.
		median_TPM_T - float Median TPM of the isoform in the tumor patients.
		median_PSI_T - float Median PSI of the isoform in the tumor patients.
		iLoopsFamily - str Loop pattern.
		gene_id - str Gene Id of the parent gene.

	Edge information:
		Id1 - str Transcript id of interactor 1.
		Id2 - str Transcript id of interactor 2.
		iLoops_prediction - bool, None iLoops predicted interaction.
		RC - float,None iLoops max RC with a prediction.
		experiment - bool,None Interaction from an experiment.
						  )
	"""
	def __init__(self):
		network.Network.__init__(self)

	def add_node(self, tx, gene_full_name):
		debug = True
		
		nameComponents 	= gene_full_name.split("|")
		if len(nameComponents) == 1:
			if debug: print("{0} not added. Unable to get gene id and symbol from {1}.".format(tx, gene_full_name))
			return False
		
		geneID = nameComponents[1]

		if "locus" in gene_full_name:
			if debug: print( "Unknown gene {0} not added.".format(gene_full_name))
			return False
		elif "locus" in gene_full_name:
			if debug: print( "Unknown gene {0} not added.".format(gene_full_name))
			return False

		self._net.add_node( 
							tx, 
							median_TPM_N	= None, #Float
							median_PSI_N	= None, #Float
							median_TPM_T	= None, #Float
							median_PSI_T	= None, #Float
							iLoopsFamily 	= None, #String
							gene_id			= geneID#String
						  )

		return True

	def update_node(self, tx, key, value):
		self._update_node(tx, key, value)

	def add_edge(self, tx1, tx2):
		self._net.add_edge( 
							tx1, 							#String
							tx2, 							#String
							iLoops_prediction 	= None, 	#Bool
							RC					= None,		#Float
							experiment 			= None 		#Bool
						  )

	def importTranscriptome(self):

		path = options.Options().qout

		for line in ut.readTable(path + "expression_normal.tsv"):
			txName 			= line[1]
			gene_full_name 	= line[0]
			if not self.add_node(txName, gene_full_name): continue
			if line[3] != "NA":	self.update_node( txName, "median_PSI_N", float(line[3]) )
			if line[5] != "NA": self.update_node( txName, "median_TPM_N", float(line[5]) )

		for line in ut.readTable(path + "expression_tumor.tsv"):
			txName 			= line[1]
			if txName not in self.nodes(): continue
			if line[3] != "NA":	self.update_node( txName, "median_PSI_T", float(line[3]) )
			if line[5] != "NA": self.update_node( txName, "median_PSI_T", float(line[5]) )

		currentLoopFamily 	= ""
		for line in ut.readTable("Data/TCGA/UnifiedFasta_" + options.Options().iLoopsVersion + "_loopFamilies.txt", header=False):
			txName = ""
			if ">" in line[0]:
				currentLoopFamily = line[0][1:]
				txName = line[1]
			else:
				txName = line.pop()
			
			if txName in self.nodes(): 
				self.update_node(txName, "iLoopsFamily", currentLoopFamily)