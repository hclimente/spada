#!/soft/devel/python-2.7/bin/python

from libs.utils import *
from network import network
from methods import method

class NetworkAnalysis(method.Method):
	def __init__(self, path, samples, candidates):
		method.Method.__init__(self, path=path, samples=samples, candidates=candidates)
		self._predicted_interactome		= network.Network()
		self._experimental_interactome	= network.Network()
		self._full_interactome 			= network.Network()

	def predInteractome(self):	return self._predicted_interactome
	def expeInteractome(self):	return self._experimental_interactome
	def fullInteractome(self):	return self._full_interactome

	def importPredictedInteractome(self, drivers):
		self._predicted_interactome.importCandidates(self._path, self._samples)
		self._predicted_interactome.importDrivers(drivers)
		for gn,info in self._candidates.iteritems():
			self._predicted_interactome.importiLoopsInteractions(self._path, gn, info["nTx"], info["tTx"])
		
	def importExperimentalInteractome(self, drivers):
		self._experimental_interactome.importCandidates(self._path, self._samples)
		self._experimental_interactome.importDrivers(drivers)
		self._experimental_interactome.importKnownInteractions()

	def predictedInteractome(self): return self._predicted_interactome
	def experimentalInteractome(self): return self._experimental_interactome 

	def joinNetworks(self):
		for node in [ x for x in self.expeInteractome().nodes() if x not in self.predInteractome().nodes() ]:
			Score 		= self._experimental_interactome.n().node[node]["score"]
			geneSymbol	= self._experimental_interactome.n().node[node]["symbol"]
			Switch		= self._experimental_interactome.n().node[node]["switch"]
			Driver		= self._experimental_interactome.n().node[node]["driver"]

			self._full_interactome.add_node(score=Score, gene_id=node, gene_symbol=geneSymbol, switch=Switch, driver=Driver)

		for node in [ x for x in self.predInteractome().nodes() if x not in self.expeInteractome().nodes() ]:
			Score 		= self._predicted_interactome.n().node[node]["score"]
			geneSymbol	= self._predicted_interactome.n().node[node]["symbol"]
			Switch		= self._predicted_interactome.n().node[node]["switch"]
			Driver		= self._predicted_interactome.n().node[node]["driver"]

			self._full_interactome.add_node(score=Score, gene_id=node, gene_symbol=geneSymbol, switch=Switch, driver=Driver)

		for node in [ x for x in self.predInteractome().nodes() if x in self.expeInteractome().nodes() ]:
			Score_exp 		= self._experimental_interactome.n().node[node]["score"]
			geneSymbol_exp	= self._experimental_interactome.n().node[node]["symbol"]
			Switch_exp		= self._experimental_interactome.n().node[node]["switch"]
			Driver_exp		= self._experimental_interactome.n().node[node]["driver"]
			Score_pre 		= self._predicted_interactome.n().node[node]["score"]
			geneSymbol_pre	= self._predicted_interactome.n().node[node]["symbol"]
			Switch_pre		= self._predicted_interactome.n().node[node]["switch"]
			Driver_pre		= self._predicted_interactome.n().node[node]["driver"]

			Score 		= Score_pre + Score_exp
			geneSymbol	= geneSymbol_pre
			if not geneSymbol_pre == geneSymbol_exp:
				print( node + ": gene symbols do not match: " + geneSymbol_exp + " and " + geneSymbol_pre)
			
			Switch		= Switch_pre
			if not Switch_pre == Switch_exp:
				print( node + ": switch information do not match." )
			Driver		= Driver_pre
			if not Driver_pre == Driver_exp:
				print( node + ": driver information do not match." )

			self._full_interactome.add_node(score=Score, gene_id=node, gene_symbol=geneSymbol, switch=Switch, driver=Driver)

		for node1, node2 in [ x for x in self.predInteractome().edges() if x not in self.expeInteractome().edges() ]:
			Weight 	= self._predicted_interactome[node1][node2]["weight"]
			Methods = self._predicted_interactome[node1][node2]["methods"]
			Sources = self._predicted_interactome[node1][node2]["sources"]

			self._predicted_interactome.add_edge(node1, node2, weight=Weight, methods=Methods, sources=Sources)

		for node1, node2 in [ x for x in self.expeInteractome().edges() if x not in self.predInteractome().edges() ]:
			Weight 	= self._experimental_interactome[node1][node2]["weight"]
			Methods = self._experimental_interactome[node1][node2]["methods"]
			Sources = self._experimental_interactome[node1][node2]["sources"]

			self._experimental_interactome.add_edge(node1, node2, weight=Weight, methods=Methods, sources=Sources)

		for node1, node2 in [ x for x in self.predInteractome().edges() if x in self.expeInteractome().edges() ]:
			Weight_exp 	= self._experimental_interactome[node1][node2]["weight"]
			Methods_exp = self._experimental_interactome[node1][node2]["methods"]
			Sources_exp = self._experimental_interactome[node1][node2]["sources"]
			Weight_pre 	= self._predicted_interactome[node1][node2]["weight"]
			Methods_pre = self._predicted_interactome[node1][node2]["methods"]
			Sources_pre = self._predicted_interactome[node1][node2]["sources"]

			Weight 	= Weight_pre + Weight_exp
			Methods = Methods_pre + Methods_exp
			Sources = Sources_pre + Sources_exp

			self._experimental_interactome.add_edge(node1, node2, weight=Weight, methods=Methods, sources=Sources)

if __name__ == '__main__':
	#path = sys.argv[1]
	path = "/home/hector/SmartAS/Results/TCGA/luad_mE-1.0/"
	samples = 57
	candidates = { "ARL1" : {"nTx": "uc001tib.2", "tTx": "uc001tic.2"} }
	driversFile = "/home/hector/Desktop/Baldo_drivers3.lst"
	
	analysis = NetworkAnalysis(path, samples, candidates)
	analysis.importExperimentalInteractome(driversFile)
	analysis.importPredictedInteractome(driversFile)
	analysis.joinNetworks()