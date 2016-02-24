#!/soft/devel/python-2.7/bin/python

from interface import out_network
from libs import options
from libs import utils
from methods import method
from network import ucsc_gene_network

from libs.guild.src import combine_scores
from libs.guild.src import create_random_networks_for_netzcore

import pandas as pd
import networkx as nx

class NetworkAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)
		self._gene_subnetworks = {}

	def getGeneSubnetwork(self,top):
		return self._gene_subnetworks[top]

	def run(self, onlyExperimental):
		self.guildAnalysis(onlyExperimental)
		self.functionalNeighbors()

	def functionalNeighbors(self):
		self.logger.info("Switches in driver neighborhoods.")
		candidatesNet = []

		for g in [ x for x,info in self._gene_network.nodes(data=True) if info["Driver"] ]:
			if self._gene_network._net.node[g]["isoformSwitches"]:
				for switch in self._gene_network._net.node[n]["isoformSwitches"]:
					nIso = switch.nTx
					tIso = switch.tTx
					nInfo = self._tx_network._net.node[nIso]
					tInfo = self._tx_network._net.node[tIso]
					if switch.cds_diff and nInfo["iLoopsFamily"] and tInfo["iLoopsFamily"] and nInfo["iLoopsFamily"] != tInfo["iLoopsFamily"]:
						candidatesNet.append(nIso)
						candidatesNet.append(tIso)
			for n in self._gene_network._net.neighbors(g):
				if self._gene_network._net.node[n]["isoformSwitches"]:
					for switch in self._gene_network._net.node[n]["isoformSwitches"]:
						nIso = switch.nTx
						tIso = switch.tTx
						nInfo = self._tx_network._net.node[nIso]
						tInfo = self._tx_network._net.node[tIso]
						if switch.cds_diff and nInfo["iLoopsFamily"] and tInfo["iLoopsFamily"] and nInfo["iLoopsFamily"] != tInfo["iLoopsFamily"]:
							candidatesNet.append(nIso)
							candidatesNet.append(tIso)

		with open(options.Options().qout+"candidatesNetwork.lst","w") as CANDIDATES:
			for iso in candidatesNet:
				CANDIDATES.write(iso+"\t0\tTo analize.\n")

		utils.cmd("scp",options.Options().qout+"candidatesNetwork.lst",
				  "hectorc@gaudi:"+options.Options().gout)

	def guildAnalysis(self,onlyExperimental):
		self.logger.info("GUILD analysis.")
		out_network.getGUILDInput(self._gene_network, onlyExperimental=onlyExperimental)

		self.guildOut = "{0}GUILD_enriched/".format(options.Options().qout)
		if onlyExperimental: self.guildOut = "{0}GUILD_experimental/".format(options.Options().qout)

		self.runGUILD( method="NetScore", reps=3, iters=2 )
		self.runGUILD( method="NetZcore", samples=100 )
		self.runGUILD( method="NetShort" )
		self.runGUILD( method="fFlow", iters=5, min_seed_score=1.0 )
		self.runGUILD( method="NetRank" )

		netComboInput = [ "{0}guild_{1}.out".format(self.guildOut,x) for x in ["NetScore", "NetZcore", "NetShort"] ]
		combine_scores.score_combined( netComboInput, "{0}guild_netCombo.out".format(self.guildOut) )
		fullComboInput = [ "{0}guild_{1}.out".format(self.guildOut,x) for x in ["NetScore", "NetZcore", "NetShort", "fFlow", "NetRank"] ]
		combine_scores.score_combined(fullComboInput, "{0}guild_fullCombo.out".format(self.guildOut) )
		
		self.readGUILDOutput("guild_fullCombo.out")
		self.extractTopSubnetwork(1)
		self.extractTopSubnetwork(5)

	def runGUILD(self, method, reps=None, iters=None, samples=None, min_seed_score=None):

		self.logger.info("Running GUILD, method {0}.".format(method))

		oneLetter = {"NetScore": "s", "NetZcore": "z", "NetShort": "d", "fFlow": "f", "NetRank": "r"}

		guild = ["{0}Pipeline/libs/guild/guild".format(options.Options().wd)]
		guild.append( "-n {0}guild_nodes.tsv".format(self.guildOut) )
		guild.append( "-e {0}guild_edges.tsv".format(self.guildOut) )
		guild.append( "-s {0}".format(oneLetter[method]) )
		guild.append( "-o {0}guild_{1}.out".format(self.guildOut,method) )

		if reps: 			guild.append( "-r {0}".format(str(reps)) )
		if iters: 			guild.append( "-i {0}".format(str(iters)) )
		if samples: 		guild.append( "-x {0}".format(str(samples)) )
		if min_seed_score:	guild.append( "-t {0}".format(str(min_seed_score)) )

		if method == "NetZcore":
			self.logger.debug("Creating 100 random networks for NetZcore.")
			create_random_networks_for_netzcore.sample_network_preserving_topology(
								"{0}guild_edges.tsv".format(self.guildOut), 
								samples, 
								"{0}guild_edges.tsv.".format(self.guildOut)
																				  )
			guild.append("-d {0}guild_edges.tsv.".format(self.guildOut) )

		utils.cmd(*guild)

	def readGUILDOutput(self, guildResults):
		guildOut = pd.DataFrame.from_csv("{0}{1}".format(self.guildOut,guildResults), sep="\t", header=None)
		guildOut.columns = ['Score']

		for Gene,row in guildOut.iterrows():
			self._gene_network.update_node("scoreG",row["Score"],gene_id=str(Gene))

	def extractTopSubnetwork(self, x):
		"""Creates a dot file on GUILD subfolder with the subnetwork top x genes, 
		ranked by GUILD score."""

		topX = 0.01 * x

		topThreshold = round(len(self._gene_network.nodes()) * topX)
		topGenes = []

		sortedNodes = sorted(self._gene_network.nodes(data=True), key=lambda (a, dct): dct['scoreG'], reverse=True)
		counter = 1
		for gene in [ y[0] for y in sortedNodes ]:
			if counter <= topThreshold: 
				topGenes.append(gene)
			else:
				break
			counter +=1

		if options.Options().annotation == "ucsc": 
			self._gene_subnetworks[x] = ucsc_gene_network.UCSCGeneNetwork()
			self._gene_subnetworks[x]._net = nx.subgraph(self._gene_network._net, topGenes)
			
		else:
			self.logger.error("Unrecognized input type {0}.".format(options.Options().annotation))
			exit()

		out_network.outTSV(self._gene_subnetworks[x],"{0}guildTop{1}".format(self.guildOut, x))

if __name__ == '__main__':
	pass