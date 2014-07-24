#!/soft/devel/python-2.7/bin/python

from interface import out_network
from libs import options
from libs import utils
from methods import method
from network import network

from libs.guild.src import combine_scores
from libs.guild.src import create_random_networks_for_netzcore

import pandas as pd
import networkx

class NetworkAnalysis(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)

	def run(self, onlyExperimental):
		self.logger.info("Network analysis.")
		out_network.getGUILDInput(self._gene_network, onlyExperimental=onlyExperimental)

		self.runGUILD( method="NetScore", reps=3, iters=2 )
		self.runGUILD( method="NetZcore", samples=100 )
		self.runGUILD( method="NetShort" )
		self.runGUILD( method="fFlow", iters=5, min_seed_score=1.0 )
		self.runGUILD( method="NetRank" )

		netComboInput = [ "{0}GUILD/guild_{1}.out".format(options.Options().qout,x) for x in ["NetScore", "NetZcore", "NetShort"] ]
		combine_scores.score_combined(netComboInput, options.Options().qout + "GUILD/guild_netCombo.out")

		fullComboInput = [ "{0}GUILD/guild_{1}.out".format(options.Options().qout,x) for x in ["NetScore", "NetZcore", "NetShort", "fFlow", "NetRank"] ]
		combine_scores.score_combined(fullComboInput, options.Options().qout + "GUILD/guild_fullCombo.out")

		self.readGUILDOutput("GUILD/guild_fullCombo.out")
		self.extractTopSubnetwork(1)
		self.extractTopSubnetwork(5)

	def runGUILD(self, method, reps=None, iters=None, samples=None, min_seed_score=None):

		self.logger.info("Running GUILD, method {0}.".format(method))

		oneLetter = {"NetScore": "s", "NetZcore": "z", "NetShort": "d", "fFlow": "f", "NetRank": "r"}

		guild = ["{0}Pipeline/libs/guild/guild".format(options.Options().wd)]
		guild.append( "-n {0}GUILD/guild_nodes.tsv".format(options.Options().qout) )
		guild.append( "-e {0}GUILD/guild_edges.tsv".format(options.Options().qout) )
		guild.append( "-s {0}".format(oneLetter[method]) )
		guild.append( "-o {0}GUILD/guild_{1}.out".format(options.Options().qout,method) )

		if reps: 			guild.append( "-r {0}".format(str(reps)) )
		if iters: 			guild.append( "-i {0}".format(str(iters)) )
		if samples: 		guild.append( "-x {0}".format(str(samples)) )
		if min_seed_score:	guild.append( "-t {0}".format(str(min_seed_score)) )

		if method == "NetZcore":
			self.logger.debug("Creating 100 random networks for NetZcore.")
			create_random_networks_for_netzcore.sample_network_preserving_topology(
															options.Options().qout + "GUILD/guild_edges.tsv", 
															samples, 
															options.Options().qout + "GUILD/guild_edges.tsv."
														)
			guild.append("-d {0}GUILD/guild_edges.tsv.".format(options.Options().qout))

		utils.cmd(*guild)

	def readGUILDOutput(self, guildResults):
		guildOut = pd.DataFrame.from_csv(options.Options().qout + guildResults, sep="\t", header=None)
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

		topNetwork = networkx.subgraph(self._gene_network._net, topGenes)

		out_network.outputDot(topNetwork, "GUILD/guildTop{0}.dot".format(x) )

if __name__ == '__main__':

	pass