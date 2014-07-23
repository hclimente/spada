#!/soft/devel/python-2.7/bin/python

from interface import out_network
from libs import options
from libs import utils
from methods import method
from network import network

from libs.guild.src import combine_scores
from libs.guild.src import create_random_networks_for_netzcore

import pandas as pd

class NetworkAnalysis(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self, __name__, gn_network, tx_network)

	def run(self):
		self.logger.info("GUILD analysis.")
		out_network.Output().getGUILDInput(self._gene_network)

		self.runGUILD( method="NetScore", reps=3, iters=2 )
		self.runGUILD( method="NetZcore", samples=100 )
 		self.runGUILD( method="NetShort" )
 		self.runGUILD( method="fFlow", iters=5, min_seed_score=1.0 )
 		self.runGUILD( method="NetRank" )

		netComboInput = [ "{0}guild_{1}.out".format(options.Options().qout,x) for x in ["NetScore", "NetZcore", "NetShort"] ]
		combine_scores.score_combined(netComboInput, options.Options().qout + "guild_netCombo.out")

		fullComboInput = [ "{0}guild_{1}.out".format(options.Options().qout,x) for x in ["NetScore", "NetZcore", "NetShort", "fFlow", "NetRank"] ]
		combine_scores.score_combined(fullComboInput, options.Options().qout + "guild_fullCombo.out")

	def runGUILD(self, method, reps=None, iters=None, samples=None, min_seed_score=None):
		guild = ["{0}Pipeline/libs/guild/guild".format(options.Options().wd)]
		guild.append( "-n {0}guild_nodes.tsv".format(options.Options().qout) )
		guild.append( "-e {0}guild_edges.tsv".format(options.Options().qout) )
		guild.append( "-s {0}".format(method) )
		guild.append( "-o {0}guild_{1}.out".format(options.Options().qout,method) )

		if reps: 			guild.append( "-r {0}".format(str(reps)) )
		if iters: 			guild.append( "-i {0}".format(str(iters)) )
		if samples: 		guild.append( "-x {0}".format(str(samples)) )
		if min_seed_score:	guild.append( "-t {0}".format(str(min_seed_score)) )

		if method == "NetZcore":
			create_random_networks_for_netzcore.sample_network_preserving_topology(
															options.Options().qout+"guild_edges.tsv", 
															samples, 
															options.Options().qout+"guild_edges.tsv."
														)

		utils.cmd(*guild)

	def readGUILDOutput(self, file):
		##CHECK SPECIALLY THIS
		guildOut = pd.DataFrame.from_csv(options.Options().qout + file, sep="\t", header=None, index_col=None)
		guildOut.columns = ['Gene', 'Score']
		guildOut.Score = guildOut.Score.astype(float)

		for Gene,row in guildOut.iterrows():
			self._gene_network._net.node[Gene]["scoreG"] = row["Score"]

		#guildOut.sort(column="Score", ascending=False, inplace=True)

		#top1 = round(len(guildOut.index) * 0.01)
		#top5 = round(len(guildOut.index) * 0.05)

		# utils.cmd("sort","-g","-r","-k","2","{0}guild_{1}.out".format(options.Options().qout,method), 
		# 		  "|","head","-114",">{0}{1}_top1.txt".format(options.Options().qout,method) )
		# utils.cmd("sort","-g","-r","-k","2","{0}guild_{1}.out".format(options.Options().qout,method), 
		# 		  "|","head","-570",">{0}{1}_top5.txt".format(options.Options().qout,method) )

if __name__ == '__main__':
	#CHECK EVERYTHING
	pass