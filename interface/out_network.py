from libs import options

import math
import networkx

class OutNetwork:
	def __init__(self):
		pass

	def interactorsInfo(self, gn_net, tx_net, gene, nIso, tIso):
		geneSym = gn_net._net.node[gene]["symbol"]
		with open("{0}InteraXv2_{1}_{2}_{3}.tsv".format(options.Options().qout, geneSym, nIso, tIso), "w" ) as INTERAX:
			INTERAX.write("#")
			INTERAX.write("\t{0}\t{1}".format( nIso, tx_net._net.node[nIso]["median_PSI_N"] ))
			INTERAX.write("\t{0}\t{1}".format( tIso, tx_net._net.node[tIso]["median_PSI_T"] ))
			INTERAX.write("\t{0}".format( gene ))
			INTERAX.write("\t{0}".format( geneSym ))
			INTERAX.write("\t{0}".format( math.log10 ( tx_net._net.node[nIso]["median_TPM_N"] + 0.0001 ) ))
			INTERAX.write("\t{0}".format( math.log10 ( tx_net._net.node[tIso]["median_TPM_N"] + 0.0001 ) ))
			INTERAX.write("\n")

			for partnerIso in set( tx_net._net.neighbors(nIso) ) | set( tx_net._net.neighbors(tIso) ):
				partnerGene = tx_net._net.node[partnerIso]["gene_id"]
				annotation = "Driver" if gn_net._net.node[partnerGene]["Driver"] else ""
				RCN = "" if partnerIso not in tx_net._net.neighbors(nIso) else tx_net._net.edge[nIso][partnerIso]["RC"]
				RCT = "" if partnerIso not in tx_net._net.neighbors(tIso) else tx_net._net.edge[tIso][partnerIso]["RC"]
				dRC = int(RCN - RCT) if RCN and RCT else ""

				INTERAX.write("{0}\t".format( partnerIso ))
				INTERAX.write("{0}\t".format( partnerGene ))
				INTERAX.write("{0}\t".format( gn_net._net.node[ partnerGene ]["symbol"] ))
				INTERAX.write("{0}\t".format( RCN ))
				INTERAX.write("{0}\t".format( RCT ))
				INTERAX.write("{0}\t".format( dRC ))
				INTERAX.write("{0}\t".format( annotation ))
				INTERAX.write("{0}\t".format( math.log10 ( tx_net._net.node[partnerIso]["median_TPM_N"] + 0.0001 ) ))
				INTERAX.write("{0}\t".format( tx_net._net.node[partnerIso]["median_PSI_N"] ))
				INTERAX.write("{0}\t".format( math.log10 ( tx_net._net.node[partnerIso]["median_TPM_T"] + 0.0001 ) ))
				INTERAX.write("{0}\t".format( tx_net._net.node[partnerIso]["median_PSI_T"] ))
				INTERAX.write("\n")

	def getGUILDInput(self):
		with open(options.Options().qout + "guild_nodes.tsv" ) as GUILD_NODES:
			for node, info in self.nodes(data=True):
				GUILD_NODES.write("{0}\t{1}\n".format(node, info["score"]))

		with open(options.Options().qout + "guild_edges.tsv" ) as GUILD_EDGES:
			for node1, node2, info in self.edges(data=True):
				GUILD_EDGES.write("{0}\t{1}\t{2}\n".format(node1, info["score"], node2))

	def outputDot(self, network, name):
		networkx.write_dot(network, options.Options().qout + name)