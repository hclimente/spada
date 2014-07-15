from libs import options

class OutNetwork:
	def __init__(self):
		pass

	def interactorsInfo(self, gn_net, tx_net, gene, nIso, tIso):
		with open(options.Options().qout + "InteraX_{0}_{1}_{2}.tsv".format(gene, nIso, tIso) ) as INTERAX:
			INTERAX.write("#")
			INTERAX.write("\t{0}\t{1}".format( nIso, tx_net._net.node[nIso]["median_PSI_N"] ))
			INTERAX.write("\t{0}\t{1}".format( tIso, tx_net._net.node[tIso]["median_PSI_T"] ))
			INTERAX.write("\t{0}".format( tx_net._net.node[nIso]["gene_id"] ))
			INTERAX.write("\t{0}".format( tx_net._net.node[nIso]["symbol"] ))
			INTERAX.write("\t{0}".format( log ( tx_net._net.node[nIso]["median_TPM_N"] + 0.0001 ) ))
			INTERAX.write("\t{0}".format( log ( tx_net._net.node[tIso]["median_TPM_N"] + 0.0001 ) ))
			INTERAX.write("\n")

			for partner in set( tx_net._net.neighbors(nIso) ) | set( tx_net._net.neighbors(tIso) ):
				gn = tx_net._net.node[partner]["gene_id"]
				
				RCN = tx_net._net.edge[nIso][partner]["RC"] if tx_net._net.edge[nIso][partner]["RC"] else ""
				RCT = tx_net._net.edge[tIso][partner]["RC"] if tx_net._net.edge[tIso][partner]["RC"] else ""
				dRC = int(RCN) - int(RCT) if RCN and RCT else ""
				annotation = "Driver" if gn_net._net.node[gn]["Driver"] else ""

				INTERAX.write("{0}\t".format( partner ))
				INTERAX.write("{0}\t".format( tx_net._net.node[partner]["gene_id"] ))
				INTERAX.write("{0}\t".format( gn_net._net.node[ gn ]["symbol"] ))
				INTERAX.write("{0}\t".format( RCN ))
				INTERAX.write("{0}\t".format( RCT ))
				INTERAX.write("{0}\t".format( dRC ))
				INTERAX.write("{0}\t".format( annotation ))
				INTERAX.write("{0}\t".format( log ( tx_net._net.node[partner]["median_TPM_N"] + 0.0001 ) ))
				INTERAX.write("{0}\t".format( tx_net._net.node[partner]["median_PSI_N"] ))
				INTERAX.write("{0}\t".format( log ( tx_net._net.node[partner]["median_TPM_T"] + 0.0001 ) ))
				INTERAX.write("{0}\t".format( tx_net._net.node[partner]["median_PSI_T"] ))
			