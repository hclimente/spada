from libs import options
from libs import utils
from methods import method

from collections import Counter
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import pandas as pd
from scipy.stats import fisher_exact

class NeighborhoodAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self.logger.info("Neighborhood Analysis.")

		#C2 Curated sets: canonical pathways
		self.searchEnrichment("canonical_pathways", "c2.cp.v4.0.entrez.gmt","two-sided")
		#C5 GO gene sets: biological process
		self.searchEnrichment("go_biological_process", "c5.bp.v4.0.entrez.gmt","two-sided")
		#C6 Oncogenic signatures
		self.searchEnrichment("oncogenic_signatures", "c6.all.v4.0.entrez.gmt","greater")

		exit()

	def searchEnrichment(self, sTag, sSetFile,H1):

		dGeneSets = {}
		sGeneSetFile = "{0}Data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(sGeneSetFile,header=False):
			sGroup = line[0]
			lGenes = line[2:]
			dGeneSets[sGroup] = {}
			dGeneSets[sGroup]["allGenes"] = lGenes

		for sGeneSet in dGeneSets:

			sInGeneSet_withSwitches		= set()
			sInGeneSet_withoutSwitches	= set()
			sOutGeneSet_withSwitches 	= set()
			sOutGeneSet_withoutSwitches = set()

			#Iterate over the genes that are expressed in the tissue
			for gene,info in [ (x,y) for x,y in self._gene_network.nodes(data=True) if y["ExpressedTranscripts"]]:
				if gene not in dGeneSets[sGeneSet]["allGenes"]:
					if info["isoformSwitches"]:
						sOutGeneSet_withSwitches.add(info["symbol"])
					else:
						sOutGeneSet_withoutSwitches.add(info["symbol"])
				else:
					if info["isoformSwitches"]:
						sInGeneSet_withSwitches.add(info["symbol"])
					else:
						sInGeneSet_withoutSwitches.add(info["symbol"])
		
			iInGeneSet_withSwitches 	= len(sInGeneSet_withSwitches)
			iInGeneSet_withoutSwitches 	= len(sInGeneSet_withoutSwitches)
			iOutGeneSet_withSwitches 	= len(sOutGeneSet_withSwitches)
			iOutGeneSet_withoutSwitches = len(sOutGeneSet_withoutSwitches)

			lContingencyTable = [
									[iInGeneSet_withSwitches,iInGeneSet_withoutSwitches],
									[iOutGeneSet_withSwitches,iOutGeneSet_withoutSwitches]
								]
			fOddsRatio,fPval = fisher_exact(lContingencyTable,H1)

			dGeneSets[sGeneSet]["pval"] = fPval
			dGeneSets[sGeneSet]["oddsRatio"] = fOddsRatio
			dGeneSets[sGeneSet]["switchGenes"] = sInGeneSet_withSwitches

		stats = importr('stats')
		p_adjust = stats.p_adjust(FloatVector([dGeneSets[x]["pval"] for x in dGeneSets]), method='BH')

		with open("{0}neighborhood_analysis/{1}.txt".format(options.Options().qout,sTag),"w") as OUT:
			OUT.write("GeneSet\tpval\tqval\tSwitchingGenes\tOddsRatio\n")
			for sGeneSet,adj_p in zip(dGeneSets,p_adjust):
				OUT.write("{0}\t{1}\t".format(sGeneSet,dGeneSets[sGeneSet]["pval"]))
				OUT.write("{0}\t{1}\t".format(adj_p,",".join(dGeneSets[sGeneSet]["switchGenes"])))
				OUT.write("{0}\n".format(dGeneSets[sGeneSet]["oddsRatio"]))

if __name__ == '__main__':
	pass