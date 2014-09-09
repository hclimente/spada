from libs import options
from libs import utils
from methods import method

from collections import Counter
import fisher
import pandas as pd

class NeighborhoodAnalysis(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self.logger.info("Neighborhood Analysis.")

		#C2 Curated sets: canonical pathways
		self.analyzeCategories("canonical_pathways", "c2.cp.v4.0.entrez.gmt")
		#C5 GO gene sets: biological process
		self.analyzeCategories("go_biological_process", "c5.bp.v4.0.entrez.gmt")
		#C6 Oncogenic signatures
		self.analyzeCategories("oncogenic_signatures", "c6.all.v4.0.entrez.gmt")

		exit()

	def analyzeCategories(self, sTag, sSetFile):

		dGeneSets = {}

		for line in utils.readTable("{0}Data/Databases/{1}".format(options.Options().wd,sSetFile),header=False):
			sGroup = line[0]
			lGenes = line[2:]
			dGeneSets[sGroup] = lGenes

		with open("{0}neighborhood_analysis/{1}.txt".format(options.Options().qout,sTag),"w") as OUT:
			for sGeneSet in dGeneSets:

				sInGeneSet_withSwitches		= set()
				sInGeneSet_withoutSwitches	= set()
				sOutGeneSet_withSwitches 	= set()
				sOutGeneSet_withoutSwitches = set()

				for gene,info in self._gene_network.nodes(data=True):
					if gene not in dGeneSets[sGeneSet]:
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

				fPValue = fisher.pvalue(iInGeneSet_withSwitches,iInGeneSet_withoutSwitches,iOutGeneSet_withSwitches,iOutGeneSet_withoutSwitches)

				OUT.write("{0}\t{1}\t{2}\n".format(sGeneSet,fPValue.two_tail,",".join(sInGeneSet_withSwitches)))

if __name__ == '__main__':
	pass