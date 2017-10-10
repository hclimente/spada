from spada import options
from spada import utils
from methods import method

from collections import Counter
from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd
from scipy.stats import fisher_exact

import pdb

class NeighborhoodAnalysis(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)
		self.logger.info("Preparing neighborhood analysis.")

		utils.cmd("mkdir","{}neighborhood_analysis".format(options.Options().qout))

		self.genesWithPPISwitch = set()
		self.allGenes = set()

		for line in utils.readTable("{}candidateList_driverEvidence.tsv".format(options.Options().qout)):
			ppi = line[7]

			if ppi == "1":
				gene = line[1]

				self.genesWithPPISwitch.add(gene)	

		for gene,info in self._gene_network.iterate_genes_byPatientNumber(onlySplicedGenes=True,onlyExpressedGenes=True,alwaysSwitchedGenes=True):
			self.allGenes.add(gene)

	def run(self):

		self.logger.info("Running neighborhood analysis.")

		#C2 Curated sets: canonical pathways
		self.searchEnrichment("canonical_pathways","c2.cp.v4.0.entrez.gmt","two-sided")
		#C5 GO gene sets: biological process
		self.searchEnrichment("go_biological_process","c5.bp.v4.0.entrez.gmt","two-sided")
		#C6 Oncogenic signatures
		self.searchEnrichment("oncogenic_signatures","c6.all.v4.0.entrez.gmt","greater")
		#Cancer hallmarks
		self.searchEnrichment("hallmarks","h.all.v5.0.entrez.gmt","greater")

	def searchEnrichment(self,sGenesetTag,sSetFile,H1):

		sets = {}
		setfile = "{0}data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(setfile,header=False):
			gs = line[0]
			genes = set(line[2:]) & self.allGenes
			sets[gs] = genes

		gsTest = {}
		for g in sets:
			gsTest[g] = {}
			genes = sets[g]
			
			sg   = len(self.genesWithPPISwitch & genes)
			nsg  = len((self.allGenes - self.genesWithPPISwitch) & genes)
			sng  = len(self.genesWithPPISwitch - genes)
			nsng = len((self.allGenes - self.genesWithPPISwitch) - genes)

			OR,pval = fisher_exact([[sg,nsg], [sng,nsng]],H1)

			gsTest[g]["pval"] = pval
			gsTest[g]["oddsRatio"] = OR
			gsTest[g]["switchGenes"] = sorted([ self._gene_network._net.node[x]["symbol"] for x in self.genesWithPPISwitch & genes ])
	
		p_adjust = multipletests([gsTest[x]["pval"] for x in gsTest],alpha=0.05,method='fdr_bh')
		p_adjust = p_adjust[1].tolist()

		with open("{}neighborhood_analysis/{}.ppi.txt".format(options.Options().qout,sGenesetTag),"w") as OUT:
			OUT.write("Tumor\tGeneSet\tpval\tfdr5\tOR\t")
			OUT.write("SwitchedGenes\n")
			for g,adj_p in zip(gsTest,p_adjust):

				OUT.write("{}\t{}\t{}\t".format(options.Options().tag,g,gsTest[g]["pval"]))
				OUT.write("{}\t{}\t".format(adj_p, gsTest[g]["oddsRatio"] ))
				OUT.write("{}\n".format(",".join(gsTest[g]["switchGenes"]) ))

if __name__ == '__main__':
	pass