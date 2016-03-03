from libs import options
from libs import utils
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

		self.switchesPerPatient = self.getPatientSwitches()

		self.genesWithAnySwitch = set()
		self.genesWithFunSwitch = set()
		self.allGenes = set()

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network,only_models=True,partialCreation=True,removeNoise=True):
			self.genesWithAnySwitch.add(gene)

			if [ x for x in self._gene_network._net.node[gene]["isoformSwitches"] if self._gene_network.createSwitch(x,self._transcript_network,True).is_functional ]:
				self.genesWithFunSwitch.add(gene)	

		for gene,info in [ (x,y) for x,y in self._gene_network.nodes(data=True) if y["expressedTxsNormal"] or y["expressedTxsTumor"] ]:
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

		geneSets = {}
		geneSetFile = "{0}data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = set(line[2:]) & self.allGenes
			geneSets[geneSet] = {}
			geneSets[geneSet]["allGenes"] = genes

		for s,switches in zip(["all","functional"],[self.genesWithAnySwitch,self.genesWithFunSwitch]):
			for g in geneSets:
				geneSets[g][s] = {}
				genes = geneSets[g]["allGenes"]
				
				geneSets[g][s]["switchIn"] = switches & genes
				geneSets[g][s]["switchOut"] = switches - genes
				geneSets[g][s]["noSwitchIn"] = (self.allGenes - switches) & genes
				geneSets[g][s]["noSwitchOut"] = (self.allGenes - switches) - genes
				geneSets[g][s]["switchGenes"] = sorted([ self._gene_network._net.node[x]["symbol"] for x in geneSets[g][s]["switchIn"] ])

				geneSets[g][s]["sg"]   = len(geneSets[g][s]["switchIn"])
				geneSets[g][s]["nsg"]  = len(geneSets[g][s]["noSwitchIn"])
				geneSets[g][s]["sng"]  = len(geneSets[g][s]["switchOut"])
				geneSets[g][s]["nsng"] = len(geneSets[g][s]["noSwitchOut"])

				lContingencyTable = [[geneSets[g][s]["sg"],geneSets[g][s]["nsg"]],
									[geneSets[g][s]["sng"],geneSets[g][s]["nsng"]]]
				OR,pval = fisher_exact(lContingencyTable,H1)

				affection = []
				for p in self.switchesPerPatient:
					affection.append(len(genes & self.switchesPerPatient[p]))

				geneSets[g][s]["affection"] = (float(sum(affection))/len(affection))/len(genes)
				geneSets[g][s]["pval"] = pval
				geneSets[g][s]["oddsRatio"] = OR
		
			p_adjust = multipletests([geneSets[x][s]["pval"] for x in geneSets],alpha=0.05,method='fdr_bh')
			p_adjust = p_adjust[1].tolist()

			with open("{}neighborhood_analysis/{}_{}{}.txt".format(options.Options().qout,sGenesetTag,s,options.Options().filetag),"w") as OUT:
				OUT.write("GeneSet\tCancer\tpval\tfdr5\tNormalizedAverageAffection\t")
				OUT.write("SwitchingGenes\tOR\tsg\tsng\tnsg\tnsng\n")
				for g,adj_p in zip(geneSets,p_adjust):

					OUT.write("{}\t{}\t{}\t".format(g,options.Options().tag,geneSets[g][s]["pval"]))
					OUT.write("{}\t{}\t{}\t".format(adj_p,geneSets[g][s]["affection"], ",".join(geneSets[g][s]["switchGenes"]) ))
					OUT.write("{}\t{}\t{}\t".format(geneSets[g][s]["oddsRatio"],geneSets[g][s]["sg"],geneSets[g][s]["sng"]))
					OUT.write("{}\t{}\n".format(geneSets[g][s]["nsg"],geneSets[g][s]["nsng"]))

	def getPatientSwitches(self):

		patients = []
		[ patients.extend(z["patients"]) for x,y in self._gene_network.iterate_genes_byPatientNumber() for z in y["isoformSwitches"] ]
		patients = list(set(patients))

		switchesPerPatient = dict([ [x,set()] for x in patients ])

		for g in self._gene_network.nodes():
			info = self._gene_network._net.node[g]

			for s in info["isoformSwitches"]:
				for p in s["patients"]:
					switchesPerPatient[p].add(g)

		return switchesPerPatient

	def findAffectedPathways(self,sTag,sSetFile):
		affectedPathway = {}
		geneSetFile = "{0}data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = line[2:]
			geneSets[geneSet] = {}
			geneSets[geneSet]["allGenes"] = genes

		for geneSet in geneSets:
			affectedPathway[geneSet] = [0]*len(options.Options().replicates)

			for gene,info,switchDict,switch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network):
				if gene not in geneSets[geneSet]: continue
				
				switchSpread = [ 1 if x in switch.patients else 0 for x in options.Options().replicates ]
				affectedPathway[geneSet] = [ x+y for x,y in zip(affectedPathway[geneSet],switchSpread)]

		with open("{0}neighborhood_analysis/{1}_patientAffectation{2}.txt".format(options.Options().qout,sTag,options.Options().filetag),"w") as OUT:
			OUT.write("GeneSet\t{0}\n".format('\t'.join(map(str,options.Options().replicates) ) ) )

			for geneSet in affectedPathway:
				OUT.write("{0}\t{1}\n".format(geneSet,'\t'.join(map(str,affectedPathway[geneSet]))))

if __name__ == '__main__':
	pass