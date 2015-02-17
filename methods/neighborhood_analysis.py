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
		#C7 Immunologic signatures
		self.searchEnrichment("immunologic_signatures", "c7.all.v4.0.entrez.gmt","greater")

		#Look for enrichment in genes known to be involved in AS in cancer
		self.searchEnrichment("asAffectedGenes_Angiogenesis","asAffected_Angiogenesis_list.csv","greater")
		self.searchEnrichment("asAffectedGenes_Genes_with_cancer_AS_events","asAffected_Genes_with_cancer_AS_events_list.csv","greater")
		self.searchEnrichment("asAffectedGenes_Apoptosis","asAffected_Apoptosis_list.csv","greater")
		self.searchEnrichment("asAffectedGenes_InvasionAndMetastasis","asAffected_Invasion_&_Metastasis_list.csv","greater")
		self.searchEnrichment("asAffectedGenes_Cancer_Therapy","asAffected_Cancer_Therapy_list.csv","greater")
		self.searchEnrichment("asAffectedGenes_Proliferation","asAffected_Proliferation_list.csv","greater")
		self.searchEnrichment("asAffectedGenes_DNA_damage","asAffected_DNA_damage_list.csv","greater")

	def searchEnrichment(self,sTag,sSetFile,H1):

		geneSets = {}
		geneSetFile = "{0}Data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = line[2:]
			geneSets[geneSet] = {}
			geneSets[geneSet]["allGenes"] = genes
			geneSets[geneSet]["switchIn"] = set()
			geneSets[geneSet]["switchOut"] = set()
			geneSets[geneSet]["noswitchIn"] = set()
			geneSets[geneSet]["noswitchOut"] = set()

		# iterate genes with isoform switches
		for gene,info,switchDict,switch in self._gene_network.iterate_relevantSwitches_ScoreWise(self._transcript_network,only_models=options.Options().onlyModels,partialCreation=True):
			#Iterate over the genes that are expressed in the tissue
			for geneSet in geneSets:
				if gene not in geneSets[geneSet]["allGenes"]:
					geneSets[geneSet]["switchOut"].add(gene)
				else:
					self._gene_network.update_node("neighborhoods",geneSet,gene_id=gene,secondKey=sTag)
					geneSets[geneSet]["switchIn"].add(gene)

		# iterate all genes and disregard previously iterated
		for gene,info in [ (x,y) for x,y in self._gene_network.nodes(data=True) if y["ExpressedTranscripts"] ]:
			for geneSet in geneSets:
				if gene in geneSets[geneSet]["switchOut"] or gene in geneSets[geneSet]["switchIn"]:
					continue
				
				if gene not in geneSets[geneSet]["allGenes"]:
					geneSets[geneSet]["noswitchOut"].add(gene)
				else:
					geneSets[geneSet]["noswitchIn"].add(gene)
		
		for geneSet in geneSets:
			iInGeneSet_withSwitches 	= len(geneSets[geneSet]["switchIn"])
			iInGeneSet_withoutSwitches 	= len(geneSets[geneSet]["noswitchIn"])
			iOutGeneSet_withSwitches 	= len(geneSets[geneSet]["switchOut"])
			iOutGeneSet_withoutSwitches = len(geneSets[geneSet]["noswitchOut"])

			lContingencyTable = [[iInGeneSet_withSwitches,iInGeneSet_withoutSwitches],
								[iOutGeneSet_withSwitches,iOutGeneSet_withoutSwitches]]
			fOddsRatio,fPval = fisher_exact(lContingencyTable,H1)

			geneSets[geneSet]["pval"] = fPval
			geneSets[geneSet]["oddsRatio"] = fOddsRatio

		stats = importr('stats')
		p_adjust = stats.p_adjust(FloatVector([geneSets[x]["pval"] for x in geneSets]), method='BH')

		with open("{0}neighborhood_analysis/{1}{2}.txt".format(options.Options().qout,sTag,options.Options().filetag),"w") as OUT:
			OUT.write("GeneSet\tpval\tqval\tSwitchingGenes\tOddsRatio\n")
			for geneSet,adj_p in zip(geneSets,p_adjust):

				switchGenes = ",".join([ self._gene_network._net.node[x]["symbol"] for x in geneSets[geneSet]["switchIn"] ])

				OUT.write("{0}\t{1}\t".format(geneSet,geneSets[geneSet]["pval"]))
				OUT.write("{0}\t{1}\t".format(adj_p,switchGenes))
				OUT.write("{0}\n".format(geneSets[geneSet]["oddsRatio"]))

	def findAffectedPathways(self,sTag,sSetFile):
		affectedPathway = {}
		geneSetFile = "{0}Data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = line[2:]
			geneSets[geneSet] = {}
			geneSets[geneSet]["allGenes"] = genes

		for geneSet in geneSets:
			affectedPathway[geneSet] = [0]*len(options.Options().replicates)

			for gene,info,switchDict,switch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network):
				if gene not in geneSets[geneSet]: continue
				
				switchSpread = [ 1 if x in switch.patients else 0 for x in options.Options().replicates ]
				affectedPathway[geneSet] = [ x+y for x,y in zip(affectedPathway[geneSet],switchSpread)]

		with open("{0}neighborhood_analysis/{1}_patientAffectation{2}.txt".format(options.Options().qout,sTag,options.Options().filetag),"w") as OUT:
			OUT.write("GeneSet\t{0}\n".format('\t'.join(map(str,options.Options().replicates) ) ) )

			for geneSet in affectedPathway:
				OUT.write("{0}\t{1}\n".format(geneSet,'\t'.join(map(str,affectedPathway[geneSet]))))


if __name__ == '__main__':
	pass