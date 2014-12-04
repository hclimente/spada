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

		for geneSet in geneSets:

			inGeneSet_withSwitches		= set()
			inGeneSet_withoutSwitches	= set()
			outGeneSet_withSwitches 	= set()
			outGeneSet_withoutSwitches 	= set()

			#Iterate over the genes that are expressed in the tissue
			for gene,info in [ (x,y) for x,y in self._gene_network.nodes(data=True) if y["ExpressedTranscripts"]]:
				if gene not in geneSets[geneSet]["allGenes"]:
					if [ x for x in info["isoformSwitches"] if x.is_relevant]:
						outGeneSet_withSwitches.add(gene)
					else:
						outGeneSet_withoutSwitches.add(gene)
				else:
					if [ x for x in info["isoformSwitches"] if x.is_relevant]:
						self._gene_network.update_node("neighborhoods",geneSet,gene_id=gene,secondKey=sTag)
						inGeneSet_withSwitches.add(gene)
					else:
						inGeneSet_withoutSwitches.add(gene)
		
			iInGeneSet_withSwitches 	= len(inGeneSet_withSwitches)
			iInGeneSet_withoutSwitches 	= len(inGeneSet_withoutSwitches)
			iOutGeneSet_withSwitches 	= len(outGeneSet_withSwitches)
			iOutGeneSet_withoutSwitches = len(outGeneSet_withoutSwitches)

			lContingencyTable = [
									[iInGeneSet_withSwitches,iInGeneSet_withoutSwitches],
									[iOutGeneSet_withSwitches,iOutGeneSet_withoutSwitches]
								]
			fOddsRatio,fPval = fisher_exact(lContingencyTable,H1)

			geneSets[geneSet]["pval"] = fPval
			geneSets[geneSet]["oddsRatio"] = fOddsRatio
			geneSets[geneSet]["switchGenes"] = inGeneSet_withSwitches

		stats = importr('stats')
		p_adjust = stats.p_adjust(FloatVector([geneSets[x]["pval"] for x in geneSets]), method='BH')

		with open("{0}neighborhood_analysis/{1}.txt".format(options.Options().qout,sTag),"w") as OUT:
			OUT.write("GeneSet\tpval\tqval\tSwitchingGenes\tOddsRatio\n")
			for geneSet,adj_p in zip(geneSets,p_adjust):

				switchGenes = ",".join([ self._gene_network._net.node[x]["symbol"] for x in geneSets[geneSet]["switchGenes"] ])

				OUT.write("{0}\t{1}\t".format(geneSet,geneSets[geneSet]["pval"]))
				OUT.write("{0}\t{1}\t".format(adj_p,switchGenes))
				OUT.write("{0}\n".format(geneSets[geneSet]["oddsRatio"]))

				# fill switch property in the positives
				if adj_p <= 0.05:
					for gene in geneSets[geneSet]["switchGenes"]:
						for switch in self._gene_network._net.node[x]["isoformSwitches"]:
							if switch._neighborhood_change is None:
								switch._neighborhood_change = set()
							switch._neighborhood_change.add((sTag,geneSet))

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

			for gene,info,switch in utils.iterate_switches_ScoreWise(self._gene_network):
				if gene not in geneSets[geneSet]: continue
				
				switchSpread = [ 1 if x in switch.patients else 0 for x in options.Options().replicates ]
				affectedPathway[geneSet] = [ x+y for x,y in zip(affectedPathway[geneSet],switchSpread)]

		with open("{0}neighborhood_analysis/{1}_patientAffectation.txt".format(options.Options().qout,sTag),"w") as OUT:
			OUT.write("GeneSet\t{0}\n".format('\t'.join(map(str,options.Options().replicates) ) ) )

			for geneSet in affectedPathway:
				OUT.write("{0}\t{1}\n".format(geneSet,'\t'.join(map(str,affectedPathway[geneSet]))))


if __name__ == '__main__':
	pass