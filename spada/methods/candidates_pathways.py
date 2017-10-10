from spada import options
from spada import utils
from methods import method

from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import fisher_exact

class CandidatesPathways(method.Method):
	def __init__(self, gn_network, tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)
		self.logger.info("Preparing neighborhood analysis.")

		utils.cmd("mkdir","{}neighborhood_analysis".format(options.Options().qout))

		self.switchesPerPatient = self.getPatientSwitches()

		self.genesWithCandidateSwitch = set()
		self.allGenes = set()

		# read recurrence
		for line in utils.readTable("{}candidateList_recurrence.tsv".format(options.Options().qout)):
			if line[5]!="NA" and float(line[5]) < 0.05:
				self.genesWithCandidateSwitch.add(line[0])

		# read co-ocurrence
		for line in utils.readTable("{}candidateList_mutationCoocurrence.tsv".format(options.Options().qout)):
			if line[8]!="NA" and float(line[8]) < 0.05:
				self.genesWithCandidateSwitch.add(line[0])

		# read mutual exclusion
		for line in utils.readTable("{}candidateList_mutationME.tsv".format(options.Options().qout)):
			if line[4]!="NA" and bool(int(line[4])):
				self.genesWithCandidateSwitch.add(line[0])

		# read switches affecting mutated features
		for line in utils.readTable("{}candidateList_mutatedFeatures.tsv".format(options.Options().qout)):
			if line[4]!="NA" and bool(int(line[4])):
				self.genesWithCandidateSwitch.add(line[0])

		for gene,info in [ (x,y) for x,y in self._gene_network.nodes(data=True) if y["expressedTxsNormal"] or y["expressedTxsTumor"] ]:
			self.allGenes.add(gene)

	def run(self):

		self.searchEnrichment("canonical_pathways","c2.cp.v4.0.entrez.gmt","greater")

	def searchEnrichment(self,sGenesetTag,sSetFile,H1):

		geneSets = {}
		geneSetFile = "{0}data/Databases/{1}".format(options.Options().wd,sSetFile)

		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = set(line[2:]) & self.allGenes
			geneSets[geneSet] = {}
			geneSets[geneSet]["allGenes"] = genes

		for g in geneSets:
			genes = geneSets[g]["allGenes"]
			
			geneSets[g]["switchIn"] = self.genesWithCandidateSwitch & genes
			geneSets[g]["switchOut"] = self.genesWithCandidateSwitch - genes
			geneSets[g]["noSwitchIn"] = (self.allGenes - self.genesWithCandidateSwitch) & genes
			geneSets[g]["noSwitchOut"] = (self.allGenes - self.genesWithCandidateSwitch) - genes
			geneSets[g]["switchGenes"] = sorted([ self._gene_network._net.node[x]["symbol"] for x in geneSets[g]["switchIn"] ])

			geneSets[g]["sg"]   = len(geneSets[g]["switchIn"])
			geneSets[g]["nsg"]  = len(geneSets[g]["noSwitchIn"])
			geneSets[g]["sng"]  = len(geneSets[g]["switchOut"])
			geneSets[g]["nsng"] = len(geneSets[g]["noSwitchOut"])

			lContingencyTable = [[geneSets[g]["sg"],geneSets[g]["nsg"]],
								[geneSets[g]["sng"],geneSets[g]["nsng"]]]
			OR,pval = fisher_exact(lContingencyTable,H1)

			affection = []
			for p in self.switchesPerPatient:
				affection.append(len(genes & self.switchesPerPatient[p]))

			geneSets[g]["affection"] = (float(sum(affection))/len(affection))/len(genes)
			geneSets[g]["pval"] = pval
			geneSets[g]["oddsRatio"] = OR
	
		p_adjust = multipletests([geneSets[x]["pval"] for x in geneSets],alpha=0.05,method='fdr_bh')
		p_adjust = p_adjust[1].tolist()

		with open("{}neighborhood_analysis/{}.candidates.{}.txt".format(options.Options().qout,sGenesetTag,options.Options().filetag),"w") as OUT:
			OUT.write("GeneSet\tCancer\tpval\tfdr5\tNormalizedAverageAffection\t")
			OUT.write("SwitchingGenes\tOR\tsg\tsng\tnsg\tnsng\n")
			for g,adj_p in zip(geneSets,p_adjust):

				OUT.write("{}\t{}\t{}\t".format(g,options.Options().tag,geneSets[g]["pval"]))
				OUT.write("{}\t{}\t{}\t".format(adj_p,geneSets[g]["affection"], ",".join(geneSets[g]["switchGenes"]) ))
				OUT.write("{}\t{}\t{}\t".format(geneSets[g]["oddsRatio"],geneSets[g]["sg"],geneSets[g]["sng"]))
				OUT.write("{}\t{}\n".format(geneSets[g]["nsg"],geneSets[g]["nsng"]))

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