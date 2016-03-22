from libs import options
from libs import utils
from methods import method

import fisher
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests

class WESMutations(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)
		
		switchPatients = []
		[ switchPatients.extend(z["patients"]) for x,y in self._gene_network.iterate_genes_byPatientNumber() for z in y["isoformSwitches"] ]
		switchPatients = set(switchPatients)

		self.mutations = {}
		self.mutations['all_mutations'] = self.readMutations('all_mutations')
		self.mutations['functional_mutations'] = self.readMutations('functional_mutations')

		mutationPatients = []
		[ mutationPatients.extend(self.mutations['all_mutations'][x]) for x in self.mutations['all_mutations'] ]
		[ mutationPatients.extend(self.mutations['functional_mutations'][x]) for x in self.mutations['functional_mutations'] ]
		mutationPatients = set(mutationPatients)

		self.patients = mutationPatients & switchPatients

		hallmarks = utils.readGeneset("h.all.v5.0.entrez.gmt")
		# pathways = utils.readGeneset("c2.cp.v4.0.entrez.gmt")

		self.hallmarks = hallmarks.copy()
		# self.hallmarks.update(pathways)

	def clean(self):
		# utils.cmd("rm","-r","{0}mutations".format(options.Options().qout))
		# utils.cmd("mkdir","{0}mutations".format(options.Options().qout))
		utils.cmd("mkdir","-p","{}mutations/hallmark_info".format(options.Options().qout))

	def run(self):

		self.logger.info("Looking for mutual exclusion between switches and mutations.")

		for mutationSet in ["all_mutations","functional_mutations"]:
			for switchSet in ["all_switches","functional_switches"]:
				self.calculateMEForGene(mutationSet,switchSet)
				self.calculateMEForSet(mutationSet,switchSet,onlyDrivers=True)
				self.calculateMEForSet(mutationSet,switchSet,onlyDrivers=False)
				self.calculateMEForPanNegative(mutationSet,switchSet)

	def readMutations(self,mutationType):

		if mutationType == 'all_mutations':
			mutFile = "{}data/{}/rawdata/{}_gene_mutation-count_full.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag)
		elif mutationType == 'functional_mutations':
			mutFile = "{}data/{}/rawdata/{}_gene_mutation-functional-count_full.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag)

		mutations = {}

		for line in utils.readTable(mutFile,header=False):

			geneID,geneSymbol = self._gene_network.nameFilter(full_name=line[3])
			mutations.setdefault(geneID,[])

			patient = line[9]

			if patient==".":
				continue
			
			patient = patient.split(";")[0][:-1]
			
			mutations[geneID].append(patient)

		return mutations

	def getContingencyTable(self,patientsWithSwitch,patientsWithMutation,allPatients):

		mutAndSwitch = 0
		onlySwitch = 0
		onlyMuts = 0
		nothing = 0
		
		for p in allPatients:

			if p in patientsWithSwitch and p in patientsWithMutation:
				mutAndSwitch += 1
			elif p in patientsWithSwitch:
				onlySwitch += 1
			elif p in patientsWithMutation:
				onlyMuts += 1
			else:
				nothing += 1

		return mutAndSwitch,onlyMuts,onlySwitch,nothing

	def calculateMEForGene(self,mutationSet,switchSet):

		# search single gene bias
		table = []
		mutations = self.mutations[mutationSet]

		if switchSet == "functional_switches":
			bOnlyFunctional = True
		elif switchSet == "all_switches":
			bOnlyFunctional = None

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network,only_models=True,relevance=bOnlyFunctional,partialCreation=True):
		
			patientsWithSwitch = set(switchDict["patients"]) & self.patients

			if gene not in mutations:
				patientsWithMutation = set()
			else:
				patientsWithMutation = set(mutations[gene]) & self.patients

			ms,m,s,n = self.getContingencyTable(patientsWithSwitch,patientsWithMutation,self.patients)

			p = fisher.pvalue(ms,m,s,n)
			try:
				j = float(ms)/(ms+m+s)
			except ZeroDivisionError:
				j = 0

			table.append({"gene": gene,"ms":ms,"m":m,"s":s,
						  "n":n,"fisher_overlap":p.right_tail,
						  "fisher_mutual_exclusion":p.left_tail,"jaccard":j,
						  "nTx":thisSwitch.nTx,"tTx":thisSwitch.tTx })

		p_adj_me = multipletests([ x["fisher_mutual_exclusion"] for x in table ],alpha=0.05,method='fdr_bh')
		p_adj_me = p_adj_me[1].tolist()
		p_adj_overlap = multipletests([ x["fisher_overlap"] for x in table ],alpha=0.05,method='fdr_bh')
		p_adj_overlap = p_adj_overlap[1].tolist()

		with open("{0}mutations/gene_{1}_{2}.txt".format(options.Options().qout,mutationSet,switchSet),"w") as OUT:
			for i in range(len(table)):
				gene = table[i]["gene"]
				ms = table[i]["ms"]
				m = table[i]["m"]
				s = table[i]["s"]
				n = table[i]["n"]
				p_me = table[i]["fisher_mutual_exclusion"]
				p_o = table[i]["fisher_overlap"]
				padj_me = p_adj_me[i]
				padj_o = p_adj_overlap[i]
				nTx = table[i]["nTx"]
				tTx = table[i]["tTx"]

				info = self._gene_network._net.node[gene]
				j = table[i]["jaccard"]

				OUT.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(options.Options().tag,gene,info["symbol"],nTx,tTx))
				OUT.write("{0}\t{1}\t{2}\t{3}\t".format(ms,m,s,n))
				OUT.write("{0}\t{1}\t{2}\t{3}\n".format(p_me,padj_me,p_o,padj_o,j))

	def calculateMEForSet(self,mutationSet,switchSet,onlyDrivers):

		mutations = self.mutations[mutationSet]
		table = []

		for hallmark in self.hallmarks:

			if onlyDrivers:
				tag = "_onlyDrivers"
				iteratedGenes = [ x for x in self.hallmarks[hallmark] if x in self._gene_network.nodes() and self._gene_network._net.node[x]["driver"] ]
			else:
				tag = ""
				iteratedGenes = [ x for x in self.hallmarks[hallmark] if x in self._gene_network.nodes() ]

			heatmapDict = {}
			for p in self.patients:
				heatmapDict[p] = dict([self._gene_network._net.node[x]["symbol"],0] for x in self.hallmarks[hallmark] if x in self._gene_network.nodes() )

			for testedGene in self.hallmarks[hallmark]:

				if testedGene not in self._gene_network.nodes():
					continue

				info = self._gene_network._net.node[testedGene]

				patientsWithMutation = []
				for g in [ x for x in iteratedGenes if x in mutations ]:
					patientsWithMutation.extend(mutations[g])
				
				patientsWithMutation = set(patientsWithMutation) & self.patients

				# search single gene bias
				patientsWithSwitch = []
				nTx = "NA"
				tTx = "NA"

				for switchDict in info["isoformSwitches"]:
					if not switchDict["noise"] and switchDict["model"]:
						if switchSet == "functional_switches" and self._gene_network.createSwitch(switchDict,self._transcript_network,True).is_functional:
							patientsWithSwitch.extend(switchDict["patients"])
							nTx = switchDict["nIso"]
							tTx = switchDict["tIso"]
						elif switchSet == "all_switches":
							patientsWithSwitch.extend(switchDict["patients"])
							nTx = switchDict["nIso"]
							tTx = switchDict["tIso"]

				patientsWithSwitch = set(patientsWithSwitch) & self.patients

				for g in [ x for x in iteratedGenes if x in mutations ]:
					for p in [ x for x in mutations[g] if x in self.patients ]:
						heatmapDict[p][self._gene_network._net.node[g]["symbol"]] = 1
				
				for p in patientsWithSwitch:
					if heatmapDict[p][info["symbol"]] == 1:
						heatmapDict[p][info["symbol"]] = 2
					else:
						heatmapDict[p][info["symbol"]] = -1
				
				ms,m,s,n = self.getContingencyTable(patientsWithSwitch,patientsWithMutation,self.patients)
				p = fisher.pvalue(ms,m,s,n)

				table.append({"gene": testedGene,"ms":ms,"m":m,"s":s,"n":n,
			  				  "p_me":p.left_tail, "hallmark": hallmark, 
			  				  "nTx":nTx, "tTx":tTx, 
			  				  "genes": ",".join([ self._gene_network._net.node[x]["symbol"] for x in iteratedGenes]) })

			with open("{}mutations/hallmark_info/{}_{}_{}{}.tsv".format(options.Options().qout,hallmark,mutationSet,switchSet,tag),"w") as HALLMARK_HEATMAP:
				HALLMARK_HEATMAP.write("patient\tgene\talteration\n")
				for p in heatmapDict:
					for g in heatmapDict[p]:
						if heatmapDict[p][g] in [1,2]:
							HALLMARK_HEATMAP.write("{0}\t{1}\tMUT\n".format(p,g))
						if heatmapDict[p][g] in [-1]:
							HALLMARK_HEATMAP.write("{0}\t{1}\tSWITCH\n".format(p,g))

		p_adj_me = multipletests([ x["p_me"] for x in table ],alpha=0.05,method='fdr_bh')
		p_adj_me = p_adj_me[1].tolist()

		with open("{0}mutations/geneset_{1}_{2}{3}.txt".format(options.Options().qout,mutationSet,switchSet,tag),"w") as OUT:
			for i in range(len(table)):
				gene = table[i]["gene"]
				hallmark = table[i]["hallmark"]
				ms = table[i]["ms"]
				m = table[i]["m"]
				s = table[i]["s"]
				n = table[i]["n"]
				p_me = table[i]["p_me"]
				genes = table[i]["genes"]
				padj_me = p_adj_me[i]
				nTx = table[i]["nTx"]
				tTx = table[i]["tTx"]

				info = self._gene_network._net.node[gene]

				OUT.write("{0}\t{1}\t".format(options.Options().tag,gene))
				OUT.write("{0}\t{1}\t{2}\t{3}\t".format(info["symbol"],nTx,tTx,hallmark))
				OUT.write("{0}\t{1}\t{2}\t{3}\t".format(ms,m,s,n))
				OUT.write("{0}\t{1}\t{2}\n".format(p_me,padj_me,genes))

	def calculateMEForPanNegative(self,mutationSet,switchSet):

		mutations = self.mutations[mutationSet]
		drivers = [ x for x,y in self._gene_network.nodes(data=True) if y["driver"] ]
		allPatientsWithAMutation = []
		[ allPatientsWithAMutation.extend(mutations[x]) for x in drivers if x in mutations ]

		if switchSet == "functional_switches":
			bOnlyFunctional = True
		elif switchSet == "all_switches":
			bOnlyFunctional = None

		table = []

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network,only_models=True,relevance=bOnlyFunctional,partialCreation=True):
		
			patientsWithSwitch = set(switchDict["patients"]) & self.patients
			patientsWithMutation = set(allPatientsWithAMutation) & self.patients

			ms,m,s,n = self.getContingencyTable(patientsWithSwitch,patientsWithMutation,self.patients)
			p = fisher.pvalue(ms,m,s,n)

			table.append({"gene": gene,"ms":ms,"m":m,"s":s,
						  "n":n,"fisher_mutual_exclusion":p.left_tail,
						  "nTx":thisSwitch.nTx,"tTx":thisSwitch.tTx })

		p_adj_me = multipletests([ x["fisher_mutual_exclusion"] for x in table ],alpha=0.05,method='fdr_bh')
		p_adj_me = p_adj_me[1].tolist()

		with open("{0}mutations/pannegative_{1}_{2}.txt".format(options.Options().qout,mutationSet,switchSet),"w") as OUT:
			for i in range(len(table)):
				gene = table[i]["gene"]
				ms = table[i]["ms"]
				m = table[i]["m"]
				s = table[i]["s"]
				n = table[i]["n"]
				p_me = table[i]["fisher_mutual_exclusion"]
				padj_me = p_adj_me[i]
				nTx = table[i]["nTx"]
				tTx = table[i]["tTx"]

				info = self._gene_network._net.node[gene]

				OUT.write("{0}\t{1}\t{2}\t".format(options.Options().tag,gene,info["symbol"]))
				OUT.write("{0}\t{1}\t{2}\t{3}\t".format(nTx,tTx,ms,m))
				OUT.write("{0}\t{1}\t{2}\t{3}\n".format(s,n,p_me,padj_me))
				

if __name__ == '__main__':
	pass
