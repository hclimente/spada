from libs import options
from libs import utils
from methods import method

from scipy.stats import fisher_exact

class MEGenesetAnalysis(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)

		# get list of all expressed genes
		self.allGenes = set([ x for x,y in self._gene_network.nodes(data=True) if y["expressedTxsNormal"] or y["expressedTxsTumor"] ])

		self.logger.info("Reading canonical pathways.")
		self.sets = {}
		geneSetFile = "{}data/Databases/c2.cp.v4.0.entrez.gmt".format(options.Options().wd)

		setNames = []
		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = set(line[2:]) & self.allGenes
			self.sets[geneSet] = genes
			setNames.append(geneSet)

		i = options.Options().parallelRange2
		self.geneset = self.sets[setNames[i]]

		self.logger.info("Reading protein-affecting mutations.")
		self.mutations = self.readFilteredWESMutations()

		# only consider the patients with an alteration in that geneset
		mutationPatients = []
		[ mutationPatients.extend(self.mutations[x]) for x in self.mutations ]
		mutationPatients = set(mutationPatients)

		# read mutations in the patients with a switch in that geneset
		switchPatients = []
		[ switchPatients.extend(z["patients"]) for x,y in self._gene_network.iterate_genes_byPatientNumber() for z in y["isoformSwitches"] and x in self.geneset ]
		switchPatients = set(switchPatients)

		self.patients = mutationPatients & switchPatients

	def clean(self):
		utils.cmd("mkdir","-p","{}mutations/genesets".format(options.Options().qout))

	def run(self):

		self.logger.info("Looking for mutual exclusion between switches and WES mutations at geneset.")
		self.calculateMEForSet(self.geneset,True)
		self.calculateMEForSet(self.geneset,False)

	def calculateMEForSet(self,geneset,onlyDrivers):

		if onlyDrivers:
			tag = "_onlyDrivers"
			iteratedGenes = [ x for x in self.sets[geneset] if x in self._gene_network.nodes() and self._gene_network._net.node[x]["specificDriver"] ]
		else:
			tag = ""
			iteratedGenes = [ x for x in self.sets[geneset] if x in self._gene_network.nodes() ]

		with open("{}mutations/genesets/{}{}.txt".format(options.Options().qout,geneset,tag),"w") as OUT:
			for testedGene in self.sets[geneset]:
				if testedGene not in self._gene_network.nodes():
					continue

				info = self._gene_network._net.node[testedGene]
				patientsWithMutation = []
				for g in [ x for x in iteratedGenes if x in self.mutations ]:
					patientsWithMutation.extend(self.mutations[g])

				patientsWithMutation = set(patientsWithMutation) & self.patients

				patientsWithSwitch = []
				nTxs = []
				tTxs = []

				for switchDict in info["isoformSwitches"]:
					if not switchDict["noise"] and switchDict["model"] and switchDict["functional"]:
						patientsWithSwitch.extend(switchDict["patients"])
						nTxs.append(switchDict["nIso"])
						tTxs.append(switchDict["tIso"])

				patientsWithSwitch = set(patientsWithSwitch) & self.patients
				nTxs = ",".join(nTxs)
				tTxs = ",".join(tTxs)

				ms,m,s,n = self.getContingencyTable(patientsWithSwitch,patientsWithMutation,self.patients)
				lContingencyTable = [[ms,m],[s,n]]
				OR,pval = fisher_exact(lContingencyTable)

				genes = ",".join([ self._gene_network._net.node[x]["symbol"] for x in iteratedGenes])

				OUT.write("{}\t{}\t".format(options.Options().tag,gene))
				OUT.write("{}\t{}\t{}\t".format(info["symbol"],nTxs,tTxs))
				OUT.write("{}\t{}\t{}\t{}\t{}\t".format(hallmark,ms,m,s,n))
				OUT.write("{}\t{}\t{}\n".format(pval,OR,genes))

	def readFilteredWESMutations(self):

		mutations = {}
		allMuts = []

		for line in utils.readTable("{}mutations/wes_mutations.txt".format(options.Options().qout)):
			tumor 				= line[0]
			geneID 				= line[1]
			symbol 				= line[2]
			transcript 			= line[3]
			patient 			= line[4]
			start 				= line[5]
			end 				= line[6]
			mutType				= line[7]
			medianExpression 	= line[8]

			if geneID not in self.geneset:
				continue

			if (patient,geneID,start,end) not in allMuts:
				mutations.setdefault(geneID,[])
				mutations[geneID].append(patient)
				allMuts.append((patient,geneID,start,end))

		return(mutations)

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

if __name__ == '__main__':
	pass