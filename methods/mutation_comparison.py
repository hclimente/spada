#!/soft/devel/python-2.7/bin/python

from libs import options
from libs import utils
from methods import method

class MutationComparison(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)
		self.mutations = self.readMutations()
		self.patients = []
		[ self.patients.extend(z["patients"]) for x,y in self._gene_network.iterate_genes_ScoreWise() for z in y["isoformSwitches"] ]
		self.patients = list(set(self.patients))

	def run(self):

		table = {}

		for gene,info in self._gene_network.iterate_genes_ScoreWise():
			table[gene] = [[0 for x in range(len(self.patients))] for x in range(2)] 
			patientsWithSwitch = []
			[ patientsWithSwitch.extend(x["patients"]) for x in info["isoformSwitches"] if not x["noise"] ]
			for p in patientsWithSwitch:
				i = self.patients.index(p)
				table[gene][0][i] = 1

			patientsWithMut = self.mutations[gene]
			for p in patientsWithMut:
				i = self.patients.index(p)
				table[gene][1][i] = 1

			informativeCols = []
			for i in range(len(self.patients)):
				if table[gene][0][i]==0 and table[gene][1][i]==0:
					continue
				informativeCols.append(i)

			import pdb
			pdb.set_trace()

	def readMutations(self):
		mutFile = "{0}Data/{1}/Rawdata/{2}_gene_mutation-functional-count_full.txt".format(options.Options().wd,options.Options().inputType,options.Options().tag)

		mutations = {}

		for line in utils.readTable(mutFile,header=False):
			patient = line[9]
			if patient==".":
				continue

			geneID,geneSymbol = self._gene_network.nameFilter(full_name=line[3])
			patient = patient.split(";")[0][:-1]
			mutations.setdefault(geneID,[])
			mutations[geneID].append(patient)

		return mutations

if __name__ == '__main__':
	pass