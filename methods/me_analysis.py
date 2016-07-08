from libs import options
from libs import utils
from methods import method

from scipy.stats import fisher_exact

class MEAnalysis(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)
		
		switchPatients = []
		[ switchPatients.extend(z["patients"]) for x,y in self._gene_network.iterate_genes_byPatientNumber() for z in y["isoformSwitches"] ]
		switchPatients = set(switchPatients)

		self.logger.info("Reading protein-affecting mutations.")
		self.mutations = self.filterWESMutations()

		mutationPatients = []
		[ mutationPatients.extend(self.mutations[x]) for x in self.mutations ]
		mutationPatients = set(mutationPatients)

		self.patients = mutationPatients & switchPatients

	def clean(self):
		utils.cmd("mkdir","-p","{}mutations".format(options.Options().qout))

	def run(self):

		self.logger.info("Looking for mutual exclusion between switches and WES mutations.")
		self.calculateMEForGene()

		# read raw mutations in the patients with a switch
		mutations = self.filterWGSMutations()

		self.logger.info("Filtering WGS mutations.")
		self.logger.info("Getting co-occurrence candidates.")
		utils.cmd("/soft/R/R-3.2.3/bin/Rscript", 
				  "pipeline/methods/mutation_coocurrence_analysis.R", 
				  "{}candidateList_info.tsv".format(options.Options().qout),
				  "{}mutations/wgs_mutations.txt".format(options.Options().qout),
				  "{}mutations/gene_wgs_mutations_all_switches.txt".format(options.Options().qout))

		self.logger.info("Getting mutual exclusion candidates.")
		utils.cmd("/soft/R/R-3.2.3/bin/Rscript", 
				  "pipeline/methods/mutual_exclusion_analysis.R", 
				  "{}mutations/gene_functional_mutations_all_switches.txt".format(options.Options().qout),
				  "{}mutations/gene_wgs_mutations_all_switches.txt".format(options.Options().qout),
				  "{}candidateList_mutationME.tsv".format(options.Options().qout))

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

	def calculateMEForGene(self):

		# search single gene bias
		table = []
		mutations = self.mutations

		with open("{}mutations/gene_functional_mutations_all_switches.txt".format(options.Options().qout),"w") as OUT:
			OUT.write("Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\t")
			OUT.write("MS\tM\tS\tN\tp\tOR\tSwitched\tMutated\n")

			for gene in self._gene_network.nodes():

				info = self._gene_network._net.node[gene]

				# get mutated patients
				if gene not in mutations:
					patientsWithMutation = set()
				else:
					patientsWithMutation = set(mutations[gene]) & self.patients

				# get switched patients
				dummySwitch = { "nIso": None, "tIso": None, "patients": []}
				switches = [dummySwitch] + info["isoformSwitches"]

				for switchDict in switches:

					patientsWithSwitch = set(switchDict["patients"]) & self.patients
					nIso = switchDict["nIso"]
					tIso = switchDict["tIso"]

					ms,m,s,n = self.getContingencyTable(patientsWithSwitch,patientsWithMutation,self.patients)
					lContingencyTable = [[ms,m],[s,n]]
					OR,pval = fisher_exact(lContingencyTable)

					OUT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
					OUT.write("{}\t{}\t".format(nIso,tIso))
					OUT.write("{}\t{}\t{}\t{}\t".format(ms,m,s,n))
					OUT.write("{}\t{}\t".format(pval,OR))
					OUT.write("{}\t".format(",".join(patientsWithSwitch)))
					OUT.write("{}\n".format(",".join(patientsWithMutation)))

	def filterWESMutations(self):

		mutFile = "{}data/{}/rawdata/{}_exon_mutation-functional-count_full.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag)

		mutations = {}
		allMuts = []

		with open("{}mutations/wes_mutations.txt".format(options.Options().qout),"w") as WES:
			WES.write("Tumor\tGene\tSymbol\tTranscript\tPatient\t")
			WES.write("Start\tEnd\tType\tMedianExpression\n")
			for line in utils.readTable(mutFile,header=False):

				patient = line[9]
				tx = line[3].split(";")[-1]

				if patient=="." or tx not in self._transcript_network.nodes():
					continue

				geneID = self._transcript_network._net.node[tx]["gene_id"]
				symbol = self._gene_network._net.node[geneID]["symbol"]
				medianExpression = self._transcript_network._net.node[tx]["median_TPM_N"]

				mutations.setdefault(geneID,[])
				
				patient,mutType = patient.split(";")

				# get genomic positions
				start = int(line[7])
				end = int(line[8])
				
				if (patient,geneID,start,end) not in allMuts:
					mutations[geneID].append(patient)
					allMuts.append((patient,geneID,start,end))

				WES.write("{}\t{}\t{}\t".format(options.Options().tag,geneID,symbol))
				WES.write("{}\t{}\t{}\t".format(tx,patient,start))
				WES.write("{}\t{}\t{}\n".format(end,mutType,medianExpression))

		return(mutations)

	def filterWGSMutations(self):

		switchPatients = set()
		for x,y in self._gene_network._net.nodes(data=True):
			for z in y["isoformSwitches"]:
				for p in z["patients"]:
					switchPatients.add(p[:4])

		mutFile = "{}data/Databases/WGS_505_TCGA/mutations.genes.{}.tsv".format(options.Options().wd,options.Options().tag)

		with open("{}mutations/wgs_mutations.txt".format(options.Options().qout),"w") as OUT:
			OUT.write("Tumor\tGene\tSymbol\tPatient\tPosition\tReference\tVariant\n")

			for line in utils.readTable(mutFile):

				patient = line[1].split("-")[2]

				# skip if variant frequency too low or bSNP138 overlap
				if float(line[8]) < 0.20 or float(line[5]) == 1:
					continue
				# skip if the patient is not among the studied patients
				elif patient not in switchPatients:
					continue

				gene = line[10]
				pos = line[4]
				ref_allele = line[6]
				var_allele = line[7]
				if gene in self._gene_network.nodes():
					symbol = self._gene_network._net.node[gene]["symbol"]
				else:
					symbol = "??"

				OUT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,symbol))
				OUT.write("{}\t{}\t{}\t".format(patient,pos,ref_allele))
				OUT.write("{}\n".format(var_allele))

if __name__ == '__main__':
	pass