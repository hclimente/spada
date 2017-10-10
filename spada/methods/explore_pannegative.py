from spada import options
from spada import utils
from methods import method

from collections import Counter
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests
import networkx

class ExplorePannegative(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

		switchPatients = []
		[ switchPatients.extend(z["patients"]) for x,y in self._gene_network.iterate_genes_byPatientNumber() for z in y["isoformSwitches"] ]
		switchPatients = set(switchPatients)

		# read functional mutations
		self.logger.info("Reading functional mutations.")
		self.functionalMutations = self.readFunctionalMutations()

		mutationPatients = []
		[ mutationPatients.extend(self.functionalMutations[x]) for x in self.functionalMutations ]
		mutationPatients = set(mutationPatients)

		self.patients = mutationPatients & switchPatients

		self.logger.info("Reading pathways.")
		self.canonicalPathways = {}
		
		geneSetFile = "{}data/Databases/c2.cp.v4.0.entrez.gmt".format(options.Options().wd)

		for line in utils.readTable(geneSetFile,header=False):
			geneSet = line[0]
			genes = set(line[2:]) & set(self._gene_network.nodes())
			self.canonicalPathways[geneSet] = genes

	def clean(self):
		utils.cmd("mkdir","-p","{}mutations".format(options.Options().qout))

	def run(self):

		drivers = [ x for x,y in self._gene_network.nodes(data=True) if y["specificDriver"] ]
		drivers = set(drivers) & set(self.functionalMutations)
		driverMutations = dict((d, self.functionalMutations[d]) for d in drivers)
		sortedDrivers = sorted(driverMutations,key=lambda a:len(driverMutations[a]),reverse=True)

		self.logger.info("Writing ordered list of mutated specific drivers")
		with open("{}mutations/driver_mutation_number.txt".format(options.Options().qout),"w") as OUT:
			OUT.write("Tumor\tGeneId\tSymbol\tSamples\n")
			
			for d in sortedDrivers:
				OUT.write("{}\t{}\t".format(options.Options().tag,d))
				OUT.write("{}\t".format(self._gene_network._net.node[d]["symbol"]))
				OUT.write("{}\n".format(",".join(driverMutations[d])))

		self.logger.info("Calculating ME with top drivers to find equivalence.")
		with open("{}mutations/mutual_exclusion_top_drivers.txt".format(options.Options().qout),"w") as OUT:
			OUT.write("Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\t")
			OUT.write("Driver\tDriverSymbol\tPathway\tDistance\tMS\tM\tS\tN\tp.me\n")
			
			for i in range(min(10, len(sortedDrivers))):
				driverName = sortedDrivers[i]
				mutatedSamples = set(self.functionalMutations[driverName]) & self.patients

				self.calculateMEWithDriver(mutatedSamples,driverName,OUT)
		
		self.logger.info("Calculating ME with pannegative patients.")
		for i in range(10):
			topDrivers = sortedDrivers[0:i]
			topDriverMutations = dict((d, self.functionalMutations[d]) for d in topDrivers)

			self.calculateMEForPanNegative(topDriverMutations,i+1)

		self.calculateMEForPanNegative(driverMutations,"all")

	def readFunctionalMutations(self):

		mutFile = "{}data/{}/rawdata/{}_exon_mutation-functional-count_full.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag)

		mutations = {}
		allMuts = []

		for line in utils.readTable(mutFile,header=False):

			patient = line[9]
			tx = line[3].split(";")[1]

			if patient==".":
				continue
			
			try:
				geneID = self._transcript_network._net.node[tx]["gene_id"]
			except KeyError:
				self.logger.debug("Transcript {} not in transcript network.".format(tx))
				continue
			except:
				self.logger.info("Unexpected error:", sys.exc_info()[0])
				raise

			mutations.setdefault(geneID,[])
			
			patient = patient.split(";")[0]

			# get genomic positions
			start = line[7]
			end = line[8]
			
			if (patient,geneID,start,end) in allMuts:
				continue
			
			mutations[geneID].append(patient)
			allMuts.append((patient,geneID,start,end))

		return mutations

	def calculateMEForPanNegative(self,driverMutations,minMuts):

		patientsWithMutation = []
		[ patientsWithMutation.extend(driverMutations[x]) for x in driverMutations ]
		patientsWithMutation = set(patientsWithMutation) & self.patients

		table = []

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network,only_models=True,partialCreation=True):
		
			patientsWithSwitch = set(switchDict["patients"]) & self.patients

			ms,m,s,n = self.getContingencyTable(patientsWithSwitch,patientsWithMutation,self.patients)

			lContingencyTable = [[ms,m],[s,n]]
			OR,pval = fisher_exact(lContingencyTable,alternative="greater")

			table.append({"gene": gene,"ms":ms,"m":m,"s":s,"n":n,
						  "fisher_mutual_exclusion":pval,
						  "nTx":thisSwitch.nTx,"tTx":thisSwitch.tTx })

		p_adj_me = multipletests([ x["fisher_mutual_exclusion"] for x in table ],alpha=0.05,method='fdr_bh')
		p_adj_me = p_adj_me[1].tolist()

		with open("{}mutations/pannegative_mutual_exclusion.top_{}_drivers.txt".format(options.Options().qout,minMuts),"w") as OUT:
			OUT.write("Tumor\tGeneId\tSymbol\t")
			OUT.write("Normal_transcript\tTumor_transcript\tMS\tM\t")
			OUT.write("S\tN\tp.me\tadjp.me\n")
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

				OUT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
				OUT.write("{}\t{}\t{}\t{}\t".format(nTx,tTx,ms,m))
				OUT.write("{}\t{}\t{}\t{}\n".format(s,n,p_me,padj_me))

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

	def calculateMEWithDriver(self,patientsWithMutation,driverId,OUT):

		driverName = self._gene_network._net.node[driverId]["symbol"]
		pathways = [ x for x in self.canonicalPathways if driverId in self.canonicalPathways[x] ]

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(self._transcript_network,only_models=True,partialCreation=True):
			patientsWithSwitch = set(switchDict["patients"]) & self.patients
			pwgene = [ x for x in self.canonicalPathways if gene in self.canonicalPathways[x] and x in pathways ]

			ms,m,s,n = self.getContingencyTable(patientsWithSwitch,patientsWithMutation,self.patients)
			lContingencyTable = [[ms,m],[s,n]]
			OR,pval = fisher_exact(lContingencyTable,alternative="greater")
			try:
				d = len(networkx.shortest_path(self._gene_network._net,gene,driverId))
			except networkx.exception.NetworkXNoPath:
				d = -1

			OUT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
			OUT.write("{}\t{}\t".format(switchDict["nIso"],switchDict["tIso"]))
			OUT.write("{}\t{}\t{}\t{}\t".format(driverId,driverName,",".join(pwgene),d))
			OUT.write("{}\t{}\t{}\t{}\t{}\n".format(ms,m,s,n,pval))