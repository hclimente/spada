from libs import options
from libs import utils
from methods import method

import os

class AnnotateSwitches(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)

		self.proteome = self.readProteome()

	def run(self):

		self.logger.info("Reading pancancer recurrence.")
		recurrence = self.readRecurrence()
		self.logger.info("Reading affection to mutated features.")
		mutAffectation = self.readSwitchesAffectingMutated()
		self.logger.info("Reading PPI effect.")
		ppi = self.readPPIAffection()
		self.logger.info("Reading pannegative switches and mutual exclusion with drivers.")
		panneg = self.readPannegative()

		self.logger.info("Annotating tumor-type switches.")
		with open("{}candidateList_driverEvidence.tsv".format(options.Options().qout),"w") as OUT:
			OUT.write("Tumor\tGeneId\tSymbol\tNormal_transcript\t")
			OUT.write("Tumor_transcript\tRecurrence\t")
			OUT.write("Affects_mutated_feature\tPPI\tPannegative\n")
			
			for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(
				self._transcript_network,only_models=True,relevance=True,partialCreation=True):

				nTx = switchDict["nIso"]
				tTx = switchDict["tIso"]
				symbol = info["symbol"]
				swt = "{}_{}".format(nTx,tTx)

				ppi.setdefault(swt,0)
				panneg.setdefault(swt,0)
			
				OUT.write("{}\t{}\t{}\t{}\t{}\t".format(options.Options().tag,gene,symbol,nTx,tTx))
				OUT.write("{}\t{}\t".format(recurrence[swt],mutAffectation[swt]))
				OUT.write("{}\t{}\n".format(ppi[swt],panneg[swt]))

	def readProteome(self):
		proteome = {}

		for line in utils.readTable("{}mutations/proteome_information.txt".format(options.Options().qout)):
			proteome[line[1]] = line[3]

		return(proteome)

	def readRecurrence(self):
		recurrent = {}
		for line in utils.readTable("{}analyses/pancancer/candidateList_recurrence.tsv".format(options.Options().wd)):
			switch = "{}_{}".format(line[2],line[3])
			recurrent[switch] = int(float(line[5]) < 0.05)

		return(recurrent)

	def readSwitchesAffectingMutated(self):
		mutAffecting = {}
		for line in utils.readTable("{}analyses/pancancer/candidateList_mutatedFeatures.tsv".format(options.Options().wd)):
			switch = "{}_{}".format(line[2],line[3])
			mutAffecting[switch] = int(line[4])

		return(mutAffecting)

	def readPPIAffection(self):
		ppi = {}
		for line in utils.readTable("{}projects_rg/eporta/raw_tables/Switched_interactions_consensus.txt".format(options.Options().wd), header=False):
			switch = "{}_{}".format(line[2],line[3])

			partner = line[4]

			if partner in self.proteome:
				allInteractions = line[6:]
				mostExpressed = [ x for x in allInteractions if self.proteome[partner] in x ]
				if [ x for x in mostExpressed if "Kept" not in x ]:
					ppi[switch] = 1

		return(ppi)

	def readPannegative(self):
		panneg = {}

		infile = "{}mutations/mutual_exclusion_top_drivers.txt".format(options.Options().qout)

		# ask for mutual exclusion with a driver from the same pathway and with at least the aggregation of 3 drivers
		if os.path.isfile(infile):
			for line in utils.readTable(infile):
				switch = "{}_{}".format(line[3],line[4])

				p = float(line[13])
				panneg.setdefault(switch,0)

				if p < 0.05:
					panneg[switch] = 1

		for i in range(3,11):
			infile = "{}mutations/pannegative_mutual_exclusion.top_{}_drivers.txt".format(options.Options().qout,i)

			if not os.path.isfile(infile):
				break
			else:
				for line in utils.readTable(infile):
					switch = "{}_{}".format(line[3],line[4])

					p = float(line[9])
					panneg.setdefault(switch,0)

					if p < 0.05 and panneg[switch] == 1:
						panneg[switch] = i

		# remove those cases where not both conditions are met
		for switch in panneg:
			if panneg[switch] < 3:
				panneg[switch] = 0

		return(panneg)

if __name__ == '__main__':
	pass