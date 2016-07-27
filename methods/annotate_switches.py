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
		self.logger.info("Reading pannegative switches.")
		panneg = self.readPannegative()
		self.logger.info("Reading mutual exclusion with drivers.")
		driverMe = self.readDriverME()

		self.logger.info("Annotating tumor-type switches.")
		with open("{}candidateList_driverEvidence.tsv".format(options.Options().qout),"w") as OUT:
			OUT.write("Tumor\tGeneId\tSymbol\tNormal_transcript\t")
			OUT.write("Tumor_transcript\tRecurrence\tAffects_mutated_feature\t")
			OUT.write("PPI\tPannegative\tDriverME\n")
			
			for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_byPatientNumber(
				self._transcript_network,only_models=True,relevance=True,partialCreation=True):

				nTx = switchDict["nIso"]
				tTx = switchDict["tIso"]
				symbol = info["symbol"]
				swt = "{}_{}".format(nTx,tTx)

				ppi.setdefault(swt,0)
				panneg.setdefault(swt,0)
				driverMe.setdefault(swt,0)
			
				OUT.write("{}\t{}\t{}\t{}\t{}\t".format(options.Options().tag,gene,symbol,nTx,tTx))
				OUT.write("{}\t{}\t".format(recurrence[swt],mutAffectation[swt]))
				OUT.write("{}\t{}\t{}\n".format(ppi[swt],panneg[swt],driverMe[swt]))

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
		for line in utils.readTable("{}projects_rg/eporta/Switched_interactions_consensus.txt".format(options.Options().wd), header=False):
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

		for i in range(2,11):
			infile = "{}mutations/pannegative_mutual_exclusion.top_{}_drivers.txt".format(options.Options().qout,i)

			if not os.path.isfile(infile):
				break
			else:
				for line in utils.readTable(infile):
					switch = "{}_{}".format(line[3],line[4])

					p = float(line[10])
					panneg.setdefault(switch,0)

					if p < 0.05:
						panneg[switch] = i

		return(panneg)

	def readDriverME(self):
		driverMe = {}
		infile = "{}mutations/mutual_exclusion_top_drivers.txt".format(options.Options().qout)

		if os.path.isfile(infile):
			for line in utils.readTable(infile):
				switch = "{}_{}".format(line[3],line[4])

				p = float(line[13])
				driverMe.setdefault(switch,0)

				if p < 0.05:
					driverMe[switch] = 1

		return(driverMe)

if __name__ == '__main__':
	pass