from libs import options
from libs import utils
from methods import method

class AnnotateSwitches(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self,__name__,gn_network,tx_network)

	def run(self):

		self.logger.info("Reading pancancer annotations.")
		recurrence = self.readRecurrence()
		# me: always false
		#me = self.readMutualExclusion()
		co = self.readCoocurrence()
		mutAffectation = self.readSwitchesAffectingMutated()

		self.logger.info("Annotating tumor-type switches.")
		with open("{}candidateList_driverEvidence.tsv".format(options.Options().qout),"w") as OUT:
			OUT.write("Tumor\tGeneId\tSymbol\tNormal_transcript\t")
			OUT.write("Tumor_transcript\tRecurrence\tAffects_mutated_feature\t")
			OUT.write("Mutual_exclusion\tCoocurrence\tPPI\n")
			for gene,info in self._gene_network.iterate_genes_byPatientNumber():
				for switchDict in info["isoformSwitches"]:
					tag = "{}_{}_{}".format(gene,switchDict["nIso"],switchDict["tIso"])
			
					OUT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
					OUT.write("{}\t{}\t".format(switchDict["nIso"],switchDict["tIso"]))
					OUT.write("{}\t{}\t".format(recurrence[tag],mutAffectation[tag]))
					#OUT.write("{}\t{}\n".format(me[tag],co[tag]))
					if tag in co:
						OUT.write("{}\t{}\t{}\n".format(0,co[tag],0))
					else:
						OUT.write("{}\t{}\t{}\n".format(0,0,0))

	def readRecurrence(self):
		recurrent = {}
		for line in utils.readTable("{}analyses/pancancer/candidateList_recurrence.tsv".format(options.Options().wd)):
			switch = "{}_{}_{}".format(line[0],line[2],line[3])
			recurrent[switch] = int(float(line[5]) < 0.05)

		return(recurrent)

	def readMutualExclusion(self):
		me = {}
		for line in utils.readTable("{}analyses/pancancer/candidateList_mutationME.tsv".format(options.Options().wd)):
			switch = "{}_{}_{}".format(line[0],line[2],line[3])
			#me[switch] = int(line[4])
			me[switch] = 0

		return(me)

	def readCoocurrence(self):
		co = {}
		for line in utils.readTable("{}analyses/pancancer/candidateList_mutationCoocurrence.tsv".format(options.Options().wd)):
			switch = "{}_{}_{}".format(line[0],line[2],line[3])
			co[switch] = int(float(line[8]) < 0.05)

		return(co)

	def readSwitchesAffectingMutated(self):
		mutAffecting = {}
		for line in utils.readTable("{}analyses/pancancer/candidateList_mutatedFeatures.tsv".format(options.Options().wd)):
			switch = "{}_{}_{}".format(line[0],line[2],line[3])
			mutAffecting[switch] = int(line[4])

		return(mutAffecting)

if __name__ == '__main__':
	pass