from interface import export_to_MSAnalysis
from interface import out_network
from libs import options
from methods import method
from network import ucsc_gene_network, ucsc_isoform_network

import cPickle

class GetSwitches(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):

		self.calculateSwitches()

		out_network.outputGTF(self._gene_network,self._transcript_network)
		out_network.outCandidateList(self._gene_network,self._transcript_network)
		export_to_MSAnalysis.Export2MSAnalysis().generateFile(self._gene_network,self._transcript_network)

		options.Options().printToFile(initialStep="get-relevant-switches")

	def calculateSwitches(self):

		import drmaa

		s = drmaa.Session()
		s.initialize()

		natSpec = ""
		natSpec += "-q normal -l 'qname=normal|long' "
		natSpec += "-cwd "
		natSpec += "-V "

		self.logger.info("Reading and summarizing input files: computing PSI values and intereplicate agreement.")
		jt = s.createJobTemplate()
		jt.remoteCommand = '/soft/R/R-3.0.0/bin/Rscript'
		jt.args = ['Pipeline/methods/explore_data.r', options.Options().qout,
				   "Data/Input/{0}/{1}/".format(options.Options().inputType,options.Options().tag)]
		jt.joinFiles=True
		jt.nativeSpecification = natSpec
		jt.nativeSpecification += "-N explore "
		jt.nativeSpecification += "-e {0}/esmartas_explore_{1}.txt ".format(options.Options().qout,options.Options().tag)
		jt.nativeSpecification += "-o {0}/osmartas_explore_{1}.txt".format(options.Options().qout,options.Options().tag)
		jobid = s.runJob(jt)

		retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)

		s.deleteJobTemplate(jt)

		self.logger.info("Extracting transcripts with high variance and high expression.")
		allPatients = options.Options().replicates.union(options.Options().unpairedReplicates)
		patientCandidates = []

		for patient in allPatients:
			jt = s.createJobTemplate()
			jt.remoteCommand = '/soft/R/R-3.0.0/bin/Rscript'
			jt.args = ['Pipeline/methods/get_candidates_for_patient.r',
					   options.Options().qout,patient]
			jt.joinFiles=True
			jt.nativeSpecification = natSpec
			jt.nativeSpecification += "-N p{0} ".format(patient)
			jt.nativeSpecification += "-e {0}/esmartas_{1}.txt ".format(options.Options().qout,patient)
			jt.nativeSpecification += "-o {0}/osmartas_{1}.txt".format(options.Options().qout,patient)

			jobid = s.runJob(jt)
			patientCandidates.append(jobid)
			s.deleteJobTemplate(jt)

		for jobid in patientCandidates:
			retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)

		self.logger.info("Filtering switches with clustering measures.")
		
		jt = s.createJobTemplate()
		jt.remoteCommand = '/soft/R/R-3.0.0/bin/Rscript'
		jt.args = ["Pipeline/methods/switch_validation.r",
				   options.Options().qout]
		jt.joinFiles=True
		jt.nativeSpecification = natSpec
		jt.nativeSpecification += "-N val{0} ".format(options.Options().tag)
		jt.nativeSpecification += "-e {0}/esmartas_validate_{1}.txt ".format(options.Options().qout,options.Options().tag)
		jt.nativeSpecification += "-o {0}/osmartas_validate_{1}.txt".format(options.Options().qout,options.Options().tag)
		jobid = s.runJob(jt)

		retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)

		s.deleteJobTemplate(jt)
		s.exit()

	def createGeneNetwork(self):

		self.logger.info("Creating gene network.")

		if options.Options().inputType == "TCGA": 
			self._gene_network = ucsc_gene_network.UCSCGeneNetwork()
		else:
			self.logger.error("Unrecognized input type {0}.".format(options.Options().inputType))
			exit()
		
		self.logger.debug("Reading gene info.")
		self._gene_network.readGeneInfo()
		self._gene_network.importDiffExpression()
		if options.Options().specificDrivers:
			self._gene_network.importSpecificDrivers()
		
		self._gene_network.importCandidates()
		self._gene_network.calculateCompatibilityTable()
		self._gene_network.importKnownInteractions()

		self._gene_network.saveNetwork("geneNetwork.pkl")

	def createTranscriptNetwork(self, recover=False):
			if options.Options().inputType == "TCGA":
				self._transcript_network = ucsc_isoform_network.UCSCIsoformNetwork()
			else:
				self.logger.error("Unrecognized input type {0}.".format(options.Options().inputType))
				exit()

			self.logger.info("Creating transcript network.")
			self._transcript_network.importTranscriptome()
			self._transcript_network.readTranscriptInfo()
			self._transcript_network.saveNetwork("txNetwork.pkl")

if __name__ == '__main__':
	pass