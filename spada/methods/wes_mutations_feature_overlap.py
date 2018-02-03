from spada.biological_entities import protein
from spada.biological_entities import transcript
from spada.methods import method
from spada import utils

import math

class WESMutationsFeatureOverlap(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self.functionalMutations = self.readFunctionalMutations()

		# txname="uc001kjc.3"
		# txinfo=self._transcript_network._net.node[txname]
		# tx = transcript.Transcript(txname,txinfo)
		# txMuts = self.getTxsMutations(tx)
		# affectedFeats = self.getFeatureMutations(tx,txMuts,"Pfam")

		# import pdb
		# pdb.set_trace()

		utils.cmd("mkdir","{}mutations".format(options.Options().qout))

		self.TXS_ALL = open("{}mutations/proteome_information.txt".format(options.Options().qout),"w")
		self.TXS_ALL.write("Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength\tasEvidence\n")

		self.TXS_SWT = open("{}mutations/switch_information.txt".format(options.Options().qout),"w")
		self.TXS_SWT.write("Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength\tasEvidence\n")

		self.MUT_ALL = open("{}mutations/proteome_mutations.txt".format(options.Options().qout),"w")
		self.MUT_ALL.write("Cancer\tGene\tSymbol\tTranscript\tAnalysis\t")
		self.MUT_ALL.write("Feature\tn\tType\tPatient\n")

		self.FT_ALL = open("{}mutations/proteome_features.txt".format(options.Options().qout),"w")
		self.FT_ALL.write("Cancer\tGene\tSymbol\tTranscript\tAnalysis\t")
		self.FT_ALL.write("Feature\tn\tFeatureLength\tStart\tEnd\n")

		self.MUT_SWT = open("{}mutations/switch_mutations.txt".format(options.Options().qout),"w")
		self.MUT_SWT.write("Cancer\tGene\tSymbol\tTranscript\tAnalysis\t")
		self.MUT_SWT.write("Feature\tn\tType\tPatient\n")

		self.FT_SWT = open("{}mutations/switch_features.txt".format(options.Options().qout),"w")
		self.FT_SWT.write("Cancer\tGene\tSymbol\tTranscript\tAnalysis\t")
		self.FT_SWT.write("Feature\tn\tFeatureLength\tStart\tEnd\n")

	def clean(self):
		utils.cmd("mkdir","-p","{}mutations".format(options.Options().qout))

	def run(self):

		for gene,info in self._gene_network.iterate_genes_byPatientNumber(onlySplicedGenes=False,onlyExpressedGenes=False):

			allTxs = set(info["expressedTxsN"]) | set(info["expressedTxsT"])
			asEvidence = bool(info["isoformSwitches"]) | (len(allTxs) >= 2)
			expressed = bool(allTxs)

			## CALCULATE DOMAIN ALTERATION FREQUENCY IN GENERAL FOR
			# get affection of prosite/pfams by mutations and their frequency
			# in the proteome (only most expressed iso per gene)
			if expressed:
				tx = self.getMostAbundantTx(gene,info,asEvidence)
				if tx:
					self.getProteinMutations(gene,info,tx,self.MUT_ALL,self.FT_ALL)

			## CALCULATE MUTATION FREQUENCY ON SWITCHED DOMAINS
			if not info["isoformSwitches"]: continue

			usedTxs = set()
			for switchDict in info["isoformSwitches"]:
				if switchDict["noise"]: continue

				thisSwitch = self._gene_network.createSwitch(switchDict,self._transcript_network,True)

				thisSwitch.readDeepRelevanceAnalysis(skipIupred=True,skipAnchor=True)

				for tx in [thisSwitch.tTranscript,thisSwitch.nTranscript]:
					# do not write the same information twice
					if tx in usedTxs:
						continue

					protein = self._transcript_network._net.node[tx.name]["proteinSequence"]
					if protein:
						proteinLength = len(protein)
					else:
						proteinLength = 0
					tpm = self._transcript_network._net.node[tx.name]["median_TPM_N"]

					self.TXS_SWT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
					self.TXS_SWT.write("{}\t{}\t{}\t{}\n".format(tx.name,tpm,proteinLength,asEvidence))

					self.getProteinMutations(gene,info,tx,self.MUT_SWT,self.FT_SWT)
					usedTxs.add(tx)

		self.TXS_ALL.close()
		self.TXS_SWT.close()
		self.MUT_ALL.close()
		self.FT_ALL.close()
		self.MUT_SWT.close()
		self.FT_SWT.close()

		utils.cmd('Rscript', 'pipeline/methods/mutated_features_analysis.R', options.Options().qout)

	def readFunctionalMutations(self):

		mutFile = "{}data/{}/rawdata/{}_exon_mutation-functional-count_full.txt".format(options.Options().wd,options.Options().annotation,options.Options().tag)

		mutations = {}
		allMuts = []

		for line in utils.readTable(mutFile,header=False):

			mutInfo = line[9]
			tx = line[3].split(";")[1]

			if mutInfo=="." or tx not in self._transcript_network.nodes():
				continue

			geneID = self._transcript_network._net.node[tx]["gene_id"]

			x = mutInfo.split(";")
			patient = x[0][:-1]
			mutType = x[1]

			# get genomic positions
			start = int(line[7])
			end = int(line[8])

			if (patient,geneID,start,end) in allMuts:
				continue

			mutations.setdefault(geneID,{})
			mutations[geneID].setdefault((start,end),[])
			mutations[geneID][(start,end)].append((patient,mutType))
			allMuts.append((patient,geneID,start,end))

		return mutations

	def getProteinMutations(self,gene,info,tx,MUTS,FTS):

		txMuts = self.getTxsMutations(tx)

		for a,featType in zip(["prosite","Pfam"],[tx._ptms,tx._pfam]):
			# get the features affected by a mutation
			affectedFeats = self.getFeatureMutations(tx,txMuts,a)

			for f in affectedFeats:
				i = 1
				for m in affectedFeats[f]:
					for patient,mutType in m:
						MUTS.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
						MUTS.write("{}\t{}\t{}\t".format(tx.name,a,f))
						MUTS.write("{}\t{}\t{}\n".format(i,mutType,patient))
					i += 1

			# only count domains from that isoform if can be altered by splicing
			for f in featType:
				i = 1
				for start,end in featType[f]:
					FTS.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
					FTS.write("{}\t{}\t{}\t".format(tx.name,a,f))
					FTS.write("{}\t{}\t{}\t".format(i,end - start,start))
					FTS.write("{}\n".format(end))

					i += 1

	def getTxsMutations(self,tx):

		txMutations = {}

		gene = self._transcript_network._net.node[tx.name]["gene_id"]

		if gene in self.functionalMutations:
			geneMutations = self.functionalMutations[gene]
		else:
			geneMutations = []

		# if there is an overlap between mutation and any transcript
		# CDS, include those mutations
		cds = [ x for x in tx.cds_ordered ]
		for m in geneMutations:
			mutationRegion = set(range(m[0],m[1]))

			affectedRegionTx = set(cds) & mutationRegion
			affectedRegionProt = set([ math.floor((cds.index(x)/3)+1) for x in affectedRegionTx ])

			if affectedRegionProt:
				txMutations[m] = (geneMutations[m],affectedRegionProt)

		return txMutations

	#def getFeaturesAffectedByMutation(self,alterationType,mutations,nTranscript,tTranscript,onlySwitched):
	def getFeatureMutations(self,tx,txMutations,alterationType):

		featureMutation = {}

		# get iterated features
		features = {}
		if alterationType=="Pfam":
			tx.readPfamDomains()
			features = tx._pfam
		elif alterationType=="prosite":
			tx.readProsite()
			features = tx._ptms
		elif alterationType=="iupred":
			pass
			#iso = protein.Protein(tx.name, self._transcript_network._net.node[tx.name])
			#iso.readIupred("long")
			#segments = iso.getSegments("disordered",minLength=5,gap=2)

			#range = [ r.num for s in segments for r in s ]

		# iterate features and find mutations affecting them
		for f in features:
			for start,end in features[f]:

				featureMutation.setdefault(f,[])
				featureProteinRange = set(range(start, end+1 ))
				affectingMuts = []

				for m in txMutations:

					thoseMutations = txMutations[m][0]
					mutProteinRange = txMutations[m][1]

					# overlap between mutations and features
					if mutProteinRange & featureProteinRange:
						affectingMuts.extend(thoseMutations)

					# check upstream mutations that affect the domain too
					if set(["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation"]) & set([ x[1] for x in set(thoseMutations) ]):
						extendedProteinRange = set(range(1,start))
						if mutProteinRange & extendedProteinRange:
							for t in thoseMutations:
								if t[1] in ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation"]:
									affectingMuts.append((t[0],t[1]+"_out"))

				if affectingMuts:
					featureMutation[f].append(affectingMuts)

		return featureMutation

	def getMostAbundantTx(self,gene,info,asEvidence):
		tx = None

		# get most-abundant, coding transcript as representative of the gene
		txs = [ (x,i["median_TPM_N"]) for x,i in self._transcript_network.nodes(data=True) if i["gene_id"]==gene and i["proteinSequence"] ]

		if txs:

			maxTpm = max([ x[1] for x in txs ])

			# only if there is any expressed isoform
			if maxTpm >= 0.1:
				txname = [ x for x,tpm in txs if tpm==maxTpm ][0]
				txinfo = self._transcript_network._net.node[txname]
				tx = transcript.Transcript(txname,txinfo)
				proteinLength = len(self._transcript_network._net.node[tx.name]["proteinSequence"])

				self.TXS_ALL.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
				self.TXS_ALL.write("{}\t{}\t{}\t{}\n".format(txname,maxTpm,proteinLength,asEvidence))

		return tx

if __name__ == '__main__':
	pass
