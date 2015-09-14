#!/soft/devel/python-2.7/bin/python

from biological_entities import protein
from biological_entities import transcript
from libs import options
from libs import utils
from methods import method

import os
import scipy

class MutationFeatureOverlap(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self.mutations = self.readMutations()

		# with open("{0}/domain_mutation.txt".format(options.Options().qout)) as OUT:
		# 	OUT.write("Domain\tScore\n")
		# 	for i in range(len(self.domainMutationFrequency)):
		# 		OUT.write("{0}\t{1}\n".format(self.domainMutationFrequency[i][0],self.domainMutationFrequency[i][1]))

	def clean(self):
		utils.cmd("mkdir","{0}mutations".format(options.Options().qout))
		utils.cmd("rm","{0}mutations/mutation_switch_feature_overlap.txt".format(options.Options().qout))

	def run(self):

		domainSwitchFrequency,totalChanges = self.getDomainSwitchFrequency()
		domainMutationFrequency,totalMuts = self.getDomainMutationFrequency()
		domainFrequency = self.getDomainFrequency()
	
		self.printResults(domainSwitchFrequency,domainMutationFrequency,domainFrequency,totalMuts,totalChanges)

		featureOverlap = []
		
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,only_models=True,relevance=True,partialCreation=True):

			if gene not in self.mutations:
				continue

			switchMutations = self.getSwitchMutations(self.mutations[gene],thisSwitch.nTranscript,thisSwitch.tTranscript)

			if thisSwitch.domainChange:
				f = self.countMutations(gene,info,switchMutations,"Pfam_domain",thisSwitch)
				featureOverlap.extend(f)
			if thisSwitch.disorderChange:
				f = self.countMutations(gene,info,switchMutations,"disorder",thisSwitch)
				featureOverlap.extend(f)
			if thisSwitch.anchorChange:
				f = self.countMutations(gene,info,switchMutations,"anchor",thisSwitch)
				featureOverlap.extend(f)
			if thisSwitch.ptmChange:
				f = self.countMutations(gene,info,switchMutations,"Prosite_pattern",thisSwitch)
				featureOverlap.extend(f)

		with open("{0}mutations/mutation_switch_feature_overlap.txt".format(options.Options().qout),"w") as F:
			# F.write("Gene\tSymbol\tCancer\tNormal_transcript\tTumor_transcript\t")
			# F.write("What\tFeatureType\tFeature\tRatio\tDriver\t")
			# F.write("MutationsInFeature\tTotalMutations\tFeatureSize\n")

			for f in featureOverlap:
				F.write("{0}\t{1}\t{2}\t".format(f["gene"],f["info"]["symbol"],options.Options().tag) )
				F.write("{0}\t{1}\t{2}\t".format(f["nTx"],f["tTx"],f["what"]) )
				F.write("{0}\t{1}\t{2}\t".format(f["featureType"],f["feature"],f["ratio"]) )
				F.write("{0}\t{1}\t{2}\t".format(f["info"]["Driver"],f["inMuts"],f["totalMuts"]) )
				F.write("{0}\n".format(f["featureSize"]) )

	def readMutations(self):

		mutFile = "{0}Data/{1}/Rawdata/{2}_gene_mutation-functional-count_full.txt".format(options.Options().wd,options.Options().inputType,options.Options().tag)

		mutations = {}

		for line in utils.readTable(mutFile,header=False):

			geneID,geneSymbol = self._gene_network.nameFilter(full_name=line[3])
			mutations.setdefault(geneID,[])

			# get patient
			patient = line[9]

			if patient==".":
				continue

			# get genomic positions
			start = int(line[7])
			end = int(line[8])

			mutInfo = {"start": start, "end": end}
			
			mutations[geneID].append(mutInfo)

		return mutations

	def getDomainMutationFrequency(self):

		totalMuts = 0
		domainMutationFrequency = {}

		for gene,info in self._gene_network.iterate_genes_ScoreWise():
			if gene not in self.mutations:
				continue

			txs = [ (x,info["median_TPM_N"]) for x,info in self._transcript_network.nodes(data=True) if info["gene_id"]==gene ]

			maxTpm = max([ x[1] for x in txs ])

			txname = [ x for x,tpm in txs if tpm==maxTpm ][0]
			txinfo = self._transcript_network._net.node[txname]

			tx = transcript.Transcript(txname,txinfo)

			mutations = self.getSwitchMutations(self.mutations[gene],tx,None)
			affectedDomains = self.getFeaturesAffectedByMutation("Pfam_domain",mutations,tx,None)
			affectedDomains = affectedDomains[txname]

			for dom in affectedDomains:
				domainMutationFrequency.setdefault(dom,0)
				domainMutationFrequency[dom] += affectedDomains[dom][0]
				totalMuts += affectedDomains[dom][0]

		return (domainMutationFrequency,totalMuts)

	def getDomainSwitchFrequency(self):

		totalChanges = 0
		domainSwitchFrequency = {}

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,partialCreation=True,removeNoise=True):

			thisSwitch.readDeepRelevanceAnalysis(skipIupred=True,skipAnchor=True,skipPtm=True)

			for change in ["Gained_in_tumor","Lost_in_tumor"]:
				for dom in thisSwitch._deep_domain_change[change]:
					domainSwitchFrequency.setdefault(dom,0)
					domainSwitchFrequency[dom] += 1
					totalChanges += 1

		return (domainSwitchFrequency,totalChanges)

	def getDomainFrequency(self):
		domainFrequency = {}
		totalCounts = 0

		for gene,info in self._gene_network.iterate_genes_ScoreWise():

			txs = [ (x,info["median_TPM_N"]) for x,info in self._transcript_network.nodes(data=True) if info["gene_id"]==gene ]

			maxTpm = max([ x[1] for x in txs ])

			txname = [ x for x,tpm in txs if tpm==maxTpm ][0]
			txinfo = self._transcript_network._net.node[txname]

			tx = transcript.Transcript(txname,txinfo)
			tx.readPfamDomains()

			for domainId in tx._pfam:
				for start,end in tx._pfam[domainId]:
					domainFrequency.setdefault(domainId,0.0)
					domainFrequency[domainId] += 1
					totalCounts += 1

		for domainId in domainFrequency:
			domainFrequency[domainId] = domainFrequency[domainId]/totalCounts

		return domainFrequency

	def getSwitchMutations(self,geneMutations,nTranscript,tTranscript):

		switchMutations = {}

		for tx in [nTranscript,tTranscript]:
			if not tx:
				continue

			txName = tx.name
			switchMutations[txName] = []

			cds = [ x for x in tx.cds_ordered ]
			for mutation in geneMutations:
				
				mutationRegion = range(mutation["start"],mutation["end"])
				mutationRegion = set(mutationRegion)

				affectedRegionTx = set(cds) & mutationRegion
				affectedRegionProt = set([ cds.index(x)/3 for x in affectedRegionTx ])

				if affectedRegionProt:
					switchMutations[txName].append(affectedRegionProt)

		return switchMutations

	def countMutations(self,gene,info,mutations,alterationType,thisSwitch):

		featureOverlap = []

		switchAffected = self.getFeaturesAffectedBySwitch(alterationType,thisSwitch)
		mutationAffected = self.getFeaturesAffectedByMutation(alterationType,mutations,thisSwitch.nTranscript,thisSwitch.tTranscript)

		for tx in switchAffected:
			for feature in switchAffected[tx]:
				if feature not in mutationAffected[tx]:
					continue

				if thisSwitch.nTx == tx:
					what = "Lost_in_tumor"
				elif thisSwitch.tTx == tx:
					what = "Gained_in_tumor"

				d = { "gene": gene, "info": info, 
					  "totalMuts": mutationAffected[tx][feature][1],
					  "inMuts": mutationAffected[tx][feature][0],
					  "nTx": thisSwitch.nTx, "tTx": thisSwitch.tTx,
					  "feature": feature, "featureType": alterationType, 
					  "ratio": mutationAffected[tx][feature][2], "what": what,
					  "featureSize": mutationAffected[tx][feature][3] }

				featureOverlap.append(d)

		return featureOverlap

	def getFeaturesAffectedBySwitch(self,alterationType,thisSwitch):
		switchAffected = { thisSwitch.nTx: [], thisSwitch.tTx: [] }

		if alterationType=="Pfam_domain":
			domainFile = "{0}structural_analysis/interpro_analysis.tsv".format(options.Options().qout)

			for line in utils.readTable(domainFile):
				skipFlag = False
				if (line[2] != thisSwitch.nTx and line[3] != thisSwitch.tTx):
					skipFlag = True
				elif line[4]=="Nothing":
					skipFlag = True
					
				if skipFlag:
					continue

				if line[4] == "Lost_in_tumor":
					tx = thisSwitch.nTx
				elif line[4] == "Gained_in_tumor":
					tx = thisSwitch.tTx
				
				switchAffected[tx].append("{0}:{1}".format(line[6],line[7]).replace(" ","_") )

		# elif alterationType=="disorder":
		# 	disorderFile = "{0}structural_analysis/iupred_analysis.tsv".format(options.Options().qout)

		# 	for line in utils.readTable(disorderFile):
		# 		skipFlag = False

		# 		# skip it its not this switch
		# 		if (line[2] != thisSwitch.nTx and line[3] != thisSwitch.tTx):
		# 			skipFlag = True

		# 		# skip if nothing happens
		# 		elif line[4]=="Nothing":
		# 			skipFlag = True

		# 		# skip if the motif is non significant
		# 		elif int(line[8]):
		# 			skipFlag = True
					
		# 		if skipFlag:
		# 			continue
				
		# 		switchAffected.append(line[5])

		elif alterationType=="anchor":
			thisSwitch.readDeepRelevanceAnalysis(skipDomain=True,skipIupred=True,skipPtm=True)

		elif alterationType=="Prosite_pattern":
			thisSwitch.readDeepRelevanceAnalysis(skipDomain=True,skipIupred=True,skipAnchor=True)

			for tx,change in zip([thisSwitch.tTx,thisSwitch.nTx],["Gained_in_tumor","Lost_in_tumor"]):
				for feat in thisSwitch._deep_ptm_change[change]:
					switchAffected[tx].append(feat)

		return switchAffected

	def getFeaturesAffectedByMutation(self,alterationType,mutations,nTranscript,tTranscript):

		mutationAffected = {}

		for tx in [nTranscript,tTranscript]:
			# skip when no transcript is passed
			if not tx:
				continue

			mutationAffected.setdefault(tx.name,{})

			if alterationType=="Pfam_domain":

				tx.readPfamDomains()

				for domainId in tx._pfam:
					for start,end in tx._pfam[domainId]:

						inMuts = 0.0
						allMuts = 0.0

						featureProteinRange = set(range(start, end+1 ))
						domainSize = float(len(featureProteinRange))/(len(tx.cds)/3)

						for mutation in mutations[tx.name]:
							if mutation & featureProteinRange:
								inMuts += 1
							
							allMuts += 1

						if inMuts == 0 and allMuts == 0:
							p = float("nan")
							ratio = float("nan")
						elif allMuts == 0:
							p = float("nan")
							ratio = float("inf")
						else:
							p = scipy.stats.binom_test(inMuts, n=allMuts,p=domainSize)
							ratio = 100 * inMuts/allMuts
						
						p = scipy.stats.binom_test(inMuts, n=allMuts,p=domainSize)

						mutationAffected[tx.name][domainId] = (inMuts,allMuts,ratio,domainSize)

			elif alterationType=="Prosite_pattern":
				
				tx.readProsite()

				for prositeId in tx._ptms:
					for start,end in tx._ptms[prositeId]:

						inMuts = 0.0
						allMuts = 0.0

						featureProteinRange = set(range(start, end+1 ))
						featureSize = 100 * len(featureProteinRange)/(len(tx.cds)/3)

						for mutation in mutations[tx.name]:
							if mutation & featureProteinRange:
								inMuts += 1
							
							allMuts += 1

						if inMuts == 0 and allMuts == 0:
							ratio = float("nan")
						elif allMuts == 0:
							ratio = float("inf")
						else:
							ratio = 100 * inMuts/allMuts
						
						mutationAffected[tx.name][prositeId] = (inMuts,allMuts,ratio,featureSize)

		return mutationAffected

	def printResults(self,domainSwitchFrequency,domainMutationFrequency,domainFrequency,totalMuts,totalSwitches):
		with open("{0}domain_enrichment.txt".format(options.Options().qout),"w") as OUT:
			OUT.write("Cancer\tDomain\tMutRatio\tSwitchRatio\tMutIn\tMutOut\tSwitchesIn\tSwitchesOut\tDomainFrequency\n")

			allDomains = domainMutationFrequency.keys()
			allDomains.extend(domainSwitchFrequency.keys())

			for domainId in set(allDomains):
				if domainId in domainMutationFrequency:
					mutIn = int(domainMutationFrequency[domainId])
				else:
					mutIn = 0

				if domainId in domainSwitchFrequency:
					switchesIn = domainSwitchFrequency[domainId]
				else:
					switchesIn = 0

				# QUITAR VER SI TIENE SENTIDO QUE NO ESTE EN ESTE DICCIONARIO
				if domainId in domainFrequency:
					freq = domainFrequency[domainId]
				else:
					freq = 0

				mutRatio = mutIn/totalMuts
				switchRatio = float(switchesIn)/totalSwitches

				OUT.write("{0}\t{1}\t{2}\t".format(options.Options().tag,domainId,mutRatio))
				OUT.write("{0}\t{1}\t{2}\t".format(switchRatio,mutIn,int(totalMuts-mutIn)))
				OUT.write("{0}\t{1}\t{2}\n".format(switchesIn,totalSwitches-switchesIn,freq))

if __name__ == '__main__':
	pass
