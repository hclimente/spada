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

	def clean(self):
		utils.cmd("mkdir","{0}mutations".format(options.Options().qout))
		utils.cmd("rm","{0}mutations/mutation_switch_feature_overlap.txt".format(options.Options().qout))

	def run(self):

		## CALCULATE DOMAIN ALTERATION FREQUENCY IN GENERAL
		# get affection of prosite/pfams by switches and mutations
		featSwitchCounts = self.getFeatureSwitchFrequency()
		featMutationCounts = self.getFeatureMutationFrequency()

		# get frequency of feature in the proteome (only most expressed iso per gene)
		featCounts,featSizeCounts,totalProteomeSize = self.getFeatureFrequency()
	
		# print frequency of affection
		self.printFeatureFreq(featSwitchCounts,featMutationCounts,featCounts,featSizeCounts,totalProteomeSize)

		# QUITAR
		exit()

		## CALCULATE MUTATION FREQUENCY ON SWITCHED DOMAINS
		featureOverlap = []
		
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,only_models=True,relevance=True,partialCreation=True):

			if gene not in self.mutations:
				continue

			switchMutations = self.getSwitchMutations(self.mutations[gene],thisSwitch.nTranscript,thisSwitch.tTranscript)

			if thisSwitch.domainChange:
				f = self.countMutations(gene,info,switchMutations,"interpro",thisSwitch)
				featureOverlap.extend(f)
			if thisSwitch.ptmChange:
				f = self.countMutations(gene,info,switchMutations,"prosite",thisSwitch)
				featureOverlap.extend(f)

		self.printOverlap(featureOverlap)

	def readMutations(self):

		mutFile = "{0}Data/{1}/Rawdata/{2}_gene_mutation-functional-count_full.txt".format(options.Options().wd,options.Options().inputType,options.Options().tag)

		mutations = {}

		for line in utils.readTable(mutFile,header=False):

			geneID,geneSymbol = self._gene_network.nameFilter(full_name=line[3])
			mutations.setdefault(geneID,{})

			# get patient
			mutInfo = line[9]

			if mutInfo==".":
				continue

			mutType = mutInfo.split(";")[1]

			# get genomic positions
			start = int(line[7])
			end = int(line[8])

			mutInfo = (start,end)
			
			mutations[geneID].setdefault(mutInfo,[])
			mutations[geneID][mutInfo].append(mutType)

		return mutations

	def getFeatureMutationFrequency(self):

		featureMutationCounts = {}

		for gene,info in self._gene_network.iterate_genes_ScoreWise():
			if gene not in self.mutations:
				continue

			# get most-abundant, coding transcript as representative of the gene
			txs = [ (x,info["median_TPM_N"]) for x,info in self._transcript_network.nodes(data=True) if info["gene_id"]==gene and info["proteinSequence"] ]

			if not txs: continue

			maxTpm = max([ x[1] for x in txs ])

			txname = [ x for x,tpm in txs if tpm==maxTpm ][0]
			txinfo = self._transcript_network._net.node[txname]

			tx = transcript.Transcript(txname,txinfo)

			# get mutations affecting said transcript
			# mutations = { (genomicPos1,genomicPos2): ([mut1,mut2...],set(cdsPos1,cdsPos2...)), ... }
			mutations = self.getSwitchMutations(self.mutations[gene],tx,None)

			for t in ["interpro","prosite"]:
				# get the features affected by a mutation
				affectedFeats = self.getFeaturesAffectedByMutation(t,mutations,tx,None,False)
				affectedFeats = affectedFeats[txname]

				for f in affectedFeats:
					for inMuts,allMuts,ratio,featSize,mutationTypes in affectedFeats[f]:
						featureMutationCounts.setdefault(f,
							{"Frame_Shift_Del": 0,"Frame_Shift_Ins": 0,
							 "In_Frame_Del": 0,"In_Frame_Ins": 0,
							 "Missense_Mutation": 0,"Nonsense_Mutation": 0,
							 "Nonstop_Mutation": 0})
						for t in mutationTypes:
							featureMutationCounts[f][t] += 1

		return featureMutationCounts

	def getFeatureSwitchFrequency(self):

		featSwitchCounts = {}

		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,partialCreation=True,only_models=True):

			thisSwitch.readDeepRelevanceAnalysis(skipIupred=True,skipAnchor=True)

			for featType in [thisSwitch._deep_domain_change,thisSwitch._deep_ptm_change]:
				for c in ["Gained_in_tumor","Lost_in_tumor"]:
					for f in featType[c]:
						featSwitchCounts.setdefault(f,{"Gained_in_tumor": 0,"Lost_in_tumor": 0})
						featSwitchCounts[f][c] += len(switchDict["patients"])
			
		return featSwitchCounts

	def getFeatureFrequency(self):
		featCounts = {}
		featSizeCounts = {}
		totalProteomeSize = 0

		for gene,info in self._gene_network.iterate_genes_ScoreWise():

			txs = [ (x,info["median_TPM_N"]) for x,info in self._transcript_network.nodes(data=True) if info["gene_id"]==gene and info["proteinSequence"] ]

			if not txs: continue

			maxTpm = max([ x[1] for x in txs ])

			txname = [ x for x,tpm in txs if tpm==maxTpm ][0]
			txinfo = self._transcript_network._net.node[txname]

			tx = transcript.Transcript(txname,txinfo)
			tx.readPfamDomains()
			tx.readProsite()
			totalProteomeSize += len(self._transcript_network._net.node[tx.name]["proteinSequence"])

			for featType in [tx._ptms,tx._pfam]:
				for f in featType:
					for start,end in featType[f]:
						featCounts.setdefault(f,0.0)
						featSizeCounts.setdefault(f,0.0)
						featCounts[f] += 1
						featSizeCounts[f] += end - start

		return (featCounts,featSizeCounts,totalProteomeSize)

	def getSwitchMutations(self,geneMutations,nTranscript,tTranscript):

		switchMutations = {}

		# if there is an overlap between mutation and any transcript 
		# CDS, include those mutations
		for tx in [nTranscript,tTranscript]:
			if not tx:
				continue

			txName = tx.name
			switchMutations[txName] = {}

			cds = [ x for x in tx.cds_ordered ]
			for m in geneMutations:

				mutationRegion = set(range(m[0],m[1]))
				
				affectedRegionTx = set(cds) & mutationRegion
				affectedRegionProt = set([ (cds.index(x)/3)+1 for x in affectedRegionTx ])

				if affectedRegionProt:
					switchMutations[txName][m] = (geneMutations[m],affectedRegionProt)

		return switchMutations

	def countMutations(self,gene,info,mutations,alterationType,thisSwitch):

		featureOverlap = []

		switchAffected = self.getFeaturesAffectedBySwitch(alterationType,thisSwitch)
		mutationAffected = self.getFeaturesAffectedByMutation(alterationType,mutations,thisSwitch.nTranscript,thisSwitch.tTranscript,True)

		for tx in switchAffected:
			for f in switchAffected[tx]:
				if switchAffected[tx][f] != len(mutationAffected[tx][f]):
					self.logger.error("Differing number of features for gene {}.".format(gene))

				for i in range(switchAffected[tx][f]):
					inMuts,allMuts,ratio,featSize,mutationTypes = mutationAffected[tx][f][i]

					if thisSwitch.nTx == tx:
						what = "Lost_in_tumor"
					elif thisSwitch.tTx == tx:
						what = "Gained_in_tumor"

					d = { "gene": gene, "info": info, 
						  "inMuts": inMuts, "totalMuts": allMuts,
						  "nTx": thisSwitch.nTx, "tTx": thisSwitch.tTx,
						  "feature": f, "featureType": alterationType, 
						  "ratio": ratio, "what": what, "featureSize": featSize,
						  "domainNumber": i }

					featureOverlap.append(d)

		return featureOverlap

	def getFeaturesAffectedBySwitch(self,alterationType,thisSwitch):

		switchAffected = { thisSwitch.nTx: {}, thisSwitch.tTx: {} }
		switchedFeature = []

		if alterationType=="interpro":
			thisSwitch.readDeepRelevanceAnalysis(skipIupred=True,skipPtm=True,skipAnchor=True)
			switchedFeature = thisSwitch._deep_domain_change

		elif alterationType=="prosite":
			thisSwitch.readDeepRelevanceAnalysis(skipDomain=True,skipIupred=True,skipAnchor=True)
			switchedFeature = thisSwitch._deep_ptm_change

		for tx,change in zip([thisSwitch.tTx,thisSwitch.nTx],["Gained_in_tumor","Lost_in_tumor"]):
			for feat in switchedFeature[change]:
				switchAffected[tx].setdefault(feat,0)
				switchAffected[tx][feat] += 1

		return switchAffected

	def getFeaturesAffectedByMutation(self,alterationType,mutations,nTranscript,tTranscript,onlySwitched):

		mutationAffected = {}

		for tx in [nTranscript,tTranscript]:
			# skip when no transcript is passed
			if not tx:
				continue

			mutationAffected.setdefault(tx.name,{})
			mutationsInFeature = {}

			# get iterated features
			features = {}
			if alterationType=="interpro":
				tx.readPfamDomains()
				features = tx._pfam
			elif alterationType=="prosite":
				tx.readProsite()
				features = tx._ptms

			mutsOnAnyFeature = set()

			# iterate features and find mutations affecting them
			for f in features:
				for start,end in features[f]:
					
					inMuts = 0
					mutationTypes = {}
					featureProteinRange = set(range(start, end+1 ))
					featSize = float(len(featureProteinRange))/(len(tx.cds)/3)

					for m in mutations[tx.name]:

						thoseMutations = mutations[tx.name][m][0]
						mutProteinRange = mutations[tx.name][m][1]
	
						# overlap between mutations and features
						## add ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation"]
						if mutProteinRange & featureProteinRange:
							inMuts += len(thoseMutations)
							mutsOnAnyFeature.add(m)
							
							for t in thoseMutations:
								mutationTypes.setdefault(t,0)
								mutationTypes[t] += 1

					mutationsInFeature.setdefault(f,[])
					mutationsInFeature[f].append((inMuts,featSize,mutationTypes))

			allMuts = sum([ len(mutations[tx.name][x][0]) for x in mutsOnAnyFeature ])

			for f in mutationsInFeature:
				mutationAffected[tx.name].setdefault(f,[])
				for inMuts,featSize,mutationTypes in mutationsInFeature[f]:
					if inMuts == 0 and allMuts == 0:
						ratio = float("nan")
					else:
						ratio = 100 * float(inMuts)/allMuts
					mutationAffected[tx.name][f].append((inMuts,allMuts,ratio,featSize,mutationTypes))

		if onlySwitched:
			filteredMutationAffected = {}

			for tx in [nTranscript,tTranscript]:
				if tx is not None: 
					filteredMutationAffected.setdefault(tx.name,{})

			for line in utils.readTable("{}structural_analysis/{}_analysis.tsv".format(options.Options().qout,alterationType)):
				if nTranscript is not None and line[2] != nTranscript.name:
					continue
				elif tTranscript is not None and line[3] != tTranscript.name:
					continue
				
				if line[4] != "Nothing":
					for x,tx in zip([6,7],[nTranscript,tTranscript]):
						if tx is None:
							continue 

						reps = [ int(y) for y in line[x].split("/") ]

						if reps[0] > reps[1]:
							continue

						f = line[5].replace(" ","_")

						filteredMutationAffected[tx.name].setdefault(f,[])
						filteredMutationAffected[tx.name][f].append(mutationAffected[tx.name][f][reps[0]-1])

			mutationAffected = filteredMutationAffected

		return mutationAffected

	def printFeatureFreq(self,featSwitchCounts,featMutationCounts,featCounts,featSizeCounts,totalProteomeSize):

		totalMuts  		= sum([ featMutationCounts[x][y] for x in featMutationCounts for y in featMutationCounts[x] ])
		totalSwitches 	= sum([ featSwitchCounts[x][y] for x in featSwitchCounts for y in featSwitchCounts[x] ])
		totalCounts 	= sum([ featCounts[x] for x in featCounts ])

		with open("{0}mutations/feature_enrichment.txt".format(options.Options().qout),"w") as OUT:
			OUT.write("Cancer\tFeature\tMutRatio\tSwitchRatio\tMutIn\tAllMuts\t")
			OUT.write("SwitchesInLost\tSwitchesInGain\tAllSwitches\tDomainCount\tAllDomains\t")
			OUT.write("DomainSize\tTotalProteomeSize\tFrame_Shift_Del\tFrame_Shift_Ins\tIn_Frame_Del\t")
			OUT.write("In_Frame_Ins\tMissense_Mutation\tNonsense_Mutation\tNonstop_Mutation\n")

			allDomains = featMutationCounts.keys()
			allDomains.extend(featSwitchCounts.keys())

			for f in set(allDomains):
				if f in featMutationCounts:
					mutIn = sum([ featSwitchCounts[f][x] for x in featSwitchCounts[f] ])
				else:
					mutIn = 0

				if f in featSwitchCounts:
					switchesInLost = featSwitchCounts[f]["Lost_in_tumor"]
					switchesInGain = featSwitchCounts[f]["Gained_in_tumor"]
				else:
					switchesInLost = 0
					switchesInGain = 0

				# if not present in the most abundant isoform add a really low value
				if f in featCounts:
					c = featCounts[f]
					k = featSizeCounts[f]
				else:
					c = 0
					k = 0

				OUT.write("{}\t{}\t{}\t".format(options.Options().tag,f,mutIn))
				OUT.write("{}\t{}\t{}\t".format(totalMuts,switchesInLost,switchesInGain))
				OUT.write("{}\t{}\t{}\t".format(totalSwitches,c,totalCounts))
				OUT.write("{}\t{}\t".format(k,totalProteomeSize))
				OUT.write("{}\t".format(featSwitchCounts[f]["Frame_Shift_Del"]))
				OUT.write("{}\t".format(featSwitchCounts[f]["Frame_Shift_Ins"]))
				OUT.write("{}\t".format(featSwitchCounts[f]["In_Frame_Del"]))
				OUT.write("{}\t".format(featSwitchCounts[f]["In_Frame_Ins"]))
				OUT.write("{}\t".format(featSwitchCounts[f]["Missense_Mutation"]))
				OUT.write("{}\t".format(featSwitchCounts[f]["Nonsense_Mutation"]))
				OUT.write("{}\n".format(featSwitchCounts[f]["Nonstop_Mutation"]))

	def printOverlap(self,featureOverlap):
		with open("{0}mutations/mutation_switch_feature_overlap.txt".format(options.Options().qout),"w") as F:
			F.write("Gene\tSymbol\tCancer\tNormal_transcript\tTumor_transcript\t")
			F.write("What\tFeatureType\tFeature\tDomainNumber\tRatio\tDriver\t")
			F.write("MutationsInFeature\tTotalMutations\tFeatureSize\n")

			for f in featureOverlap:
				F.write("{0}\t{1}\t{2}\t".format(f["gene"],f["info"]["symbol"],options.Options().tag) )
				F.write("{0}\t{1}\t{2}\t".format(f["nTx"],f["tTx"],f["what"]) )
				F.write("{0}\t{1}\t{2}\t".format(f["featureType"],f["feature"],f["domainNumber"]) )
				F.write("{0}\t{1}\t{2}\t".format(f["ratio"],f["info"]["Driver"],f["inMuts"]) )
				F.write("{0}\t{1}\n".format(f["totalMuts"],f["featureSize"]) )


if __name__ == '__main__':
	pass
