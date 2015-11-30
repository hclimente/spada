#!/soft/devel/python-2.7/bin/python

from biological_entities import protein
from biological_entities import transcript
from libs import options
from libs import utils
from methods import method

class MutationFeatureOverlap(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self.mutations = self.readMutations()

	def clean(self):
		utils.cmd("mkdir","{0}mutations".format(options.Options().qout))
		utils.cmd("rm","{0}mutations/mutation_switch_feature_overlap.txt".format(options.Options().qout))
		utils.cmd("rm","{0}mutations/debug_mutations.txt".format(options.Options().qout))

	def run(self):

		## CALCULATE DOMAIN ALTERATION FREQUENCY IN GENERAL
		# get affection of prosite/pfams by mutations and their frequency
		# in the proteome (only most expressed iso per gene)
		self.studyFeaturesMutationsAndSize()

		## CALCULATE MUTATION FREQUENCY ON SWITCHED DOMAINS
		featureOverlap = []
		
		for gene,info,switchDict,thisSwitch in self._gene_network.iterate_switches_ScoreWise(self._transcript_network,only_models=True,relevance=True,partialCreation=True):

			if gene not in self.mutations:
				continue

			switchMutations = self.getSwitchMutations(self.mutations[gene],thisSwitch.nTranscript,thisSwitch.tTranscript)

			if thisSwitch.domainChange:
				f = self.countMutations(gene,info,switchMutations,"Pfam",thisSwitch)
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

			mutations[geneID].setdefault((start,end),[])
			mutations[geneID][mutInfo].append(mutType)

		return mutations

	def studyFeaturesMutationsAndSize(self):

		featureMutationCounts = {}

		with open("{0}mutations/mutations_enrichment.txt".format(options.Options().qout),"w") as MUT, open("{0}mutations/features_information.txt".format(options.Options().qout),"w") as FT, open("{0}mutations/transcripts_information.txt".format(options.Options().qout),"w") as TXS:

			MUT.write("Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\t")
			MUT.write("TPM\tFrame_Shift_Del\tFrame_Shift_Ins\tIn_Frame_Del\t")
			MUT.write("In_Frame_Ins\tMissense_Mutation\tNonsense_Mutation\tNonstop_Mutation\t")
			MUT.write("Frame_Shift_Del_out\tFrame_Shift_Ins_out\tNonsense_Mutation_out\n")

			FT.write("Cancer\tGene\tSymbol\tTranscript\tAnalysis\t")
			FT.write("Feature\tn\tFeatureLength\tStart\tEnd\n")

			TXS.write("Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength\n")
						
			for gene,info in self._gene_network.iterate_genes_ScoreWise():
				if gene not in self.mutations:
					continue

				# get most-abundant, coding transcript as representative of the gene
				txs = [ (x,i["median_TPM_N"]) for x,i in self._transcript_network.nodes(data=True) if i["gene_id"]==gene and i["proteinSequence"] ]

				if not txs: continue

				maxTpm = max([ x[1] for x in txs ])

				txname = [ x for x,tpm in txs if tpm==maxTpm ][0]
				txinfo = self._transcript_network._net.node[txname]

				tx = transcript.Transcript(txname,txinfo)
				tx.readPfamDomains()
				tx.readProsite()

				proteinLength = len(self._transcript_network._net.node[tx.name]["proteinSequence"])

				# get mutations affecting said transcript
				# mutations = { (genomicPos1,genomicPos2): ([mut1,mut2...],set(cdsPos1,cdsPos2...)), ... }
				mutations = self.getSwitchMutations(self.mutations[gene],tx,None)

				self.printMutationInfo(gene,txname,mutations)

				TXS.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
				TXS.write("{}\t{}\t{}\n".format(txname,maxTpm,proteinLength))

				for a,featType in zip(["prosite","Pfam"],[tx._ptms,tx._pfam]):
					# get the features affected by a mutation
					affectedFeats = self.getFeaturesAffectedByMutation(a,mutations,tx,None,False)
					affectedFeats = affectedFeats[txname]

					for f in affectedFeats:
						i = 1
						for inMuts,allMuts,ratio,featSize,mutationTypes in affectedFeats[f]:
							MUT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
							MUT.write("{}\t{}\t{}\t".format(txname,a,f))
							MUT.write("{}\t{}\t".format(i,maxTpm))
							MUT.write("{}\t".format(mutationTypes["Frame_Shift_Del"]))
							MUT.write("{}\t".format(mutationTypes["Frame_Shift_Ins"]))
							MUT.write("{}\t".format(mutationTypes["In_Frame_Del"]))
							MUT.write("{}\t".format(mutationTypes["In_Frame_Ins"]))
							MUT.write("{}\t".format(mutationTypes["Missense_Mutation"]))
							MUT.write("{}\t".format(mutationTypes["Nonsense_Mutation"]))
							MUT.write("{}\t".format(mutationTypes["Nonstop_Mutation"]))
							MUT.write("{}\t".format(mutationTypes["Frame_Shift_Del_out"]))
							MUT.write("{}\t".format(mutationTypes["Frame_Shift_Ins_out"]))
							MUT.write("{}\n".format(mutationTypes["Nonsense_Mutation_out"]))

							i += 1

					for f in featType:
						i = 1
						for start,end in featType[f]:
							FT.write("{}\t{}\t{}\t".format(options.Options().tag,gene,info["symbol"]))
							FT.write("{}\t{}\t{}\t".format(txname,a,f))
							FT.write("{}\t{}\t{}\t".format(i,end - start,start))
							FT.write("{}\n".format(end))

							i += 1

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

		if alterationType=="Pfam":
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


			mutsOnAnyFeature = set()

			# iterate features and find mutations affecting them
			for f in features:
				for start,end in features[f]:
					
					inMuts = 0
					mutationTypes = {"Frame_Shift_Del": 0,"Frame_Shift_Ins": 0,
									 "In_Frame_Del": 0,"In_Frame_Ins": 0,
									 "Missense_Mutation": 0,"Nonsense_Mutation": 0,
									 "Nonstop_Mutation": 0,"Frame_Shift_Del_out": 0, 
									 "Frame_Shift_Ins_out": 0, "Nonsense_Mutation_out": 0}

					featureProteinRange = set(range(start, end+1 ))
					featSize = float(len(featureProteinRange))/(len(tx.cds)/3)

					for m in mutations[tx.name]:

						thoseMutations = mutations[tx.name][m][0]
						mutProteinRange = mutations[tx.name][m][1]
	
						# overlap between mutations and features
						if mutProteinRange & featureProteinRange:
							inMuts += len(thoseMutations)
							mutsOnAnyFeature.add(m)
							
							for t in thoseMutations:
								mutationTypes.setdefault(t,0)
								mutationTypes[t] += 1

						# check upstream mutations that affect the domain too
						elif set(["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation"]) & set(thoseMutations):
							extendedProteinRange = set(range(1,end+1))

							if mutProteinRange & extendedProteinRange:
								for t in thoseMutations:
									if t in ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation"]:
										mutationTypes.setdefault(t+"_out",0)
										mutationTypes[t+"_out"] += 1

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

			if alterationType=="Pfam":
				a = "interpro"
			elif alterationType=="prosite": 
				a = "prosite"

			for line in utils.readTable("{}structural_analysis/{}_analysis.tsv".format(options.Options().qout,a)):
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

	def printMutationInfo(self,gene,tx,mutations):
		with open("{0}mutations/debug_mutations.txt".format(options.Options().qout),"a") as F:
			#F.write("gene\tstart\tend\tmut\tcdsPosStart\tcdsPosEnd\n")

			# mutations = { (genomicPos1,genomicPos2): ([mut1,mut2...],set(cdsPos1,cdsPos2...)), ... }

			for start,end in mutations[tx]:
				muts = mutations[tx][(start,end)][0]
				cdsPos = sorted(list(mutations[tx][(start,end)][1]))

				F.write("{}\t{}\t{}\t{}\t".format(gene,tx,start,end))
				F.write("{}\t{}\t{}\n".format(";".join(muts),cdsPos[0],cdsPos[-1]))


if __name__ == '__main__':
	pass
