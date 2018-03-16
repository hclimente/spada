from spada.io import io
from spada.methods import method

from itertools import groupby
import numpy as np
from operator import itemgetter

class Summary(method.Method):
	def __init__(self,gn_network,tx_network):
		method.Method.__init__(self, __name__,gn_network,tx_network)

		self.proteinStats = []
		self.exonStats  = []
		self.alternativeSplicingStats  = []

		self.featuresTable = []

	def run(self):

		self.logger.info("Summarizing results.")
		txDict = self._txs.nodes(data=True)

		io.printSwitches(self._genes, self._txs)

		self.logger.info("Taking measures on switches.")
		for gene,info,thisSwitch in self._genes.iterate_switches_byPatientNumber(self._txs, only_models=True, removeNoise=True):
			self.logger.debug("Getting statistics for switch {} - {}.".format(thisSwitch.nTx,thisSwitch.tTx))

			# general protein, switch and gene info
			self.exonOverview(gene, info, thisSwitch)
			self.proteinOverview(txDict, thisSwitch)

			# structural info
			self.changedStructuralFeatures(gene, info, thisSwitch)

		self.printSplicingInfo()
		self.printStructutalInfo()

	def proteinOverview(self,txDict,thisSwitch):

		nIsoLength = 0
		tIsoLength = 0
		nIsoSpecificLength = 0
		tIsoSpecificLength = 0
		nTxLength = 0
		tTxLength = 0
		switch = "{}_{}".format(thisSwitch.nTx, thisSwitch.tTx)

		if thisSwitch.nIsoform:
			nIsoLength = len(thisSwitch.nIsoform.seq)
			nIso = thisSwitch.nIsoform.getSegments('isoform-specific')
			nSp = []
			[nSp.extend(x) for x in nIso]
			nIsoSpecificLength = len(nSp)
		if thisSwitch.tIsoform:
			tIsoLength = len(thisSwitch.tIsoform.seq)
			tIso = thisSwitch.tIsoform.getSegments('isoform-specific')
			tSp = []
			[tSp.extend(x) for x in tIso]
			tIsoSpecificLength = len(tSp)

		self.proteinStats.append((switch,nIsoLength,tIsoLength,nIsoSpecificLength,tIsoSpecificLength))

	def exonOverview(self, gene, info, thisSwitch):
		nTx = thisSwitch.nTranscript
		tTx = thisSwitch.tTranscript

		nSpecificCds = set([ x for x in nTx.cds if x not in tTx.cds ])
		nSpecificUtr = set([ x for x in nTx.utr if x not in tTx.utr ])
		tSpecificCds = set([ x for x in tTx.cds if x not in nTx.cds ])
		tSpecificUtr = set([ x for x in tTx.utr if x not in nTx.utr ])

		for specificCds,specificUtr,cds,origin in zip([nSpecificCds,tSpecificCds],[nSpecificUtr,tSpecificUtr],[nTx.cds,tTx.cds],["nIso","tIso"]):

			specificRegions = sorted(list(specificUtr | specificCds))
			exons = [ set(map(itemgetter(1),g)) for k,g in groupby(enumerate(specificRegions), lambda x: x[0]-x[1]) ]

			for exon in exons:

				exonInfo = {}
				exonInfo["switch"] = "{}_{}_{}".format(gene, thisSwitch.nTx, thisSwitch.tTx)
				exonInfo["length"] = len(exon)
				exonInfo["origin"] = origin

				if exon & specificCds:
					exonicCds = sorted(list(exon & specificCds))

					exonInfo["cdsLength"] = len(exonicCds)
					exonInfo["keepORF"] = True if len(exonicCds)%3==0 else False

					firstPos = exonicCds[0]
					pos = [ i for i,x in enumerate(cds) if x==firstPos ][0]
					exonInfo["position"] = float(pos)/len(cds)
					exonInfo["relativeSize"] = float(exonInfo["cdsLength"])/len(cds)

					if exon & specificUtr:
						exonInfo["role"] = "CDS-UTR"
					else:
						exonInfo["role"] = "CDS"

				elif exon & specificUtr:
					exonInfo["role"] = "UTR"
					exonInfo["cdsLength"] = 0
					exonInfo["keepORF"] = "NA"
					exonInfo["position"] = "NA"
					exonInfo["relativeSize"] = "NA"

				self.exonStats.append(exonInfo)

		txCorresp = thisSwitch.analyzeSplicing()
		orfChange = 0

		for i in range(len(txCorresp)):
			nVersion = txCorresp[i][0]
			tVersion = txCorresp[i][1]

			if nVersion == tVersion:
				tag = "COMMON"
			elif i == 0:
				tag = "BEGINNING"
			elif i == len(txCorresp) - 1:
				tag = "ENDING"
			else:
				tag = "MIDDLE"

			exon = {}

			exon["gene"] = gene
			exon["symbol"] = info["symbol"]
			exon["nTranscript"] = thisSwitch.nTx
			exon["tTranscript"] = thisSwitch.tTx
			exon["nVersion"] = 0 if nVersion is None else len(nVersion)
			exon["tVersion"] = 0 if tVersion is None else len(tVersion)
			exon["tag"] = tag

			orfChange = orfChange + exon["nVersion"]%3 - exon["tVersion"]%3
			if abs(orfChange) >= 3:
				orfChange = orfChange - 3*np.sign(orfChange)
			exon["orfChange"] = orfChange

			self.alternativeSplicingStats.append(exon)

	def printSplicingInfo(self):

		with open("isoform_length.tsv".format(), "w" ) as F:
			F.write("Experiment\tNormal_transcript\tTumor_transcript\tnIsoLength")
			F.write("\ttIsoLength\tnIsoSpecificLength\ttIsoSpecificLength\n")
			for switch,nlen,tlen,nsplen,tsplen in self.proteinStats:
				switchIsoforms = switch.split("_")
				nIso = switchIsoforms[0]
				tIso = switchIsoforms[1]
				F.write("{}\t{}\t{}\t".format(self._genes.tumor,nIso,tIso))
				F.write("{}\t{}\t{}\t{}\n".format(nlen,tlen,nsplen,tsplen))

		with open("exons.tsv", "w" ) as F:
			F.write("Experiment\tSwitch\tOrigin\tType\tLength\tCDSLength\t")
			F.write("CDSRelativeSize\tPosition\tKeepOrf\n")
			for exon in self.exonStats:
				F.write("{}\t".format(self._genes.tumor))
				F.write("{}\t{}\t".format(exon["switch"],exon["role"]))
				F.write("{}\t".format(exon["origin"]))
				F.write("{}\t{}\t".format(exon["length"],exon["cdsLength"]))
				F.write("{}\t{}\t".format(exon["relativeSize"],exon["position"]))
				F.write("{}\n".format(exon["keepORF"]))

		with open("exons_new.tsv", "w" ) as F:
			F.write("Experiment\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\t")
			F.write("Tag\tOrfChange\tnormalSegment\ttumorSegment\n")
			for exon in self.alternativeSplicingStats:
				F.write("{}\t".format(self._genes.tumor))
				F.write("{}\t{}\t".format(exon["gene"],exon["symbol"]))
				F.write("{}\t{}\t".format(exon["nTranscript"],exon["tTranscript"]))
				F.write("{}\t{}\t".format(exon["tag"],exon["orfChange"]))
				F.write("{}\t{}\n".format(exon["nVersion"],exon["tVersion"]))

	def printStructutalInfo(self):

		with open("structural_features.tsv", "w" ) as F:
			F.write("Experiment\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\t")
			F.write("Analysis\tWhatsHappenning\tFeature\tDriver\tASDriver\tDriverType\n")

			for tag,featureDict in self.featuresTable:
				switchElements = tag.split("_")
				gene = switchElements[0]
				symbol = switchElements[1]
				nTx = switchElements[2]
				tTx = switchElements[3]

				for analysis in ["Pfam","idr","Prosite"]:
					for data in featureDict[analysis]:
						F.write("{}\t{}\t".format(self._genes.tumor, gene))
						F.write("{}\t{}\t".format(symbol,nTx))
						F.write("{}\t{}\t".format(tTx,analysis))
						F.write("{}\t{}\t".format(data[1],data[0].replace(" ","_")))
						F.write("{}\n".format(featureDict["driver"]))

	def changedStructuralFeatures(self, gene, info, thisSwitch):

		tag = "{}_{}_{}_{}".format(gene,info["symbol"],thisSwitch.nTx,thisSwitch.tTx)

		switchFeatures = {}

		pfam = []
		prosite = []
		idr = []

		if thisSwitch.isFunctional:
			if thisSwitch._pfamChange:
				for x in thisSwitch.analyzeDomains('Pfam'):
					pfam.append((x['what'], x['feature']))

			if thisSwitch._idrChange:
				for x in thisSwitch.analyzeIDR(0.2):
					idr.append((x['what'], x['feature']))

			if thisSwitch._prositeChange:
				for x in thisSwitch.analyzeDomains('Prosite'):
					prosite.append((x['what'], x['feature']))

		switchFeatures["Pfam"] = pfam
		switchFeatures["idr"] = idr
		switchFeatures["Prosite"] = prosite

		switchFeatures["driver"] = self._genes.isDriver(gene)
		switchFeatures["Functional"] = int(thisSwitch.isFunctional)
		switchFeatures["Model"] = int(thisSwitch.isMain)
		switchFeatures["Noise"] = int(thisSwitch.isNoise)

		self.featuresTable.append((tag,switchFeatures))
