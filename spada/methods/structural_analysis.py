from spada.methods import method

class StructuralAnalysis(method.Method):
	def __init__(self, annotation = 'annotation.pklz', search_pfam = True, 
				 search_prosite = True, search_idr = True):

		method.Method.__init__(self, __name__, annotation)
		self.search_pfam 	= search_pfam
		self.search_prosite = search_prosite
		self.search_idr 	= search_idr

	def run(self):

		self.featureAnalysis()

		if self.search_pfam:
			self.ppiAnalysis()

	def featureAnalysis(self):

		self.logger.info("Feature analysis.")

		if self.search_pfam:
			PFAM = open("pfam_analysis.tsv", "w")
			self.writeDomainsHeader(PFAM)
		if self.search_prosite:
			PROSITE = open("prosite_analysis.tsv", "w")
			self.writeDomainsHeader(PROSITE)
		if self.search_idr:
			IDR = open("idr_analysis.tsv", "w")
			self.writeIDRHeader(IDR)

		for gene,_,thisSwitch in self._genes.switches(self._txs):
			if self.search_pfam:
				pfam_change	= thisSwitch.analyzeDomains("Pfam")
				self.writeDomains(PFAM, 'Pfam', gene, thisSwitch, pfam_change)

			if self.search_prosite:
				prosite_change = thisSwitch.analyzeDomains("Prosite")
				self.writeDomains(PROSITE, 'Prosite', gene, thisSwitch, prosite_change)
			
			if self.search_idr:
				idr_change = thisSwitch.analyzeIDR(0.2)
				self.writeIDR(IDR, gene, thisSwitch, idr_change)

	def ppiAnalysis(self):

		self.logger.info("PPI analysis.")

		with open("ppi_analysis.tsv", "w") as PPI:

			self.writePPIHeader(PPI)

			for gene,info,thisSwitch in self._genes.switches(self._txs):

				DDIchanges = self.analyzeDDIs(thisSwitch)
				self.writePPI(PPI, gene, thisSwitch, DDIchanges)

	def analyzeDDIs(self, thisSwitch):

		DDIchanges = {}
		PPIs = { thisSwitch.ctrl: {}, thisSwitch.case: {} }

		for tx in [thisSwitch.ctrl, thisSwitch.case]:
			ppis = self._txs._net.edges(tx, data=True)

			for tx1, tx2, eDict in ppis:
				pairs = [ (next(iter(x)), next(iter(x))) for x in eDict['ddi'] if len(x) == 1 ] # case for homodomain interactions
				pairs.extend([ tuple(p) for p in eDict['ddi'] if len(p) > 1 ])					# case for heterodomain interactions
				domainsTx1 = self._txs[tx1]['Pfam'].keys()
				domainsTx2 = self._txs[tx2]['Pfam'].keys()

				PPIs[tx1][tx2] = set()
				for dA, dB in pairs:
					if dA in domainsTx1 and dB in domainsTx2:
						PPIs[tx1][tx2].add(dA + "@" + dB)
					if dA in domainsTx2 and dB in domainsTx1:
						PPIs[tx1][tx2].add(dB + "@" + dA)

		partners = set(PPIs[thisSwitch.ctrl].keys()) | set(PPIs[thisSwitch.case].keys())

		for p in partners:

			nDDI = PPIs[thisSwitch.ctrl].get(p, set())
			tDDI = PPIs[thisSwitch.case].get(p, set())
			DDIchanges[p] = { "nDDIs": nDDI - tDDI,
							  "tDDIs": tDDI - nDDI,
							  "bothDDIs": nDDI & tDDI,
							  "what": "Unaffected" }

			if DDIchanges[p]["bothDDIs"]:
				if DDIchanges[p]["nDDIs"] or DDIchanges[p]["tDDIs"]:
					DDIchanges[p]["what"] = "Affected"
			elif DDIchanges[p]["nDDIs"] and not DDIchanges[p]["tDDIs"]:
				DDIchanges[p]["what"] = "Lost_in_cases"
			elif not DDIchanges[p]["nDDIs"] and DDIchanges[p]["tDDIs"]:
				DDIchanges[p]["what"] = "Gained_in_cases"

		return DDIchanges

	def writeDomainsHeader(self, OUT):
		OUT.write("Experiment\tGeneId\tControl_transcript\tCase_transcript\t")
		OUT.write("Feature_type\tFeature\tWhat\tIndex\tControl_start\tControl_end\t")
		OUT.write("Case_start\tCase_end\tControl_MacroScore\tControl_MicroScore\t")
		OUT.write("Control_Jaccard\tCase_MacroScore\tCase_MicroScore\tCase_Jaccard\n")

	def writeDomains(self, OUT, featureType, gene, thisSwitch, changes):
		for c in changes:
			OUT.write("{}\t{}\t{}\t".format(self._genes._name, gene, thisSwitch.ctrl))
			OUT.write("{}\t{}\t{}\t".format(thisSwitch.case, featureType, c["feature"]))
			OUT.write("{}\t{}\t{}\t".format(c["what"], c["index"], c['nStart']))
			OUT.write("{}\t{}\t{}\t".format(c["nEnd"], c["tStart"], c['tEnd']))
			OUT.write("{}\t{}\t{}\t".format(c["nM"], c["nm"], c["nJ"]))
			OUT.write("{}\t{}\t{}\n".format(c["tM"], c["tm"], c["tJ"]))

	def writeIDRHeader(self, OUT):
		OUT.write("Experiment\tGeneId\tControl_transcript\tCase_transcript\t")
		OUT.write("What\tSequence\tStartPos\tEndPos\t")
		OUT.write("microScore\tmacroScore\tJaccard\n")

	def writeIDR(self, OUT, gene, thisSwitch, changes):
		for c in changes:
			OUT.write("{}\t".format( self._genes._name ))
			OUT.write("{}\t{}\t{}\t".format(gene, thisSwitch.ctrl, thisSwitch.case))
			OUT.write("{}\t{}\t{}\t".format(c["what"], c["feature"], c["start"]))
			OUT.write("{}\t{}\t{}\t{}\n".format(c["end"], c["M"], c["m"], c["J"]))

	def writePPIHeader(self, OUT):
		OUT.write("Experiment\tGeneId\tControl_transcript\tCase_transcript\t")
		OUT.write("Other_gene\tOther_symbol\tOther_transcript\tWhat\t")
		OUT.write("#nDDIs\t#tDDIs\t#BothDDIs\tnDDIs\ttDDIs\tBothDDIs\n")

	def writePPI(self, OUT, gene, thisSwitch, DDIchanges):
		for tx, ddis in DDIchanges.items():

			other_gene = self._txs[tx]["gene_id"]

			OUT.write("{}\t{}\t".format( self._genes._name, gene ))
			OUT.write("{}\t{}\t".format( thisSwitch.ctrl, thisSwitch.case ))
			OUT.write("{}\t{}\t".format( other_gene, self._genes[other_gene]["symbol"] ))
			OUT.write("{}\t{}\t".format( tx, ddis["what"] ))
			OUT.write("{}\t{}\t".format( len(ddis["nDDIs"]), len(ddis["tDDIs"]) ))
			OUT.write("{}\t{}\t".format( len(ddis["bothDDIs"]), ";".join(ddis["nDDIs"]) ))
			OUT.write("{}\t{}\n".format( ";".join(ddis["tDDIs"]), ";".join(ddis["bothDDIs"]) ))

if __name__ == '__main__':
	pass
