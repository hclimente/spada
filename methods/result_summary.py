#!/soft/devel/python-2.7/bin/python

from libs import options
from libs import utils
from methods import method

import pandas as pd

class ResultSummary(method.Method):
	def __init__(self, gn_network, tx_network, gn_subnetwork):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		self.logger.info("Summarizing results.")

		self.generalSwitchInfo()
		self.structuralInfo()

	def generalSwitchInfo(self):

		stats = {}

		stats["hasCds"] = {}
		stats["hasCds"].setdefault("both", 0)
		stats["hasCds"].setdefault("onlyN", 0)
		stats["hasCds"].setdefault("onlyT", 0)
		stats["hasCds"].setdefault("none", 0)

		stats["hasCdsChange"] = {}
		stats["hasCdsChange"].setdefault("y", 0)
		stats["hasCdsChange"].setdefault("n", 0)

		stats["hasUtrChange"] = {}
		stats["hasUtrChange"].setdefault("y", 0)
		stats["hasUtrChange"].setdefault("n", 0)

		stats["isDriver"] = {}
		stats["isDriver"].setdefault("y", 0)
		stats["isDriver"].setdefault("n", 0)

		total = 0
		
		for gene,info,switch in utils.iterate_switches_ScoreWise(self._gene_network):
			nIso = switch.nTranscript
			tIso = switch.tTranscript

			if nIso.cds and tIso.cds: 	stats["hasCds"]["both"]  += 1
			elif nIso.cds: 				stats["hasCds"]["onlyN"] += 1
			elif tIso.cds: 				stats["hasCds"]["onlyT"] += 1
			else: 						stats["hasCds"]["none"]  += 1

			if switch.cds_diff: stats["hasCdsChange"]["y"] += 1
			else: 				stats["hasCdsChange"]["n"] += 1

			if switch.utr_diff: stats["hasUtrChange"]["y"] += 1
			else: 				stats["hasUtrChange"]["n"] += 1

			if info["Driver"]:  stats["isDriver"]["y"]	+=1
			else: 				stats["isDriver"]["n"]	+=1

			total += 1

		with open("{0}resultSummary_switches.tsv".format(options.Options().qout), "w" ) as SUMMARY:
			SUMMARY.write("Has_CDS\tBoth\tOnly_nIso\tOnly_tIso\tNone\n")
			SUMMARY.write("\t{0}\t{1}\t".format(stats["hasCds"]["both"]/total*100,
												stats["hasCds"]["onlyN"]/total*100) )
			SUMMARY.write("{0}\t{1}\n".format(stats["hasCds"]["onlyT"]/total*100,
											  stats["hasCds"]["none"]/total*100) )

			SUMMARY.write("Has_CDS_change\tYes\tNo\n")
			SUMMARY.write("\t{0}\t{1}\n".format(stats["hasCdsChange"]["y"]/total*100,
											  stats["hasCdsChange"]["n"]/total*100) )

			SUMMARY.write("Has_UTR_change\tYes\tNo\n")
			SUMMARY.write("\t{0}\t{1}\n".format(stats["hasUtrChange"]["y"]/total*100,
												stats["hasUtrChange"]["n"]/total*100) )

			SUMMARY.write("Is_driver\tYes\tNo\n")
			SUMMARY.write("{0}\t{1}\n".format(stats["isDriver"]["y"]/total*100,
											  stats["isDriver"]["n"]/total*100) )

	def structuralInfo(self):
		stats = {}

		stats["loopsChange"] = {}
		stats["loopsChange"].setdefault("different", 0)
		stats["loopsChange"].setdefault("same", 0)
		stats["loopsChange"].setdefault("onlyN", 0)
		stats["loopsChange"].setdefault("onlyT", 0)
		stats["loopsChange"].setdefault("noLoops", 0)

		total = 0

		for gene,info,switch in utils.iterate_switches_ScoreWise(self._gene_network):
			nIso = switch.nTranscript
			tIso = switch.tTranscript

			nLoops = tx_network._net.node[nIso.name]["iLoopsFamily"]
			tLoops = tx_network._net.node[tIso.name]["iLoopsFamily"]

			if nLoops or tLoops:
				if nLoops == tLoops: 
					stats["loopsChange"]["same"] += 1
				elif nLoops:
					if tLoops:
						stats["loopsChange"]["different"] += 1
					else: 
						stats["loopsChange"]["onlyN"]	 += 1
				elif tLoops:
					stats["loopsChange"]["onlyT"]+= 1
			else:
				stats["loopsChange"]["noLoops"]	+=1

			total += 1

		with open("{0}resultSummary_structural_loops.tsv".format(options.Options().qout), "w" ) as SUMMARY:
			SUMMARY.write("Mapped_loops\tDifferent\tSame\tOnly_nIso\tOnly_tIso\tNone\n")
			SUMMARY.write("\t{0}\t{1}\t".format(stats["loopsChange"]["different"]/total*100,
												stats["loopsChange"]["same"]/total*100) )
			SUMMARY.write("{0}\t{1}\t".format(stats["loopsChange"]["onlyN"]/total*100,
											  stats["loopsChange"]["onlyT"]/total*100) )
			SUMMARY.write("{0}\n".format(stats["loopsChange"]["noLoops"]/total*100) )


		allMotifs = []
		[ allMotifs.extend(z.disorderChange) for x,y,z in utils.iterate_switches_ScoreWise(self._gene_network) ]
		meanLength = float(sum([ len(x) for x in allMotifs ])/len(allMotifs)) if len(allMotifs) > 0 else float('nan')

		with open("{0}resultSummary_structural_disorder.tsv".format(options.Options().qout), "w" ) as SUMMARY:
			pass

	def neighborhoodInfo(self):
		neighborhoodAnalysis = [ y["neighborhoods"] for x,y,z in utils.iterate_switches_ScoreWise(self._gene_network) ]
		expressedGenes = [ (x,y) for x,y in self._gene_network.nodes(data=True) if y["ExpressedTranscripts"]]

		patients = options.Options().replicates
		if options.Options().unpairedReplicates:
			patients = options.Options().unpairedReplicates

		for analysis in neighborhoodAnalysis:

			geneSets = []
			[ geneSets.extend(y["neighborhoods"][analysis]) for x,y,z in utils.iterate_switches_ScoreWise(self._gene_network) ]
			geneSets = set(geneSets)

			with open("{0}resultSummary_neighborhoods_genes_{1}.tsv".format(options.Options().qout,analysis), "w" ) as SUMMARY:
				SUMMARY.write("GeneSet")
				[ SUMMARY.write("\t{0}".format(gene)) for gene in expressedGenes ]
				SUMMARY.write("\n")

				for geneSet in geneSets:
					SUMMARY.write("{0}".format(geneSet))
					for x,y in expressedGenes:
						if geneSet in ["neighborhoods"][analysis]:
							SUMMARY.write("\t1")
						else:
							SUMMARY.write("\t0")
					SUMMARY.write("\n")

			with open("{0}resultSummary_neighborhoods_patients_{1}.tsv".format(options.Options().qout,analysis), "w" ) as SUMMARY:
				SUMMARY.write("GeneSet")
				[ SUMMARY.write("\t{0}".format(patient)) for patient in patients ]
				SUMMARY.write("\n")

				for geneSet in geneSets:
					SUMMARY.write("{0}".format(geneSet))
					for patient in patients:
						for gene,info in self._gene_network.nodes(data=True):
							if patient in [ x.patients for x in info["isoformSwitches"] ]:
								SUMMARY.write("\t1")
							else:
								SUMMARY.write("\t0")
						SUMMARY.write("\n")