from libs import options

import math
import networkx
import logging

class OutNetwork:
	def __init__(self):
		pass

	def interactorsInfo(self, gn_net, tx_net, gene, nIso, tIso):
		geneSym = gn_net._net.node[gene]["symbol"]
		logging.info("Writing InteraX file for gene {0}({1}), isoforms {2} and {3}.".format(
							gene, geneSym, nIso, tIso ))
		with open("{0}InteraXv2_{1}_{2}_{3}.tsv".format(options.Options().qout, geneSym, nIso, tIso), "w" ) as INTERAX:
			INTERAX.write("#")
			INTERAX.write("\t{0}\t{1}".format( nIso, tx_net._net.node[nIso]["median_PSI_N"] ))
			INTERAX.write("\t{0}\t{1}".format( tIso, tx_net._net.node[tIso]["median_PSI_T"] ))
			INTERAX.write("\t{0}".format( gene ))
			INTERAX.write("\t{0}".format( geneSym ))
			INTERAX.write("\t{0}".format( math.log10 ( tx_net._net.node[nIso]["median_TPM_N"] + 0.0001 ) ))
			INTERAX.write("\t{0}".format( math.log10 ( tx_net._net.node[tIso]["median_TPM_N"] + 0.0001 ) ))
			INTERAX.write("\n")

			for partnerIso in set( tx_net._net.neighbors(nIso) ) | set( tx_net._net.neighbors(tIso) ):
				partnerGene = tx_net._net.node[partnerIso]["gene_id"]
				annotation = "Driver" if gn_net._net.node[partnerGene]["Driver"] else ""
				RCN = "" if partnerIso not in tx_net._net.neighbors(nIso) else tx_net._net.edge[nIso][partnerIso]["RC"]
				RCT = "" if partnerIso not in tx_net._net.neighbors(tIso) else tx_net._net.edge[tIso][partnerIso]["RC"]
				dRC = int(RCN - RCT) if RCN and RCT else ""

				INTERAX.write("{0}\t".format( partnerIso ))
				INTERAX.write("{0}\t".format( partnerGene ))
				INTERAX.write("{0}\t".format( gn_net._net.node[ partnerGene ]["symbol"] ))
				INTERAX.write("{0}\t".format( RCN ))
				INTERAX.write("{0}\t".format( RCT ))
				INTERAX.write("{0}\t".format( dRC ))
				INTERAX.write("{0}\t".format( annotation ))
				INTERAX.write("{0}\t".format( math.log10 ( tx_net._net.node[partnerIso]["median_TPM_N"] + 0.0001 ) ))
				INTERAX.write("{0}\t".format( tx_net._net.node[partnerIso]["median_PSI_N"] ))
				INTERAX.write("{0}\t".format( math.log10 ( tx_net._net.node[partnerIso]["median_TPM_T"] + 0.0001 ) ))
				INTERAX.write("{0}\t".format( tx_net._net.node[partnerIso]["median_PSI_T"] ))
				INTERAX.write("\n")

	def getGUILDInput(self, gn_net):
		with open(options.Options().qout + "guild_nodes.tsv" ) as GUILD_NODES:
			for node, info in gn_net.nodes(data=True):
				GUILD_NODES.write("{0}\t{1}\n".format(node, info["score"]))

		with open(options.Options().qout + "guild_edges.tsv" ) as GUILD_EDGES:
			for node1, node2, info in gn_net.edges(data=True):
				GUILD_EDGES.write("{0}\t{1}\t{2}\n".format(node1, info["score"], node2))

	def outputDot(self, network, name):
		networkx.write_dot(network, options.Options().qout + name)

	def outputGTF(self, nodesOfInterest, tx_network):
		with open(options.Options().qout + "/candidates_normal.gtf", 'w') as nGTF, \
			 open(options.Options().qout + "/candidates_tumor.gtf", 'w') as tGTF, \
			 open("Data/" + options.Options().inputType + "/annotation.gtf", "r") as ALLTRANSCRIPTS:
		
			switchesInfo = []

			for gene, properties in sortedNodes:
				if not properties["isoformSwitches"]: continue

				for switch, score in properties["isoformSwitches"]:
					nIso 	= switch[0]
					tIso 	= switch[1]

					switchesInfo.append((nIso, tIso), score)

			for line in ALLTRANSCRIPTS:
				for switch in switchesInfo:
					if switch[0][0] in line:
						nGTF.write("{0};patients_affected={1}\n".format(line.strip(), switch[1] ))
					elif switch[0][1] in line:
						tGTF.write("{0};patients_affected={1}\n".format(line.strip(), switch[1] ))

	def outCandidateList(self, gn_network, tx_network):
		sortedNodes = sorted(gn_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
		with open(options.Options().qout + "candidateList_v2.tsv") as cList:
			cList.write("GeneId\tSymbol\tNormal_transcript\tTumor_transcript\t")
			cList.write("Affected_patients\tDriver\tEpigenetic_factor\tRBP\t")
			cList.write("CDS\tCDS_change\tUTR_change\n")
			for gene, properties in sortedNodes:
				if not properties["isoformSwitches"]: continue

				symbol 	= properties["symbol"]
				driver 	= properties["Driver"]
				epiFac 	= properties["EpiFactor"]
				rbp 	= properties["RBP"]

				for switch, score in properties["isoformSwitches"]:
					nIso 	= switch[0]
					tIso 	= switch[1]
					nInfo 	= tx_network._net.node[nIso]
					tInfo 	= tx_network._net.node[tIso]

					cds 		= False
					cdsChange 	= False
					utrChange 	= False

					if nInfo["cdsCoords"][0] != nInfo["cdsCoords"][1] or tInfo["cdsCoords"][0] != tInfo["cdsCoords"][1]:
						cds = True
					if nInfo["txCoords"][0] != tInfo["txCoords"][0] or nInfo["txCoords"][1] != tInfo["txCoords"][1]:
						utrChange = True
					elif nInfo["cdsCoords"][0] != tInfo["cdsCoords"][0] or nInfo["cdsCoords"][1] != tInfo["cdsCoords"][1]:
						utrChange = True
					if nInfo["exonStructure"] != tInfo["exonStructure"]:
						cdsChange = True

					cList.write("{0}\t{1}\t".format( gene, symbol ))
					cList.write("{0}\t{1}\t".format( nIso, tIso ))
					cList.write("{0}\t{1}\t".format( score, driver ))
					cList.write("{0}\t{1}\t".format( epiFac, rbp ))
					cList.write("{0}\t{1}\t".format( cds, cdsChange ))
					cList.write("{0}\n".format( utrChange ))

# #CDS and UTR change
# if txInfo[candidate1]["cdsStart"] == txInfo[candidate1]["cdsEnd"] or txInfo[candidate2]["cdsStart"] == txInfo[candidate2]["cdsEnd"]:
# 	aCandidate["CDS?"] = "No"

# if txInfo[candidate1]["exonStarts"] != txInfo[candidate2]["exonStarts"] or txInfo[candidate1]["exonEnds"] != txInfo[candidate2]["exonEnds"]:
# 	if txInfo[candidate1]["cdsStart"] != txInfo[candidate2]["cdsStart"] or txInfo[candidate1]["cdsEnd"] != txInfo[candidate2]["cdsEnd"]:
# 		aCandidate["CDS Change"] = "Yes"
# 	if txInfo[candidate1]["txStart"] != txInfo[candidate2]["txStart"] or txInfo[candidate1]["txEnd"] != txInfo[candidate2]["txEnd"]:
# 		aCandidate["UTR Change"] = "Yes"

# #Check if each of the exons at the CDS are the same.
# for anExonStart_1, anExonEnd_1 in zip(txInfo[candidate1]["exonStarts"], txInfo[candidate1]["exonEnds"] ):
# 	exonMatch = False
	
# 	for anExonStart_2, anExonEnd_2 in zip(txInfo[candidate2]["exonStarts"], txInfo[candidate2]["exonEnds"] ):
# 		if anExonStart_1 == anExonStart_2 and anExonEnd_1 == anExonEnd_2:
# 			exonMatch = True
# 			break
	
# 	if not exonMatch:
# 		if anExonStart_1 < txInfo[candidate1]["cdsStart"] or anExonStart_1 > txInfo[candidate1]["cdsEnd"]:
# 			aCandidate["UTR Change"] = "Yes"
# 		else:
# 			aCandidate["CDS Change"] = "Yes"

# Uniprot_N = "None"
# Uniprot_T = "None"
# emptySeqRecord = SeqRecord(Seq("", IUPAC.protein), id="", name="", description="")

# for uniprot_iso, uniprot_seq in Uniprot_dict.iteritems():
# 	if str(Protein_dict.get(candidate1, emptySeqRecord).seq) == str(uniprot_seq.seq):
# 		if Uniprot_N == "None":
# 			Uniprot_N = uniprot_iso.split("|")[1]
# 		else:
# 			Uniprot_N += ";" + uniprot_iso.split("|")[1]
# 	if str(Protein_dict.get(candidate2, emptySeqRecord).seq) == str(uniprot_seq.seq):
# 		if Uniprot_T == "None":
# 			Uniprot_T = uniprot_iso.split("|")[1]
# 		else:
# 			Uniprot_T += ";" + uniprot_iso.split("|")[1]