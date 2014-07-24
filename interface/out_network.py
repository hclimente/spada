from biological_entities import transcript
from libs import options

import math
import networkx
import logging

def interactorsInfo(gn_net, tx_net, gene, nIso, tIso):
	geneSym = gn_net._net.node[gene]["symbol"]
	logging.info("Writing InteraX file for gene {0}({1}), isoforms {2} and {3}.".format(
						gene, geneSym, nIso, tIso ))
	with open("{0}/{1}/InteraXv2_{2}_{3}_{4}.tsv".format(options.Options().qout, options.Options().iLoopsVersion, geneSym, nIso, tIso), "w" ) as INTERAX:
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

			nTPM = None
			if tx_net._net.node[partnerIso]["median_TPM_N"] is not None:
				nTPM = math.log10 ( tx_net._net.node[partnerIso]["median_TPM_N"] + 0.0001 )
			tTPM = None
			if tx_net._net.node[partnerIso]["median_TPM_T"] is not None:
				tTPM = math.log10 ( tx_net._net.node[partnerIso]["median_TPM_T"] + 0.0001 )

			INTERAX.write("{0}\t".format( partnerIso ))
			INTERAX.write("{0}\t".format( partnerGene ))
			INTERAX.write("{0}\t".format( gn_net._net.node[ partnerGene ]["symbol"] ))
			INTERAX.write("{0}\t".format( RCN ))
			INTERAX.write("{0}\t".format( RCT ))
			INTERAX.write("{0}\t".format( dRC ))
			INTERAX.write("{0}\t".format( annotation ))
			INTERAX.write("{0}\t".format( nTPM ))
			INTERAX.write("{0}\t".format( tx_net._net.node[partnerIso]["median_PSI_N"] ))
			INTERAX.write("{0}\t".format( tTPM ))
			INTERAX.write("{0}\t".format( tx_net._net.node[partnerIso]["median_PSI_T"] ))
			INTERAX.write("\n")

def getGUILDInput(gn_net, onlyExperimental=False):
	logging.info("Writing GUILD input files.")
	with open(options.Options().qout + "GUILD/guild_nodes.tsv", "w" ) as GUILD_NODES:
		for node,info in gn_net.nodes(data=True):
			GUILD_NODES.write("{0} {1}\n".format(node, info["score"]))

	with open(options.Options().qout + "GUILD/guild_edges.tsv", "w" ) as GUILD_EDGES:
		for node1,node2,info in gn_net.edges(data=True):
			if onlyExperimental:
				if info["experimental"]:
					GUILD_EDGES.write("{0} {1} {2}\n".format(node1, info["score"], node2))
			else:
				GUILD_EDGES.write("{0} {1} {2}\n".format(node1, info["score"], node2))

def outputDot(network, name):
	networkx.write_dot(network, options.Options().qout + name)

def outputGTF(nodesOfInterest, tx_network):
	logging.info("Writing GTF files.")
	with open(options.Options().qout + "/candidates_normal.gtf", 'w') as nGTF, \
		 open(options.Options().qout + "/candidates_tumor.gtf", 'w') as tGTF, \
		 open("Data/" + options.Options().inputType + "/annotation.gtf", "r") as ALLTRANSCRIPTS:
	
		switchesInfo = []

		for gene, properties in nodesOfInterest:
			if not properties["isoformSwitches"]: continue

			for switch,score,patients in properties["isoformSwitches"]:
				nIso 	= switch[0]
				tIso 	= switch[1]

				switchesInfo.append([(nIso, tIso), score])

		for line in ALLTRANSCRIPTS:
			for switch in switchesInfo:
				if switch[0][0] in line:
					nGTF.write("{0};patients_affected={1}\n".format(line.strip(), switch[1] ))
				elif switch[0][1] in line:
					tGTF.write("{0};patients_affected={1}\n".format(line.strip(), switch[1] ))

def outCandidateList(gn_network, tx_network):
	logging.info("Writing candidateList_v2.")
	sortedNodes = sorted(gn_network.nodes(data=True), key=lambda (a, dct): dct['score'], reverse=True)
	with open(options.Options().qout + "candidateList_v2.tsv", "w") as cList:
		cList.write("GeneId\tSymbol\tNormal_transcript\tTumor_transcript\tNormal_protein\t")
		cList.write("Tumor_protein\tPatient_percentage\tDriver\tEpigenetic_factor\tRBP\t")
		cList.write("CDS\tCDS_change\tUTR_change\tPatients_affected\n")
		
		for gene, geneProperties in sortedNodes:
			if not geneProperties["isoformSwitches"]: continue

			for switch,score,patients in geneProperties["isoformSwitches"]:

				nIso = transcript.Transcript( switch[0], tx_network._net.node[switch[0]] )
				tIso = transcript.Transcript( switch[1], tx_network._net.node[switch[1]] )

				nUniprot 	= tx_network._net.node[switch[0]]["Uniprot"]
				tUniprot 	= tx_network._net.node[switch[1]]["Uniprot"]
				cds 		= False
				cdsChange 	= False
				utrChange 	= False

				if nIso.cds or tIso.cds: 	cds 		= True
				if nIso.get_cdsDiff(tIso): 	cdsChange 	= True
				if nIso.get_utrDiff(tIso): 	utrChange 	= True

				cList.write("{0}\t{1}\t".format( gene, geneProperties["symbol"] ))
				cList.write("{0}\t{1}\t".format( nIso, tIso ))
				cList.write("{0}\t{1}\t".format( nUniprot, tUniprot ))
				cList.write("{0}\t{1}\t".format( score, geneProperties["Driver"] ))
				cList.write("{0}\t{1}\t".format( geneProperties["EpiFactor"], geneProperties["RBP"] ))
				cList.write("{0}\t{1}\t".format( cds, cdsChange ))
				cList.write("{0}\t{1}\n".format( utrChange, ",".join(patients) ))