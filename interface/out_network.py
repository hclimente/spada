from libs import options

import math
import networkx
import logging

def interactorsInfo(gn_net, tx_net, gene, nIso, tIso):
	geneSym = gn_net._net.node[gene]["symbol"]
	logging.info("Writing InteraX file for gene {0}({1}), isoforms {2} and {3}.".format(
						gene, geneSym, nIso, tIso ))
	with open("{0}iLoops/{1}/InteraXv2_{2}_{3}_{4}.tsv".format(options.Options().qout, options.Options().iLoopsVersion, geneSym, nIso, tIso), "w" ) as INTERAX:
		
		nTPM_N = None
		if tx_net._net.node[nIso]["median_TPM_N"] is not None:
			nTPM_N = math.log10 ( tx_net._net.node[nIso]["median_TPM_N"] + 0.0001 )

		tTPM_T = None
		if tx_net._net.node[nIso]["median_TPM_N"] is not None:
			tTPM_T = math.log10 ( tx_net._net.node[tIso]["median_TPM_T"] + 0.0001 )

		INTERAX.write("#\t")
		INTERAX.write("{0}\t{1}\t".format( nIso, tx_net._net.node[nIso]["median_PSI_N"] ))
		INTERAX.write("{0}\t{1}\t".format( tIso, tx_net._net.node[tIso]["median_PSI_T"] ))
		INTERAX.write("{0}\t{1}\t".format( gene, geneSym ))
		INTERAX.write("{0}\t{1}\n".format( nTPM_N, tTPM_T ))

		INTERAX.write("Partner\tGene\tSymbol\tRC_n\tRC_t\tdeltaRC\t")
		INTERAX.write("Annotation\tTPM_n\tPSI_n\tTPM_t\tPSI_t\n")

		for partnerIso in set( tx_net._net.neighbors(nIso) ) | set( tx_net._net.neighbors(tIso) ):
			partnerGene = tx_net._net.node[partnerIso]["gene_id"]
			annotation = "Driver" if gn_net._net.node[partnerGene]["Driver"] else ""
			RCN = None if partnerIso not in tx_net._net.neighbors(nIso) else tx_net._net.edge[nIso][partnerIso]["RC"]
			RCT = None if partnerIso not in tx_net._net.neighbors(tIso) else tx_net._net.edge[tIso][partnerIso]["RC"]
			dRC = int(RCN - RCT) if RCN and RCT else None

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
	outputFolder = "{0}GUILD_enriched/".format(options.Options().qout)
	if onlyExperimental: outputFolder = "{0}GUILD_experimental/".format(options.Options().qout)

	with open("{0}guild_nodes.tsv".format(outputFolder), "w" ) as GUILD_NODES:
		for node,info in gn_net.nodes(data=True):
			GUILD_NODES.write("{0} {1}\n".format(node, info["score"]))

	with open("{0}guild_edges.tsv".format(outputFolder), "w" ) as GUILD_EDGES:
		for node1,node2,info in gn_net.edges(data=True):
			if onlyExperimental:
				if info["experimental"]:
					GUILD_EDGES.write("{0} {1} {2}\n".format(node1, info["score"], node2))
			else:
				GUILD_EDGES.write("{0} {1} {2}\n".format(node1, info["score"], node2))

def outputGTF(gn_network, tx_network):
	logging.info("Writing GTF files.")
	with open(options.Options().qout + "/candidates_normal.gtf", 'w') as nGTF, \
		 open(options.Options().qout + "/candidates_tumor.gtf", 'w') as tGTF, \
		 open("Data/" + options.Options().inputType + "/annotation.gtf", "r") as ALLTRANSCRIPTS:
	
		switchesInfo = [ [(z.nTx,z.tTx),z.score] for x,y,z in utils.iterate_switches_ScoreWise(gn_network) ]

		for line in ALLTRANSCRIPTS:
			for switch in switchesInfo:
				if switch[0][0] in line:
					nGTF.write("{0};patients_affected={1}\n".format(line.strip(), switch[1] ))
				elif switch[0][1] in line:
					tGTF.write("{0};patients_affected={1}\n".format(line.strip(), switch[1] ))

def outCandidateList(gn_network, tx_network):
	logging.info("Writing candidateList_v2.")
	with open(options.Options().qout + "candidateList_v2.tsv", "w") as cList:
		cList.write("GeneId\tSymbol\tNormal_transcript\tTumor_transcript\tNormal_protein\t")
		cList.write("Tumor_protein\tPatient_percentage\tDriver\tEpigenetic_factor\tRBP\t")
		cList.write("CDS\tCDS_change\tUTR_change\tPatients_affected\n")
		
		for gene,info,switch in utils.iterate_switches_ScoreWise(gn_network):
			nIso = switch.nTranscript
			tIso = switch.tTranscript

			nUniprot 	= tx_network._net.node[switch.nTx]["Uniprot"]
			tUniprot 	= tx_network._net.node[switch.tTx]["Uniprot"]
			cds 		= False
			cdsChange 	= False
			utrChange 	= False

			if nIso.cds or tIso.cds: 	cds 		= True
			if switch.cds_diff: 	 	cdsChange 	= True
			if switch.utr_diff: 		utrChange 	= True

			cList.write("{0}\t{1}\t".format( gene, info["symbol"] ))
			cList.write("{0}\t{1}\t".format( nIso.name, tIso.name ))
			cList.write("{0}\t{1}\t".format( nUniprot, tUniprot ))
			cList.write("{0}\t{1}\t".format( switch.score, info["Driver"] ))
			cList.write("{0}\t{1}\t".format( info["EpiFactor"], info["RBP"] ))
			cList.write("{0}\t{1}\t".format( cds, cdsChange ))
			cList.write("{0}\t{1}\n".format( utrChange, ",".join(switch.patients) ))

def outTSV(network,path):
	with open(path + "_nodes.tsv", "w") as NODES:
		columns = network.nodes(data=True)[0][1].keys()
		NODES.write("Node")
		for col in columns:
			NODES.write("\t"+col)
		NODES.write("\n")

		for node,properties in network.nodes(data=True):
			NODES.write(node)
			for key in properties:
				value = ""
				if isinstance(properties[key],set) or isinstance(properties[key],list):
					for thing in properties[key]:
						if isinstance(thing,list) or isinstance(thing,set):
							value += "#" + ";".join(thing)
						else:
							value += str(thing) + ";"
				else:
					value = str(properties[key])
				NODES.write("\t"+value)
			NODES.write("\n")

	with open(path + "_edges.tsv", "w") as EDGES:
		columns = network.edges(data=True)[0][2].keys()
		EDGES.write("Node_1\tNode_2")
		for col in columns:
			EDGES.write("\t"+col)
		EDGES.write("\n")
		
		for node1,node2,properties in network.edges(data=True):
			EDGES.write(node1 + "\t" + node2)
			for key in properties:
				EDGES.write("\t"+str(properties[key]))
			EDGES.write("\n")