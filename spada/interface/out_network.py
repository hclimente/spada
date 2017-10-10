from spada import utils

import logging

def outputGTF(gn_network, tx_network, annotation):
	logging.info("Writing GTF files.")
	with open("candidates_normal.gtf", 'w') as nGTF, \
		 open("candidates_tumor.gtf", 'w') as tGTF, \
		 open("data/{}/annotation.gtf".format(annotation), "r") as ALLTRANSCRIPTS:

		switchesInfo = [ (z.nTx,z.tTx) for w,x,y,z in gn_network.iterate_switches_byPatientNumber(tx_network,partialCreation=True) ]

		for line in ALLTRANSCRIPTS:
			for switch in switchesInfo:
				if switch[0] in line:
					nGTF.write(line)
				elif switch[1] in line:
					tGTF.write(line)

def outCandidateList(gn_network,tx_network,filename="candidateList_info.tsv"):
	logging.info("Writing candidateList.")
	with open(filename, "w") as cList:
		cList.write("GeneId\tSymbol\tNormal_transcript\tTumor_transcript\t")
		cList.write("Normal_protein\tTumor_protein\tAnnotation\tDriverAnnotation\t")
		cList.write("NotNoise\tIsModel\tIsFunctional\tDriver\tSpecificDriver\tDruggable\t")
		cList.write("CDS_Normal\tCDS_Tumor\tCDS_change\tUTR_change\tPatients_affected\n")

		hallmarksDict = utils.readGeneset("h.all.v5.0.entrez.gmt")
		bpDict = utils.readGeneset("c5.bp.v4.0.entrez.gmt")

		for gene,info,switchDict,switch in gn_network.iterate_switches_byPatientNumber(tx_network,partialCreation=True,removeNoise=False):
			nIso = switch.nTranscript
			tIso = switch.tTranscript

			nUniprot 	= tx_network._net.node[switch.nTx]["Uniprot"]
			tUniprot 	= tx_network._net.node[switch.tTx]["Uniprot"]
			cdsChange 	= False
			utrChange 	= False

			if switch.cds_diff: 	 	cdsChange 	= True
			if switch.utr_diff: 		utrChange 	= True

			try:
				relevance = int(switch.is_functional)
			except Exception:
				relevance = None

			annotation,driverAnnotation = gn_network.getGeneAnnotation(gene,hallmarksDict,bpDict)

			cList.write("{}\t{}\t".format( gene, info["symbol"] ))
			cList.write("{}\t{}\t".format( nIso.name, tIso.name ))
			cList.write("{}\t{}\t".format( nUniprot, tUniprot ))
			cList.write("{}\t{}\t".format( annotation,driverAnnotation ))
			cList.write("{}\t{}\t".format( int(not switchDict["noise"]), int(switchDict["model"]) ))
			cList.write("{}\t{}\t".format( relevance, int(info["driver"]) ))
			cList.write("{}\t".format( int(info["specificDriver"]) ))
			cList.write("{}\t{}\t".format( int(info["druggable"]), int(bool(nIso.cds)) ))
			cList.write("{}\t{}\t".format( int(bool(tIso.cds)), int(cdsChange), ))
			cList.write("{}\t{}\n".format( int(utrChange), ",".join(switch.patients) ))

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
