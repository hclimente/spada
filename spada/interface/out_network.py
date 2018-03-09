from spada import utils

import logging

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
