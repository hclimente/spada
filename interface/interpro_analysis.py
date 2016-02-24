from libs import options
from libs import utils

import os

class InterproAnalysis:
    def __init__(self):
        pass

    def launchAnalysis(self,tx,seq):
        out = "{}data/{}/InterPro/{}.tsv".format(options.Options().wd, options.Options().annotation, tx)

        if not os.path.isfile(out):
            proc = utils.cmdOut("{}pipeline/libs/interpro/interproscan.sh".format(options.Options().wd),
                                "--disable-precalc","--outfile",out,"--formats","TSV",
                                "--input","{0}protein.fa".format(options.Options().qout))
        return out

    def readInterpro(self,interproOut,protein):

        for cols in utils.readTable(interproOut, header=False):

            acceptedAnalysis = ["Pfam"]

            # https://code.google.com/p/interproscan/wiki/OutputFormats#Tab-separated_values_format_%28TSV%29
            protein_accession   = cols[0] #Protein Accession
            md5_digest          = cols[1] #Sequence MD5 digest
            seq_length          = cols[2] #Sequence Length
            analysis            = cols[3] #Analysis
            signature_accession = cols[4] #Signature Accession
            signature_descript  = cols[5] #Signature Description
            start               = int(cols[6]) #Start location
            stop                = int(cols[7]) #Stop location
            try:
                score           = float(cols[8]) #Score - is the e-value of the match reported by member database method
            except ValueError:
                score           = None
            status              = True if cols[9] == "T" else False #Status - is the status of the match (T: true)
            date                = cols[10] #Date - is the date of the run
            if len(cols) > 11:
                interpro_annotation = cols[11] #(InterPro annotations - accession)
            if len(cols) > 12:
                interpro_descript   = cols[12] #(InterPro annotations - description)
            if len(cols) > 13:
                go_annotation       = cols[13] #(GO annotations)
            if len(cols) > 14:
                pathway_annotation  = cols[14] #(Pathways annotations)

            tx=protein_accession.strip().split("#")[0]

            if score and score > 0.01: 
                continue
            elif analysis not in acceptedAnalysis: 
                continue
            elif tx != protein.tx:
                raise Exception("Error reading InterPro features for {0}. Invalid identifier {1} found.".format(protein.tx,tx))
            
            isoSpecificRes = set([ x._num for x in protein._structure if x.isoformSpecific ])
            featureRes = set(range(start,stop+1))
            intersection = float(len(isoSpecificRes & featureRes))

            featInfo = {}
            featInfo["region"]          = [start,stop]
            featInfo["accession"]       = signature_accession
            featInfo["description"]     = signature_descript
            featInfo["analysis"]        = analysis
            featInfo["percentAffected"] = intersection/len(featureRes) * 100

            yield featInfo
