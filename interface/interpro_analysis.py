from libs import options
from libs import utils

import base64, os, SOAPpy, time

class InterproAnalysis:
    def __init__(self,mode="offline"):
        self._mode = mode
        if self._mode=="online":

            SOAPpy.Client.SOAPUserAgent = self.SOAPUserAgent

            # Create the service interface
            self._server = SOAPpy.WSDL.Proxy('http://www.ebi.ac.uk/Tools/services/soap/iprscan5?wsdl')

            # Fix message namespace (not set from the WSDL).
            for method in self._server.methods:
                if self._server.methods[method].namespace == None:
                    self._server.methods[method].namespace = 'http://soap.jdispatcher.ebi.ac.uk'

            # Configure HTTP proxy from OS environment (e.g. http_proxy="http://proxy.example.com:8080")
            if os.environ.has_key('http_proxy'):
                http_proxy_conf = os.environ['http_proxy'].replace('http://', '')
            elif os.environ.has_key('HTTP_PROXY'):
                http_proxy_conf = os.environ['HTTP_PROXY'].replace('http://', '')
            else:
                http_proxy_conf = None
            self._server.soapproxy.http_proxy = http_proxy_conf

    def launchAnalysis(self,tx,seq):
        out = "{0}Data/{1}/InterPro/{2}.tsv".format(options.Options().wd, options.Options().annotation, tx)

        # asume that if the file is not present, 
        # its because no feature could be mapped
        if os.path.isfile(out):
            return out
        else:
            return None

        if self._mode=="online":
            params = {}
            params['sequence']  = seq
            params['goterms']   = 1
            params['pathways']  = 1

            jobId = self._server.run(email="hector.climente@upf.edu", title=tx, parameters=params)
            self.checkStatus(jobId)

            for resultType in self._server.getResultTypes(jobId=jobId)["type"]:
                if resultType['identifier'] != "tsv": continue
                
                resultBase64 = self._server.getResult(jobId=jobId, type=resultType['identifier'])
                result = base64.decodestring(resultBase64)

                with open(out, 'w') as interPro_Out:
                    interPro_Out.write(result)
        else:

            proc = utils.cmdOut("{0}Pipeline/libs/interpro/interproscan.sh".format(options.Options().wd),
                                "--disable-precalc","--outfile",out,"--formats","TSV",
                                "--input","{0}protein.fa".format(options.Options().qout))
        return out

    def SOAPUserAgent(self):
        return "SmartAS"

    def checkStatus(self, jobId):
        while True:
            result = self._server.getStatus(jobId = jobId)
            if result == 'RUNNING' or result == 'PENDING':
                time.sleep(15)
            else:
                return True

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
