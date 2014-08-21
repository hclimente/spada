from libs import options

import base64, os, SOAPpy, time

import pdb

class InterproAnalysis:
    def __init__(self):
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
        out = "{0}InterPro/{1}/{2}.tsv.txt".format(options.Options().wd, options.Options().inputType, tx)

        if os.path.isfile(out):
            return out

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