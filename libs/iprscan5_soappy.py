#!/usr/bin/env python
# $Id: iprscan5_soappy.py 2760 2014-04-10 15:24:31Z hpm $
# ======================================================================
# 
# Copyright 2008-2014 EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# ======================================================================
# InterProScan 5 (SOAP) service, Python client using SOAPpy.
#
# Tested with:
#   Python 2.5.2 with SOAPpy 0.12.0 (Ubuntu 8.04 LTS)
#   Python 2.6.5 with SOAPpy 0.12.0 (Ubuntu 10.04 LTS)
#   Python 2.7.3 with SOAPpy 0.12.0 (Ubuntu 12.04 LTS)
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_soap
# http://www.ebi.ac.uk/Tools/webservices/tutorials/06_programming/python
# ======================================================================
# WSDL URL for service
wsdlUrl = 'http://www.ebi.ac.uk/Tools/services/soap/iprscan5?wsdl'

# Load libraries
import base64, platform, os, SOAPpy, sys, time
import warnings
from SOAPpy import WSDL

# Suppress all deprecation warnings (not recommended for development)
warnings.simplefilter('ignore', DeprecationWarning)

# Set interval for checking status
checkInterval = 3
# Output level
outputLevel = 1
# Debug level
debugLevel = 0

# Tool specific options
appl = #', help='signature methods to use, see --paramDetail appl')
crc = #', action="store_true", help='enable InterProScan Matches look-up (ignored)')
nocrc = #', action="store_true", help='disable InterProScan Matches look-up (ignored)')
goterms = #', action="store_true", help='enable inclusion of GO terms')
nogoterms = #', action="store_true", help='disable inclusion of GO terms')
pathways = #', action="store_true", help='enable inclusion of pathway terms')
nopathways = #', action="store_true", help='disable inclusion of pathway terms')
sequence = #', help='input sequence file name')
# General options
email = #', help='e-mail address')
title = #', help='job title')
outfile = #', help='file name for results')
outformat = #', help='output format for results')
async = #', action='store_true', help='asynchronous mode')
jobid = #', help='job identifier')
polljob = #', action="store_true", help='get job result')
status = #', action="store_true", help='get job status')
resultTypes = #', action='store_true', help='get result types')
params = #', action='store_true', help='list input parameters')
paramDetail = #', help='get details for parameter')
quiet = #', action='store_true', help='decrease output level')
verbose = #', action='store_true', help='increase output level')
trace = #', action="store_true", help='show SOAP messages')
WSDL = #', default=wsdlUrl, help='WSDL URL for service')
debugLevel = #', type='int', default=debugLevel, help='debug output level')

# Increase output level
if options.verbose:
    outputLevel += 1

# Decrease output level
if options.quiet:
    outputLevel -= 1

# Debug level
if options.debugLevel:
    debugLevel = options.debugLevel

# Debug print
def printDebugMessage(functionName, message, level):
    if(level <= debugLevel):
        print >>sys.stderr, '[' + functionName + '] ' + message

# Get input parameters list
def serviceGetParameters():
    printDebugMessage('serviceGetParameters', 'Begin', 1)
    result = server.getParameters()
    printDebugMessage('serviceGetParameters', 'End', 1)
    return result

# Get input parameter information
def serviceGetParameterDetails(paramName):
    printDebugMessage('serviceGetParameterDetails', 'Begin', 1)
    result= server.getParameterDetails(parameterId=paramName)
    printDebugMessage('serviceGetParameterDetails', 'End', 1)
    return result

# Submit job
def serviceRun(email, title, params):
    printDebugMessage('serviceRun', 'Begin', 1)
    jobid = server.run(email=email, title=title, parameters=params)
    printDebugMessage('serviceRun', 'End', 1)
    return jobid

# Get job status
def serviceCheckStatus(jobId):
    printDebugMessage('serviceCheckStatus', 'jobId: ' + jobId, 1)
    result = server.getStatus(jobId = jobId)
    return result

# Get available result types for job
def serviceGetResultTypes(jobId):
    printDebugMessage('serviceGetResultTypes', 'Begin', 1)
    result = server.getResultTypes(jobId=jobId)
    printDebugMessage('serviceGetResultTypes', 'End', 1)
    return result['type']

# Get result
def serviceGetResult(jobId, type):
    printDebugMessage('serviceGetResult', 'Begin', 1)
    printDebugMessage('serviceGetResult', 'jobId: ' + jobId, 1)
    printDebugMessage('serviceGetResult', 'type: ' + type, 1)
    resultBase64 = server.getResult(jobId=jobId, type=type)
    result = base64.decodestring(resultBase64)
    printDebugMessage('serviceGetResult', 'End', 1)
    return result

# Client-side poll
def clientPoll(jobId):
    printDebugMessage('clientPoll', 'Begin', 1)
    result = 'PENDING'
    while result == 'RUNNING' or result == 'PENDING':
        result = serviceCheckStatus(jobId)
        print >>sys.stderr, result
        if result == 'RUNNING' or result == 'PENDING':
            time.sleep(15)
    printDebugMessage('clientPoll', 'End', 1)

# Get result for a jobid
def getResult(jobId):
    printDebugMessage('getResult', 'Begin', 1)
    printDebugMessage('getResult', 'jobId: ' + jobId, 1)
    # Check status and wait if necessary
    clientPoll(jobId)
    # Get available result types
    resultTypes = serviceGetResultTypes(jobId)
    for resultType in resultTypes:
        # Get the result
        result = serviceGetResult(jobId, resultType['identifier'])
        # Derive the filename for the result
        if options.outfile:
            filename = options.outfile + '.' + resultType['identifier'] + '.' + resultType['fileSuffix']
        else:
            filename = jobId + '.' + resultType['identifier'] + '.' + resultType['fileSuffix']
        # Write a result file
        if not options.outformat or options.outformat == resultType['identifier']:
            fh = open(filename, 'w');
            fh.write(result)
            fh.close()
            print filename
    printDebugMessage('getResult', 'End', 1)

# Read a file
def readFile(filename):
    printDebugMessage('readFile', 'Begin', 1)
    fh = open(filename, 'r')
    data = fh.read()
    fh.close()
    printDebugMessage('readFile', 'End', 1)
    return data

# Output parameter details.
def printGetParameterDetails(paramName):
    printDebugMessage('printGetParameterDetails', 'Begin', 1)
    paramDetail = serviceGetParameterDetails(paramName)
    print paramDetail['name'], "\t", paramDetail['type']
    print paramDetail['description']
    for value in paramDetail['values']['value']:
        print value['value'],
        if(value['defaultValue'] == 'true'):
            print '(default)',
        print
        print "\t", value['label']
        if(hasattr(value, 'properties')):
            if(isinstance(value['properties']['property'], (list, tuple))):
                for wsProperty in value['properties']['property']:
                    print "\t", wsProperty['key'], "\t", wsProperty['value']
            else:
                print "\t", value['properties']['property']['key'], "\t", value['properties']['property']['value']
    printDebugMessage('printGetParameterDetails', 'End', 1)

# Output available result types for job.
def printGetResultTypes(jobId):
    printDebugMessage('printGetResultTypes', 'Begin', 1)
    for resultType in serviceGetResultTypes(jobId):
        print resultType['identifier']
        if(hasattr(resultType, 'label')):
            print "\t", resultType['label']
        if(hasattr(resultType, 'description')):
            print "\t", resultType['description']
        if(hasattr(resultType, 'mediaType')):
            print "\t", resultType['mediaType']
        if(hasattr(resultType, 'fileSuffix')):
            print "\t", resultType['fileSuffix']
    printDebugMessage('printGetResultTypes', 'End', 1)

# Set the client user-agent.
clientRevision = '$Revision: 2760 $'
clientVersion = '0'
if len(clientRevision) > 11:
    clientVersion = clientRevision[11:-2] 
userAgent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
    clientVersion, os.path.basename( __file__ ),
    platform.python_version(), platform.system(),
    SOAPpy.Client.SOAPUserAgent()
)
# Function to return User-agent.
def SOAPUserAgent():
    return userAgent
# Redefine default User-agent function to return custom User-agent.
SOAPpy.Client.SOAPUserAgent = SOAPUserAgent
printDebugMessage('main', 'User-agent: ' + SOAPpy.Client.SOAPUserAgent(), 1)

# Create the service interface
printDebugMessage('main', 'WSDL: ' + options.WSDL, 1)
server = WSDL.Proxy(options.WSDL)

# Fix message namespace (not set from the WSDL).
for method in server.methods:
    if server.methods[method].namespace == None:
        server.methods[method].namespace = 'http://soap.jdispatcher.ebi.ac.uk'

# Configure HTTP proxy from OS environment (e.g. http_proxy="http://proxy.example.com:8080")
if os.environ.has_key('http_proxy'):
    http_proxy_conf = os.environ['http_proxy'].replace('http://', '')
elif os.environ.has_key('HTTP_PROXY'):
    http_proxy_conf = os.environ['HTTP_PROXY'].replace('http://', '')
else:
    http_proxy_conf = None
server.soapproxy.http_proxy = http_proxy_conf

# If required enable SOAP message trace
if options.trace:
    server.soapproxy.config.dumpSOAPOut = 1
    server.soapproxy.config.dumpSOAPIn = 1

# No options... print help.
if numOpts < 2:
    parser.print_help()
# List parameters
elif options.params:
    for paramName in serviceGetParameters()['id']:
        print paramName
# Get parameter details
elif options.paramDetail:
    printGetParameterDetails(options.paramDetail)
# Submit job
elif options.email and not options.jobid:
    params = {}
    if len(args) > 0:
        if os.access(args[0], os.R_OK): # Read file into content
            params['sequence'] = readFile(args[0])
        else: # Argument is a sequence id
            params['sequence'] = args[0]
    elif options.sequence: # Specified via option
        if os.access(options.sequence, os.R_OK): # Read file into content
            params['sequence'] = readFile(options.sequence)
        else: # Argument is a sequence id
            params['sequence'] = options.sequence
    # Booleans need to be represented as 1/0 rather than True/False
    #if options.crc:
    #    params['nocrc'] = 0
    #if options.nocrc:
    #    params['nocrc'] = 1
    if options.goterms:
        params['goterms'] = 1
    if options.nogoterms:
        params['goterms'] = 0
    if options.pathways:
        params['pathways'] = 1
    if options.nopathways:
        params['pathways'] = 0
    # Add the other options (if defined)
    if options.appl:
        params['appl'] = {'string':options.appl}
    # Submit the job
    jobid = serviceRun(options.email, options.title, params)
    if options.async: # Async mode
        print jobid
    else: # Sync mode
        if outputLevel > 0:
            print jobid
        time.sleep(5)
        getResult(jobid)
# Get job status
elif options.status and options.jobid:
    status = serviceCheckStatus(options.jobid)
    print status
# List result types for job
elif options.resultTypes and options.jobid:
    printGetResultTypes(options.jobid)
# Get results for job
elif options.polljob and options.jobid:
    getResult(options.jobid)
else:
    print 'Error: unrecognised argument combination'
    parser.print_help()
