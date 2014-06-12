#!/soft/devel/python-2.7/bin/python

import pandas as pd
import argparse
from libs import utils as ut

class Options(object):
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Options, cls).__new__(cls, *args, **kwargs)
        return cls._instance

    def __init__(self):
        options = self.parseOptions()

        self._config_file           = options.config_file
        self._initial_step          = options.initial_step
        self._wd                    = options.wd
        self._gwd                   = options.gwd
        self._minExpression         = options.minExpression
        self._inputType             = options.inputType
        self._replicates            = options.replicates
        self._iLoopsVersion         = options.iLoopsVersion
        self._unpairedReplicates    = options.unpairedReplicates
        self._tag                   = options.tag
        self._external              = options.external
        if self._external:
            if self._unpairedReplicates == 0:
                self._out           = self._inputType + "/" + self._tag + "/"
        else:
            self._out               = self._inputType + "/" + self._tag + "_mE" + str(self._minExpression) + "/"
        self._gOut                  = self._gwd + "/Results/" + self._out
        self._quickOut              = self._wd + "/Results/" + self._out

    def configFile(self):           return self._config_file
    def initialStep(self):          return self._initial_step
    def wd(self):                   return self._wd
    def gwd(self):                  return self._gwd
    def minExpression(self):        return self._minExpression
    def inputType(self):            return self._inputType
    def replicates(self):           return self._replicates
    def iLoopsVersion(self):        return self._iLoopsVersion
    def unpairedReplicates(self):   return self._unpairedReplicates
    def tag(self):                  return self._tag
    def external(self):             return self._external
    def out(self):                  return self._out
    def gout(self):                 return self._gOut
    def qout(self):                 return self._quickOut

    def parseOptions(self, *args, **kwds):
        parser = Parser(prog="SmartAS.py", description = "Find significant alternative splicing switches. Analyze their functional impact.", epilog= "Hector Climente, 2014", fromfile_prefix_chars='@')

        parser.add_argument('-f','--config-file', dest='config_file', action = 'store', 
                            default='Parameters.cfg', help = 'Input file with the parameters')
        parser.add_argument('-s', '--initial-step', dest='initial_step', action='store', default='0',
                            type=int, help='Where should SmartAS start.')
        parser.add_argument('-wd', '--working-directory', dest='wd', action='store', default='/home/hector/SmartAS/',
                            help='Root file of SmartAS folder in the current machine.')
        parser.add_argument('-gwd', '--gaudi-wd', dest='gwd', action='store', 
                            default='/sbi/users/hectorc/SmartAS', help='Root file of SmartAS folder in Gaudi.')
        parser.add_argument('-m', '--minimum-expression', dest='minExpression', action='store', default='0',
                            type=float, help='Minimum expression to consider a transcript not residual.')
        parser.add_argument('-i', '--input-type', dest='inputType', action='store', default='TCGA',
                            help='Origin of the data.')
        parser.add_argument('-r', '--replicates', dest='replicates', action='store', default='0',
                            type=int, help='Number of patients or biological replicates.')
        parser.add_argument('-v', '--iloops-version', dest='iLoopsVersion', action='store', 
                            default='iLoops13', help='Version of iLoops to run.')
        parser.add_argument('-u', '--unpaired-replicates', dest='unpairedReplicates', action='store', 
                            default='0', type=int, help='Number of unpaired samples.')
        parser.add_argument('-t', '--tag', dest='tag', action='store', default='20',
                            help='Tag of the files.')
        parser.add_argument('-e', '--external', dest='external', action='store', default='',
                            help='Path of external data.')
        parser.add_argument('-o', '--output', dest='out', action='store', default='',
                            help='Path of output data, under the Results/ directory.')
        parser.add_argument('-go', '--goutput', dest='gout', action='store', default='',
                            help='Path of output data in Gaudi.')

        #"Conditions" : ["N", "T"]

        config_opt  = parser.parse_args()
        options     = parser.parse_args(["@" + config_opt.config_file])
        return options

    def printToFile(self):
        with open("Parameters.cfg", "w") as CONFIG:
            if self._initial_step:          CONFIG.write("initial-step=" + str(self._initial_step) + "\n")
            if self._wd:                    CONFIG.write("working-directory=" + self._wd + "\n")
            if self._gwd:                   CONFIG.write("gaudi-wd=" + self._gwd + "\n")
            if self._minExpression:         CONFIG.write("minimum-expression=" + str(self._minExpression) + "\n")
            if self._inputType:             CONFIG.write("input-type=" + self._inputType + "\n")
            if self._replicates:            CONFIG.write("replicates=" + str(self._replicates) + "\n")
            if self._iLoopsVersion:         CONFIG.write("iloops-version=" + self._iLoopsVersion + "\n")
            if self._unpairedReplicates:    CONFIG.write("unpaired-replicates=" + str(self._unpairedReplicates) + "\n")
            if self._tag:                   CONFIG.write("tag=" + self._tag + "\n")
            if self._external:              CONFIG.write("external=" + self._external + "\n")
            #if self._out:                   CONFIG.write("out=" + self._out + "\n")
            #if self._gOut:                  CONFIG.write("gOut=" + self._gOut + "\n")

class Parser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        modified_arg_line = "--" + arg_line
        for arg in modified_arg_line.split("="):
            if not arg.strip():
                continue
            yield arg

if __name__ == '__main__':
    s1=Options()
    s2=Options()
    if(id(s1)==id(s2)):
        print "Same"
    else:
        print "Different"

    s1.printToFile()