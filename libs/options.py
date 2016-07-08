#!/soft/devel/python-2.7/bin/python

import argparse
import os.path

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
        self._wd                    = options.wd if options.wd[-1]== "/" else options.wd + "/"
        self._minExpression         = options.minExpression
        self._annotation            = options.annotation
        self._tag                   = options.tag 
        self._only_models           = options.onlyModels
        self._specificDrivers       = "{}data/{}/specific_drivers/{}Drivers.txt".format(
                                        self._wd,self._annotation,self._tag)
        self._quickOut              = "{}analyses/{}/".format(self._wd,self._tag)
        self._parallelRange         = options.parallelRange
        self._parallelRange2         = options.parallelRange2
        self._externalSwitchesFile  = options.externalSwitchesFile
        self._parentTag             = options.parentTag
        # step for gene packs to parallelize
        self._step                  = 200

    # Getters ##
    @property
    def configFile(self):           return self._config_file
    @property
    def initialStep(self):          return self._initial_step
    @property
    def wd(self):                   return self._wd
    @property
    def minExpression(self):        return self._minExpression
    @property
    def annotation(self):            return self._annotation
    @property
    def tag(self):                  return self._tag
    @property
    def onlyModels(self):           return self._only_models
    @property
    def qout(self):                 return self._quickOut
    @property
    def specificDrivers(self):      return self._specificDrivers
    @property
    def filetag(self):
        filetag = ""

        if self.onlyModels:
            filetag += "_onlyModels"
        if self.parallelRange:
            filetag += "_{0}".format(self.parallelRange)

        return filetag

    @property
    def parallelRange(self):    return self._parallelRange
    @property
    def parallelRange2(self):    return self._parallelRange2
    @property
    def externalSwitchesFile(self): return self._externalSwitchesFile
    @property
    def parentTag(self): return self._parentTag
    @property
    def step(self): return self._step

    def parseOptions(self, *args, **kwds):
        parser = Parser(prog="smartas.py", description = "Find significant alternative splicing switches. Analyze their functional impact.", epilog= "Hector Climente, 2016", fromfile_prefix_chars='@')

        parser.add_argument('-f','--config-file', dest='config_file', action = 'store', 
                            help = 'Input file with the parameters')
        parser.add_argument('-s', '--initial-step', dest='initial_step', action='store', 
                            type=str, help='Where should SmartAS start.')
        parser.add_argument('-wd', '--working-directory', dest='wd', action='store', 
                            default='/home/hector/smartas/', 
                            help='Root file of SmartAS folder in the current machine.')
        parser.add_argument('-m', '--minimum-expression', dest='minExpression', action='store', 
                            default='-1', type=float, 
                            help='Minimum expression to consider a transcript not residual.')
        parser.add_argument('-i', '--annotation', dest='annotation', action='store', default='ucsc',
                            help='Used annotation.')
        parser.add_argument('-a', '--all-switches', dest='onlyModels', action='store_false',
                            help='Only use the model switches.')
        parser.add_argument('-t', '--tag', dest='tag' ,action='store',
                            help='Identifier of the analysis.')
        parser.add_argument('-d', '--specific-drivers', dest='specificDrivers', action='store', 
                            default='',help='Path of the specific drivers for the cancer type.')
        parser.add_argument('-p', '--parallel-range', dest='parallelRange', action='store', default='0',
                            type=int,help='Range of nodes if parallel.')
        parser.add_argument('-p2', '--scd-parallel-range', dest='parallelRange2', action='store', default='0',
                            type=int,help='2nd range of nodes if parallel.')
        parser.add_argument('-e', '--external-switches', dest='externalSwitchesFile', action='store', 
                            default=None,type=str, 
                            help='File containing switches calculated with other methods.')
        parser.add_argument('-x', '--parent-tag', dest='parentTag', action='store', default=None,
                            type=str, help='Tag of another experiment to partially use the information.')

        parser.set_defaults(onlyModels=True)
        options = parser.parse_args()
        if options.config_file!=None:
            if os.path.isfile(options.config_file):
                options = parser.parse_args(["@" + options.config_file])
        return options

    def printToFile(self,filename="",initialStep=None, wd=None, gwd=None, minExpression=None, annotation=None, tag=None, specificDrivers=None,parallelRange=None,onlyModels=None):
        """Print the config to a new file, only those values that are different 
        than the default ones. Overwrite those that are passed as arguments."""
        if not filename:
            cfgFilename = "{0}{1}.cfg".format(self.wd,self._tag)
        else:
            cfgFilename = "{0}{1}.cfg".format(self.wd,filename)

        with open(cfgFilename, "w") as CONFIG:
            if initialStep:
                CONFIG.write("initial-step=" + str(initialStep) + "\n")
            elif self._initial_step != 0:
                CONFIG.write("initial-step=" + str(self._initial_step) + "\n")
            
            if wd:
                CONFIG.write("working-directory=" + wd + "\n")
            elif self._wd != '/home/hector/smartas/':
                CONFIG.write("working-directory=" + self._wd + "\n")
            
            if minExpression:
                CONFIG.write("minimum-expression=" + str(minExpression) + "\n")
            elif self._minExpression != -1.0:         
                CONFIG.write("minimum-expression=" + str(self._minExpression) + "\n")
            
            if annotation:
                CONFIG.write("annotation=" + annotation + "\n")
            elif self._annotation != "ucsc":
                CONFIG.write("annotation=" + self._annotation + "\n")
            
            if tag:
                CONFIG.write("tag=" + tag + "\n")
            elif self._tag:
                CONFIG.write("tag=" + self._tag + "\n")

            if parallelRange:
                CONFIG.write("parallel-range={0}\n".format(parallelRange) )

            if onlyModels is not None and not onlyModels:
                CONFIG.write("all-switches\n" )
            
            if specificDrivers:
                CONFIG.write("specific-drivers=" + specificDrivers + "\n")
            elif self._specificDrivers:
                CONFIG.write("specific-drivers=" + self._specificDrivers + "\n")

        return cfgFilename

    def getCommandLineParameters(self,initialStep=None, wd=None, gwd=None, minExpression=None, annotation=None, tag=None, specificDrivers=None,parallelRange=None,onlyModels=None):
        """Print the config to a new file, only those values that are different 
        than the default ones. Overwrite those that are passed as arguments."""

        cfg = []

        if initialStep:
            cfg.extend(["--initial-step",str(initialStep)])
        elif self._initial_step != 0:
            cfg.extend(["--initial-step",str(self._initial_step)])
        
        if wd:
            cfg.extend(["--working-directory",wd])
        elif self._wd != '/home/hector/smartas/':
            cfg.extend(["--working-directory",self._wd])
        
        if minExpression:
            cfg.extend(["--minimum-expression",str(minExpression)])
        elif self._minExpression != -1.0:         
            cfg.extend(["--minimum-expression",str(self._minExpression)])
        
        if annotation:
            cfg.extend(["--annotation",annotation])
        elif self._annotation != "ucsc":
            cfg.extend(["--annotation",self._annotation])
        
        if tag:
            cfg.extend(["--tag",tag])
        elif self._tag:
            cfg.extend(["--tag",self._tag])

        if parallelRange:
            cfg.extend(["--parallel-range",str(parallelRange)])

        if onlyModels is not None and not onlyModels:
            cfg.extend(["--all-switches"])
        
        if specificDrivers:
            cfg.extend(["--specific-drivers",specificDrivers])
        elif self._specificDrivers:
            cfg.extend(["--specific-drivers",self._specificDrivers])

        return cfg

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
        print("Same")
    else:
        print("Different")

    s1.printToFile()