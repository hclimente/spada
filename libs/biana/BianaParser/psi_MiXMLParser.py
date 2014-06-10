"""
File        : psi_MiXMLParser.py
Author      : Emre Guney & Javier Garcia 
Creation    : December 2007
Contents    : 
                - Psi_MiXMLParser object to parse an XML file in PSI-MI (PSI molecular interaction format) 
                - Psi_MiEntry object which corresponds to the "entry" in the file (container for interactors & interactions)
                - Psi_MiInteractor object 
                - Psi_MiInteraction object
                - Psi_MiParticipant object
                - Psi_MiNames object
                - Psi_MiXref object
                - DBReference object
                - Psi_MiHandler object (Content Handler for XML parsing with SAX)
                
Called from : psi_MiFormattedDB2piana.py

=======================================================================================================

Generic parser aiming to read data available in PSI-MI XML format in various databases
 
"""

#import xml.dom.minidom
from xml.sax import saxutils, handler, make_parser

#import time
#TIME_CHECK = 0 # time.clock()

"""
=======================================================================================================
                    Psi_MiXMLParser OBJECT
""" 

class Psi_MiXMLParser:
    """
    Class for parsing individual XML files obeying PSI-MI standarts
    """  
    
    def __init__(self, flagVerbose=False): #, fileName=None, listEntry=None):
        self.fileName = None
        self.file = None
        self.listEntry = []
        self.handler = Psi_MiHandler(flagVerbose)
        self.saxParser = make_parser()
        self.saxParser.setContentHandler(self.handler)
        return
    
    def __del__(self):
        if self.file is not None and not self.file.closed:
            self.file.close()
        return
    
    def __str__(self):
        return  "" 
    
    def parseFile(self, fileName=None):
        self.__init__() # first reset old contents
        if fileName is not None:
            self.fileName = fileName
            self.file = open(fileName)
        self.saxParser.parse(self.fileName)
        #objPsi_MiEntry = self._parseEntry(node)
        #self.addEntry(objPsi_MiEntry)
        self.listEntry = self.handler.listEntry
        if not self.file.closed:
            self.file.close()    
        return
    
    def addEntry(self, objPsi_MiEntry):
        self.listEntry.append(objPsi_MiEntry)
        return
    
    def getEntries(self):
        return self.listEntry

"""
=======================================================================================================
                    Psi_MiEntry OBJECT
""" 
class Psi_MiHandler(handler.ContentHandler):
    """
    Class to handle content in PSI-MI XML files
    """
    def __init__(self, flagVerbose):
        handler.ContentHandler.__init__(self)
        self.listTagStack = []
        self.listObjectStack = []
        #self.tagCurrentParent = None # either: experiment, interaction, interactor, ?? participant
        #self.objPsi_MiCurrent = None
        self.strCurrent = None
        self.tagCurrent = None
        self.strAttributeCurrent = None
        self.listEntry = []
        self.flagVerbose = flagVerbose
        return    

    def _printObjectStack(self):
        print self.listObjectStack
        return

    def _getAttribute(self, attrs, nameAttribute, flagReportInExistance=False):
        strAttribute = None
        for (name, value) in attrs.items():
            if name == nameAttribute:
                strAttribute = saxutils.escape(value)
        if strAttribute is None and flagReportInExistance:
            if self.flagVerbose:
                print "Warning: Attribute not found - %s", nameAttribute
        return strAttribute

    # ContentHandler methods     
    def startDocument(self):
        return

    def startElement(self, name, attrs):
        ##print name
        if name == Psi_MiEntry.PSI_MI_TAG_ENTRY: 
            self.listObjectStack.append(Psi_MiEntry())
        elif name == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION: 
            id = self._getAttribute(attrs, Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_ATTRIBUTE_ID, True)
            self.listObjectStack.append(Psi_MiExperiment(id))         
        elif name == Psi_MiInteractor.PSI_MI_TAG_INTERACTOR:
            id = self._getAttribute(attrs, Psi_MiInteractor.PSI_MI_TAG_INTERACTOR_ATTRIBUTE_ID, True)
            #print id,
            self.listObjectStack.append(Psi_MiInteractor(id)) 
        elif name == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
            id = self._getAttribute(attrs, Psi_MiInteraction.PSI_MI_TAG_INTERACTION_ATTRIBUTE_ID, True)
            self.listObjectStack.append(Psi_MiInteraction(id))
        elif name == Psi_MiParticipant.PSI_MI_TAG_PARTICIPANT:
            if self.listTagStack[-1] ==  Psi_MiInteraction.PSI_MI_TAG_PARTICIPANT_LIST:
                id = self._getAttribute(attrs, Psi_MiParticipant.PSI_MI_TAG_PARTICIPANT_ATTRIBUTE_ID, True)
                self.listObjectStack.append(Psi_MiParticipant(id))
        elif name == Psi_MiNames.PSI_MI_TAG_NAMES:
            # inserting many unnecessary objects (names inside which are not used but the other way is checking all coupled tag relationships)
            self.listObjectStack.append(Psi_MiNames()) 
        elif name == Psi_MiNames.PSI_MI_TAG_ALIAS:
            if self.listTagStack[-1] ==  Psi_MiNames.PSI_MI_TAG_NAMES:
                self.strAttributeCurrent = self._getAttribute(attrs, Psi_MiNames.PSI_MI_TAG_ALIAS_ATTRIBUTE_TYPE)
        elif name == Psi_MiXRef.PSI_MI_TAG_XREF:
            self.listObjectStack.append(Psi_MiXRef())
        elif name == Psi_MiXRef.PSI_MI_TAG_REF_PRIMARY:
            #if self.listTagStack[-1] ==  Psi_MiXRef.PSI_MI_TAG_XREF:
            objPsi_MiXRef = self.listObjectStack[-1]
            db = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_DB, True)
            id = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_ID, True)
            type = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_TYPE)
            secondary = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_SECONDARY)
            objPsi_MiXRef.refPrimary = DBReference(db, id, type, secondary)
        elif name == Psi_MiXRef.PSI_MI_TAG_REF_SECONDARY:
            objPsi_MiXRef = self.listObjectStack[-1]
            db = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_DB)
            id = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_ID)
            type = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_TYPE)
            secondary = self._getAttribute(attrs, Psi_MiXRef.PSI_MI_TAG_REF_ATTRIBUTE_SECONDARY)
            objPsi_MiXRef.listRefSecondary.append(DBReference(db, id, type, secondary))
        elif name == Psi_MiInteractor.PSI_MI_TAG_ORGANISM:
            if self.listTagStack[-1] == Psi_MiInteractor.PSI_MI_TAG_INTERACTOR:
                objPsi_MiInteractor = self.listObjectStack[-1]
                objPsi_MiInteractor.taxId = self._getAttribute(attrs, Psi_MiInteractor.PSI_MI_TAG_ORGANISM_ATTRIBUTE_TAX_ID, True)
        #self._printObjectStack()
        self.listTagStack.append(name)
        self.tagCurrent = name # to be able to remember the last tag inserted into stack and combine strings splitted inside the same tag by characters function  
        return

    def endElement(self, name):
        tagPrevious = self.listTagStack.pop()
        if name != tagPrevious:
            print "Warning: Tag inconsistency in stack", name          
        if name == Psi_MiEntry.PSI_MI_TAG_ENTRY: 
            objPsi_MiEntry = self.listObjectStack.pop()
            self.listEntry.append(objPsi_MiEntry)
        elif name == Psi_MiInteractor.PSI_MI_TAG_INTERACTOR:
            if self.listTagStack[-1] == Psi_MiEntry.PSI_MI_TAG_INTERACTOR_LIST and self.listTagStack[-2] == Psi_MiEntry.PSI_MI_TAG_ENTRY:
                objPsi_MiInteractor = self.listObjectStack.pop()
                objPsi_MiEntry = self.listObjectStack[-1]
                objPsi_MiEntry.addInteractor(objPsi_MiInteractor)
            elif self.listTagStack[-1] == Psi_MiParticipant.PSI_MI_TAG_PARTICIPANT and self.listTagStack[-2] == Psi_MiInteraction.PSI_MI_TAG_PARTICIPANT_LIST and self.listTagStack[-3] == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
                objPsi_MiInteractor = self.listObjectStack.pop()
                objPsi_MiParticipant = self.listObjectStack[-1]
                #objPsi_MiInteraction = self.listObjectStack[-2]
                objPsi_MiEntry = self.listObjectStack[-3]
                objPsi_MiEntry.addInteractor(objPsi_MiInteractor)
                objPsi_MiParticipant.interactorId = objPsi_MiInteractor.id
        elif name == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION:
            if self.listTagStack[-1] == Psi_MiEntry.PSI_MI_TAG_EXPERIMENT_LIST and self.listTagStack[-2] == Psi_MiEntry.PSI_MI_TAG_ENTRY:
                objPsi_MiExperiment = self.listObjectStack.pop()
                objPsi_MiEntry = self.listObjectStack[-1]
                objPsi_MiEntry.addExperiment(objPsi_MiExperiment)
            elif self.listTagStack[-1] == Psi_MiEntry.PSI_MI_TAG_EXPERIMENT_LIST and self.listTagStack[-2] == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
                objPsi_MiExperiment = self.listObjectStack.pop()
                objPsi_MiInteraction = self.listObjectStack[-1]
                objPsi_MiEntry = self.listObjectStack[-2]
                objPsi_MiEntry.addExperiment(objPsi_MiExperiment)
                objPsi_MiInteraction.listExperimentId.append(objPsi_MiExperiment.id)
        elif name == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
            objPsi_MiInteraction = self.listObjectStack.pop()
            objPsi_MiEntry = self.listObjectStack[-1]
            objPsi_MiEntry.addInteraction(objPsi_MiInteraction)
        elif name == Psi_MiParticipant.PSI_MI_TAG_PARTICIPANT:
            if self.listTagStack[-1] == Psi_MiInteraction.PSI_MI_TAG_PARTICIPANT_LIST:
                objPsi_MiParticipant = self.listObjectStack.pop()
                objPsi_MiInteraction = self.listObjectStack[-1]
                objPsi_MiInteraction.addParticipant(objPsi_MiParticipant)
        elif name == Psi_MiNames.PSI_MI_TAG_NAMES:
            #if isinstance(objPsi_MiCurrent, Psi_MiInteractor):
            objPsi_MiNames = self.listObjectStack.pop()
            if self.listTagStack[-1] ==  Psi_MiInteractor.PSI_MI_TAG_INTERACTOR:
                objPsi_MiInteractor = self.listObjectStack[-1]
                objPsi_MiInteractor.name = objPsi_MiNames
                #print objPsi_MiNames.name, objPsi_MiNames.label
            elif self.listTagStack[-1] == Psi_MiInteractor.PSI_MI_TAG_TYPE:
                objPsi_MiInteractor = self.listObjectStack[-1]
                objPsi_MiInteractor.type = objPsi_MiNames 
                #print objPsi_MiNames.name, objPsi_MiNames.label
            elif self.listTagStack[-1] == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION:
                objPsi_MiExperiment = self.listObjectStack[-1]
                objPsi_MiExperiment.description = objPsi_MiNames
            elif self.listTagStack[-1] == Psi_MiExperiment.PSI_MI_TAG_INTERACTION_DETECTION_METHOD and self.listTagStack[-2] == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION:
                objPsi_MiExperiment = self.listObjectStack[-1]
                objPsi_MiExperiment.nameMethodInteraction = objPsi_MiNames
            elif self.listTagStack[-1] == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
                objPsi_MiInteraction = self.listObjectStack[-1]
                objPsi_MiInteraction.name = objPsi_MiNames
            elif self.listTagStack[-1] == Psi_MiParticipant.PSI_MI_TAG_BIOLOGICAL_ROLE and self.listTagStack[-2] == Psi_MiParticipant.PSI_MI_TAG_PARTICIPANT:
                objPsi_MiParticipant = self.listObjectStack[-1]
                objPsi_MiParticipant.nameRoleBiological = objPsi_MiNames
            elif self.listTagStack[-1] == Psi_MiParticipant.PSI_MI_TAG_EXPERIMENTAL_ROLE and self.listTagStack[-2] == Psi_MiParticipant.PSI_MI_TAG_EXPERIMENTAL_ROLE_LIST and self.listTagStack[-3] == Psi_MiParticipant.PSI_MI_TAG_PARTICIPANT:
                objPsi_MiParticipant = self.listObjectStack[-1]
                objPsi_MiParticipant.listNameRoleExperimental.append(objPsi_MiNames)
        elif name == Psi_MiNames.PSI_MI_TAG_NAME:
            if self.listTagStack[-1] ==  Psi_MiNames.PSI_MI_TAG_NAMES:
                objPsi_MiNames = self.listObjectStack[-1]
                objPsi_MiNames.name = self.strCurrent
        elif name == Psi_MiNames.PSI_MI_TAG_LABEL:
            if self.listTagStack[-1] == Psi_MiNames.PSI_MI_TAG_NAMES:
                objPsi_MiNames = self.listObjectStack[-1]
                objPsi_MiNames.label = self.strCurrent
        elif name == Psi_MiNames.PSI_MI_TAG_ALIAS:
            if self.listTagStack[-1] == Psi_MiNames.PSI_MI_TAG_NAMES:
                objPsi_MiNames = self.listObjectStack[-1]
                objPsi_MiNames.listAlias.append((self.strAttributeCurrent, self.strCurrent))
        elif name == Psi_MiXRef.PSI_MI_TAG_XREF:
            objPsi_MiXRef = self.listObjectStack.pop()            
            if self.listTagStack[-1] == Psi_MiInteractor.PSI_MI_TAG_INTERACTOR:
                objPsi_MiInteractor = self.listObjectStack[-1]
                objPsi_MiInteractor.xRef = objPsi_MiXRef
            elif self.listTagStack[-1] == Psi_MiExperiment.PSI_MI_TAG_BIB_REF and self.listTagStack[-2] == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION:
                objPsi_MiExperiment = self.listObjectStack[-1]
                objPsi_MiExperiment.xRefBib = objPsi_MiXRef
            elif self.listTagStack[-1] == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION:
                objPsi_MiExperiment = self.listObjectStack[-1]
                objPsi_MiExperiment.xRef = objPsi_MiXRef
            elif self.listTagStack[-1] == Psi_MiExperiment.PSI_MI_TAG_INTERACTION_DETECTION_METHOD and self.listTagStack[-2] == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION:
                objPsi_MiExperiment = self.listObjectStack[-1]
                objPsi_MiExperiment.xRefMethodInteraction = objPsi_MiXRef
            elif self.listTagStack[-1] == Psi_MiExperiment.PSI_MI_TAG_PARTICIPANT_IDENTIFICATION_METHOD and self.listTagStack[-2] == Psi_MiExperiment.PSI_MI_TAG_EXPERIMENT_DESCRIPTION:
                objPsi_MiExperiment = self.listObjectStack[-1]
                objPsi_MiExperiment.xRefMethodParticipant = objPsi_MiXRef
            elif self.listTagStack[-1] == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
                objPsi_MiInteraction = self.listObjectStack[-1]
                objPsi_MiInteraction.xRef = objPsi_MiXRef               
            elif self.listTagStack[-1] == Psi_MiInteraction.PSI_MI_TAG_TYPE:
                objPsi_MiInteraction = self.listObjectStack[-1]
                objPsi_MiInteraction.type = objPsi_MiXRef
        elif name == Psi_MiInteractor.PSI_MI_TAG_SEQUENCE:
            if self.listTagStack[-1] == Psi_MiInteractor.PSI_MI_TAG_INTERACTOR:
                objPsi_MiInteractor = self.listObjectStack[-1]
                objPsi_MiInteractor.sequence = self.strCurrent
        elif name == Psi_MiInteraction.PSI_MI_TAG_EXPERIMENT_REFERENCE:
            if self.listTagStack[-1] == Psi_MiInteraction.PSI_MI_TAG_EXPERIMENT_LIST and self.listTagStack[-2] == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
                objPsi_MiInteraction = self.listObjectStack[-1]
                objPsi_MiInteraction.listExperimentId.append(self.strCurrent)
        elif name == Psi_MiInteraction.PSI_MI_TAG_NEGATIVE:
            if self.listTagStack[-1] == Psi_MiInteraction.PSI_MI_TAG_INTERACTION:
                objPsi_MiInteraction = self.listObjectStack[-1]
                objPsi_MiInteraction.negative = True
        elif name == Psi_MiParticipant.PSI_MI_TAG_INTERACTOR_REFERENCE:
            if self.listTagStack[-1] == Psi_MiParticipant.PSI_MI_TAG_PARTICIPANT:
                objPsi_MiParticipant = self.listObjectStack[-1]
                objPsi_MiParticipant.interactorId = self.strCurrent
        self.strCurrent = None
        #self._printObjectStack()
        return

    def characters(self, content):
        #if content != "" or content != " " or content != "\t":
        #    print "char:", saxutils.escape(content)
        if self.tagCurrent == self.listTagStack[-1]:
            #if content != "\n" and content != "\t" and content != " " and content != "      ": 
            #    self.strCurrent += content            
            if content.lstrip() != "":
                if self.strCurrent is None:
                    self.strCurrent = content.lstrip()
                else: 
                    self.strCurrent += content.lstrip()
        else:
            #self.strCurrent = content
            if content.lstrip() != "":                        
                self.strCurrent = content.lstrip()        
        return

    def ignorableWhitespace(self, content):
        #print content
        return
        
    def processingInstruction(self, target, data):
        #self._out.write('<?%s %s?>' % (target, data))
        print "Skipping Processing Instruction: %s %s" % (target, data)
        return
    
    def skippedEntity(self, name):
        print "Skipping Entity: %s" % name
        return


"""
=======================================================================================================
                    Psi_MiEntry OBJECT
""" 
class Psi_MiEntry:
    """
    Class representing an entry
    """
    PSI_MI_TAG_ENTRY = "entry"
    PSI_MI_TAG_EXPERIMENT_LIST = "experimentList"
    PSI_MI_TAG_INTERACTOR_LIST = "interactorList"
    PSI_MI_TAG_INTERACTION_LIST = "interactionList"
    
    def __init__(self, dictExperiment=None, dictInteractor=None, dictInteraction=None):
        if dictExperiment is None:
            self.dictExperiment = {}
        else:
            self.dictExperiment = dictExperiment
        if dictInteractor is None:
            self.dictInteractor = {}
        else:
            self.dictInteractor = dictInteractor
        if dictInteraction is None:
            self.dictInteraction = {}
        else:
            self.dictInteraction = dictInteraction
        return
    
    def __del__(self):
        return
    
    def __str__(self):
        return  "%s" % [k for k in self.dictInteractor.iterkeys()] 
    
    def addExperiment(self, objPsi_MiExperiment):
        self.dictExperiment[objPsi_MiExperiment.id] = objPsi_MiExperiment
        return True
    
    def addInteractor(self, objPsi_MiInteractor):
        self.dictInteractor[objPsi_MiInteractor.id] = objPsi_MiInteractor 
        return True
    
    def addInteraction(self, objPsi_MiInteraction):
        self.dictInteraction[objPsi_MiInteraction.id] = objPsi_MiInteraction 
        return True

    def getExperiments(self):
        return self.dictExperiment
    
    def getInteractors(self):
        return self.dictInteractor
    
    def getInteractions(self):
        return self.dictInteraction

"""
=======================================================================================================
                    Psi_MiExperiment OBJECT
"""
class Psi_MiExperiment:
    """
    Class representing an interactor
    """
    
    PSI_MI_TAG_EXPERIMENT_DESCRIPTION = "experimentDescription"
    PSI_MI_TAG_EXPERIMENT_ATTRIBUTE_ID = "id" 
    PSI_MI_TAG_BIB_REF = "bibref"
    PSI_MI_TAG_INTERACTION_DETECTION_METHOD = "interactionDetectionMethod"
    PSI_MI_TAG_PARTICIPANT_IDENTIFICATION_METHOD = "participantIdentificationMethod"
    
    def __init__(self, id=None, objPsi_MiNames=None, objPsi_MiXRefBib=None, objPsi_MiXRef=None, objPsi_MiXRefMethodInteraction=None, objPsi_MiXRefMethodParticipant=None, objPsi_MiNamesType=None):
        """
        id: in PSI_MI_TAG_EXPERIMENT_ATTRIBUTE_ID
        objPsi_MiXRefBib: inside PSI_MI_TAG_BIBREF inside PSI_MI_TAG_XREF  
        objPsi_MiXRef: inside PSI_MI_TAG_XREF
        objPsi_MiXRefMethodInteraction: inside PSI_MI_TAG_INTERACTION_DETECTION_METHOD
        objPsi_MiXRefMethodParticipant: inside PSI_MI_TAG_PARTICIPANT_IDENTIFICATION_METHOD
        
        """
        self.id = id # An integer represented as string
        self.description = objPsi_MiNames  
        self.xRefBib = objPsi_MiXRefBib 
        self.xRef = objPsi_MiXRef
        self.xRefMethodInteraction = objPsi_MiXRefMethodInteraction    # Interaction method as xRef
        self.xRefMethodParticipant = objPsi_MiXRefMethodParticipant
        self.nameMethodInteraction = objPsi_MiNamesType    # Interaction method as name
        return
    
    def __del__(self):
        return

    def __str__(self):
        return  "{%s: %s, %s, %s, %s}" % (self.id, self.xRefBib, self.xRef, self.xRefMethodInteraction, self.xRefMethodParticipant)


"""
=======================================================================================================
                    Psi_MiInteractor OBJECT
"""
class Psi_MiInteractor:
    """
    Class representing an interactor
    """
    
    PSI_MI_TAG_INTERACTOR = "interactor"
    PSI_MI_TAG_INTERACTOR_ATTRIBUTE_ID = "id" 
    PSI_MI_TAG_TYPE = "interactorType"
    PSI_MI_TAG_ORGANISM = "organism"
    PSI_MI_TAG_ORGANISM_ATTRIBUTE_TAX_ID = "ncbiTaxId"
    PSI_MI_TAG_SEQUENCE = "sequence"
    
    def __init__(self, id=None, objPsi_MiNames=None, objPsi_MiXRef=None, objPsi_MiNamesType=None, taxId=None, sequence=None): #, label=None, name=None, listAlias=None):
        """
        id: in PSI_MI_TAG_INTERACTOR_ATTRIBUTE_ID
        objPsi_MiNames: inside PSI_MI_TAG_NAMES 
        objPsi_MiXRef: inside PSI_MI_TAG_XREF
        type: inside PSI_MI_TAG_TYPE
        taxId: in PSI_MI_TAG_ORGANISM_ATTRIBUTE_TAX_ID
        sequence: inside PSI_MI_TAG_SEQUENCE
        """
        self.id = id # An integer represented as string 
        self.name = objPsi_MiNames 
        self.xRef = objPsi_MiXRef
        self.type = objPsi_MiNamesType
        self.taxId = taxId # An integer represented as string 
        self.sequence = sequence # A sequence of letters represented as string
        return
    
    def __del__(self):
        return

    def __str__(self):
        return  "{%s: %s, %s, %s, %s, %s}" % (self.id, self.name, self.xRef, self.type, self.taxId, self.sequence)

"""
=======================================================================================================
                    Psi_MiInteraction OBJECT
"""
class Psi_MiInteraction:
    """
    Class representing an interaction
    """
    PSI_MI_TAG_INTERACTION = "interaction"
    PSI_MI_TAG_INTERACTION_ATTRIBUTE_ID = "id" 
    PSI_MI_TAG_EXPERIMENT_LIST = "experimentList"
    PSI_MI_TAG_EXPERIMENT_REFERENCE = "experimentRef"
    PSI_MI_TAG_PARTICIPANT_LIST = "participantList"
    PSI_MI_TAG_TYPE = "interactionType"
    PSI_MI_TAG_NEGATIVE = "negative"
    
    def __init__(self, id=None, objPsi_MiNames=None, objPsi_MiXRef=None, listExperimentId=None, listObjPsi_MiParticipant=None, objPsi_MiXRefType=None, flagNegative=False):
        """
        id: in PSI_MI_TAG_INTERACTION_ATTRIBUTE_ID
        objPsi_MiNames: inside PSI_MI_TAG_NAMES 
        objPsi_MiXRef: inside PSI_MI_TAG_XREF
        listExperimentId: inside PSI_MI_TAG_EXPERIMENT_LIST inside PSI_MI_TAG_EXPERIMENT_REFERENCE
        listObjPsi_MiParticipant: inside PSI_MI_TAG_PARTICIPANT_LIST
        type: inside PSI_MI_TAG_TYPE
        """
        self.id = id # An integer represented as string 
        self.name = objPsi_MiNames 
        self.xRef = objPsi_MiXRef
        if listExperimentId is None: # A list of integers (represented as string) corresponding to experiment ids
            self.listExperimentId = []
        else:
            self.listExperimentId = listExperimenId
        if listObjPsi_MiParticipant is None:
            self.listParticipant = []
        else:
            self.listParticipant = listObjPsi_MiParticipant
        self.type = objPsi_MiXRefType    # Interaction type as xRef
        self.negative = flagNegative
        return
    
    def __del__(self):
        return

    def __str__(self):
        return  "{%s: %s, %s, %s}" % (self.id, self.name, self.xRef, self.type, self.negative)

    def addParticipant(self, objPsi_MiParticipant):
        self.listParticipant.append(objPsi_MiParticipant)
        return

    
"""
=======================================================================================================
                    Psi_MiParticipant OBJECT
"""
class Psi_MiParticipant:
    """
    Class representing a participant
    """
    PSI_MI_TAG_PARTICIPANT = "participant"
    PSI_MI_TAG_PARTICIPANT_ATTRIBUTE_ID = "id" 
    PSI_MI_TAG_INTERACTOR_REFERENCE = "interactorRef"
    PSI_MI_TAG_BIOLOGICAL_ROLE = "biologicalRole"
    PSI_MI_TAG_EXPERIMENTAL_ROLE_LIST = "experimentalRoleList"
    PSI_MI_TAG_EXPERIMENTAL_ROLE = "experimentalRole"
    
    def __init__(self, id=None, interactorReference=None, objPsi_MiNamesRoleBiological=None, listObjPsi_MiNamesRoleExperimental=None):
        """
        id: in PSI_MI_TAG_INTERACTOR_ATTRIBUTE_ID
        interactorId: inside PSI_MI_TAG_INTERACTOR_REFERENCE 
        objPsi_MiXRefRoleBiological: inside PSI_MI_TAG_BIOLOGICAL_ROLE inside PSI_MI_TAG_XREF
        listObjPsi_MiXRefExperimentalRole: inside PSI_MI_TAG_EXPERIMENTAL_ROLE_LIST inside PSI_MI_TAG_EXPERIMENTAL_ROLE inside PSI_MI_TAG_XREF 
        """
        self.id = id # An integer represented as string 
        self.interactorId = interactorReference 
        self.nameRoleBiological = objPsi_MiNamesRoleBiological
        if listObjPsi_MiNamesRoleExperimental is None: # A list of Psi_MiXRef objects corresponding to experiment roles
            self.listNameRoleExperimental = []
        else:
            self.listNameRoleExperimental = listObjPsi_MiNamesRoleExperimental
        return
    
    def __del__(self):
        return

    def __str__(self):
        return  "{%s: %s, %s}" % (self.id, self.interactorId, self.nameRoleBiological)

    
"""
=======================================================================================================
                    Psi_MiNames OBJECT
""" 
class Psi_MiNames:
    """
    Class representing information encapsulated within names XML tags
    """
    PSI_MI_TAG_NAMES = "names"
    PSI_MI_TAG_LABEL = "shortLabel"
    PSI_MI_TAG_NAME = "fullName"
    PSI_MI_TAG_ALIAS = "alias"
    PSI_MI_TAG_ALIAS_ATTRIBUTE_TYPE = "type"
    
    def __init__(self, label=None, name=None, listAlias=None):
        """ 
        shortLabel: inside PSI_MI_TAG_LABEL 
        fullName: inside PSI_MI_TAG_NAME
        listAlias: [(in PSI_MI_TAG_ALIAS_ATTRIBUTE_TYPE, inside PSI_MI_TAG_ALIAS), ] ---> i.e. [('gene name', 'Caf1'), ]
        PSI_MI_TAG_ALIAS_ATTRIBUTE_TYPE can be: gene name | gene name snonym | orf name  
          
        """
#        self.dictName[self.PSI_MI_TAG_LABEL] = label 
#        self.dictName[self.PSI_MI_TAG_NAME] = name
#        self.dictName[self.PSI_MI_TAG_ALIAS] = listAlias
        self.label = label 
        self.name = name
        if listAlias is None:
            self.listAlias = []
        else:
            self.listAlias = listAlias    # A list of tuples in the format (type_string, id_string)
        return
    
    def __del__(self):
        return
    
    def __str__(self):
        #return  "%s" % self.dictName[PSI_MI_TAG_LABEL]
        return  "%s" % self.label
        #return  "%s %s %s" % (self.label, self.label, self.listAlias)
      

"""
=======================================================================================================
                    Psi_MiXref OBJECT
""" 
class Psi_MiXRef:
    """
    Class representing information encapsulated within xref XML tags
    """
    PSI_MI_TAG_XREF = "xref"
    PSI_MI_TAG_REF_PRIMARY = "primaryRef"
    PSI_MI_TAG_REF_SECONDARY = "secondaryRef"
    PSI_MI_TAG_REF_ATTRIBUTE_DB = "db"
    PSI_MI_TAG_REF_ATTRIBUTE_ID = "id"
    PSI_MI_TAG_REF_ATTRIBUTE_TYPE = "refType"
    PSI_MI_TAG_REF_ATTRIBUTE_SECONDARY = "secondary"
    
    def __init__(self, objDBReferenceRefPrimary=None, listObjDBReferenceRefSecondary=None):
        """
        refPrimary: (in PSI_MI_TAG_REF_ATTRIBUTE_DB, in PSI_MI_TAG_REF_ATTRIBUTE_ID) 
        listRefSecondary: [(in PSI_MI_TAG_REF_ATTRIBUTE_DB, in PSI_MI_TAG_REF_ATTRIBUTE_ID), ] ---> i.e. [('go', 'GO:0035098')]
        """
        self.refPrimary = objDBReferenceRefPrimary
        if listObjDBReferenceRefSecondary is None:
            self.listRefSecondary = []
        else:
            self.listRefSecondary = listObjDBReferenceRefSecondary
        return
    
    def __del__(self):
        return
    
    def __str__(self):
        return "%s" % (self.refPrimary)
    
class DBReference:
    """
    Class representing reference information to a database with db, id and type fields
    """
    def __init__(self, db=None, id=None, type=None, secondary=None):
        self.db = db
        self.id = id
        self.type = type
        self.secondary = secondary  
        return
    
    def __del__(self):
        return
    
    def __str__(self):
        return "(%s, %s)" % (self.db, self.id)
              
  
