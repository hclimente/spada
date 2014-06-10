"""
    BIANA: Biologic Interactions and Network Analysis
    Copyright (C) 2009  Javier Garcia-Garcia, Emre Guney, Baldo Oliva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

"""
File        : psi_MiFormattedDB2biana.py
Author      : Javier Garcia & Emre Guney
Creation    : December 2007
Contents    : inserts information from a PSI-MI formatted XML file database into biana
Called from : 

=======================================================================================================

This file implements a program that fills up tables in database biana with information from a PSI-MI formatted database

--> Databases must be in PSI-MI XML format

"""

import os
from bianaParser import *
from psi_MiXMLParser import *

DICT_METHOD_CONVERSION_GRID_TO_PSI_MI = { 'Biochemical Activity': "MI:0401",
                                          'Co-crystal Structure': "MI:0114",
                                          'FRET': "MI:0055",
                                          'Co-localization': "MI:0403",
                                          'Co-purification': "MI:0025",
                                          'Invitro': "MI:0492",
                                          'Two-hybrid': "MI:0018",
                                          'Far Western': "MI:0047",
                                          'Invivo': "MI:0493",
                                          'Phenotypic Enhancement': "MI:0802",
                                          'Phenotypic Suppression': "MI:0796",
                                          'Affinity Capture-Western': "MI:0004", #i "MI:0113",
                                          'Co-fractionation': "MI:0027",
                                          'Affinity Capture-RNA': "MI:0004", #i ~"MI:0709"
                                          'Affinity Capture-MS': "MI:0004", #i "MI:0427"
                                          'Synthetic Rescue': "MI:0262",
                                          'Reconstituted Complex': "MI:0492",
                                          'Dosage Rescue': "MI:0261",
                                          'Protein-peptide': "MI:0084",
                                          'Affinity Capture-Luminescence': "MI:0004", #i ~"MI:0729" 
                                          'Protein-RNA': "MI:0316",
                                          'Synthetic Lethality': "MI:0441",
                                          'Dosage Growth Defect': "MI:0274",
                                          'Dosage Lethality': "MI:0441",
                                          'Synthetic Growth Defect': "MI:0274",
                                          'PCA': "MI:0090", # Protein Complementation assy
                                          'AffinityCapture-MS': "MI:0004" } # Only assigned affinity cromatography... should be inserted MS too?

MAX_NAME_LENGTH = 100

class Psi_MiFormattedDBParser(BianaParser):
    """
    PSI-MI formatted DB Parser Class
    """

    name = "psi_mi_2.5"
    description = "This parser inserts psi-mi 2.5 formated information to biana database"
    external_entity_definition = "Each relation participant is considered as a distinct External Entity"
    external_entity_relations = "External Entity Relations"
    
    
    dictDBNameToPrefix = {}
    #dictPrefixToDBName = {}
    def __init__(self):
        # Start with the default values
        BianaParser.__init__(self, default_db_description = "PSI-MI formatted protein-protein interaction database",
                             default_script_name = "psi_Mi25Parser.py",
                             default_script_description = "",
                             additional_compulsory_arguments = [("default-attribute=",None,"Name of the default identifier that this database gives (such as intact/mint/biogrid/dip/hprd/bind/mpact...)")])


        return

    def parse_database(self):
        """
        Method that implements the specific operations of PSI-MI formatted database parser
        """

        self.not_recognized_cross_refs = set()
        
        #directoryData = self.input_file[:self.input_file.rfind("/")+1]
        directoryData = os.path.dirname(self.input_file)
        command = None 
        onlyOneFileFlag = False 

        if os.path.isdir(self.input_file):
            command = None
            directoryData = os.path.dirname(self.input_file+os.sep)+os.sep
        elif os.path.isfile(self.input_file):
            directoryData = os.path.dirname(self.input_file)+os.sep
            if self.input_file.endswith(".zip"):
                command = "unzip"
            elif self.input_file.endswith(".gz"):
                command = "gunzip"
            elif self.input_file.endswith(".xml"):
                command = None
                onlyOneFileFlag = True
        else:
            sys.stderr.write("Warning: Input file extension (%s) not recognized by parser\n" % self.input_file[-3:])
            return
             
        if command is not None:
            os.chdir(directoryData)
            os.system("%s %s" % (command, self.input_file))
        
        if onlyOneFileFlag:
            listFileName = [self.input_file[self.input_file.rfind(os.sep)+1:]]
        else:
            #print directoryData
            listFileName = os.listdir(directoryData)

	self.file_number = 0

        parser = Psi_MiXMLParser(self.verbose)

        
        flagContinuePointReached = False


        for fileName in listFileName:
            
            sys.stderr.write("Parsing file %s\n" %fileName)

            if not (fileName.endswith(".xml") or fileName.endswith(".xsd.xml") or fileName.endswith(".mif25")):
                sys.stderr.write("Ignoring file: %s\n" % fileName)
                continue 
#            if not flagContinuePointReached:
#                if fileName == "BIOGRID-ORGANISM-Caenorhabditis_elegans-2.0.37.psi25.xml":
#                    flagContinuePointReached = True
#                continue

	    if self.time_control:
		if self.file_number%10==0:
			sys.stderr.write("%s files done in %s seconds\n" %(self.file_number, time.time()-self.initial_time))

	    self.file_number += 1
            
            if self.verbose:
		sys.stderr.write("\n------- %s\n" %fileName)

            # continue # to print just names
            try:
                parser.parseFile(directoryData+fileName)
            except Exception, inst:
                sys.stderr.write("%s\n" %inst)
            listEntry = parser.getEntries()
            
            psi_MiFormatted_object_number = 0
            for objEntry in listEntry:
                dictIdInteractorToIdExternal = {}
                dictExperiment = objEntry.getExperiments()
                dictInteractor = objEntry.getInteractors()
                dictInteraction = objEntry.getInteractions()
                
                if self.verbose:
                    sys.stderr.write("\nInteractors:\n")

                # Create external entities for interactors 
                for objInteractor in dictInteractor.itervalues():
                    if self.verbose:
                        sys.stderr.write("%s\n" %objInteractor.id)

                    ###if objInteractor.id != "350":
                    ###    continue
                    # Start new entry 
                    #print objInteractor.type.label
                    interactorType = self.decideInteractorTypeSpecificConversions(objInteractor.type.label)
                    if interactorType is None:
                        interactorType = self.decideInteractorTypeSpecificConversions(objInteractor.type.name)
                    psi_MiFormatted_object = ExternalEntity( source_database = self.database, type=interactorType) # "protein")
                    psi_MiFormatted_object_number += 1
                    # Fill the new entry
                    # Fill name
                    self.addNameAttributesToExternalEntityObject(objInteractor.name, psi_MiFormatted_object)
                    # Fill xRef
                    self.addXRefAttributesToExternalEntityObject(objInteractor.xRef, psi_MiFormatted_object)                    
                    # Fill taxId
                    if objInteractor.taxId is not None and int(objInteractor.taxId) >= 0:
                        psi_MiFormatted_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "taxid", value = objInteractor.taxId, type = "cross-reference"))
                    # Fill sequence
                    if objInteractor.sequence is not None:
                        sequenceType = self.decideSequenceTypeSpecificConversions(objInteractor.type.label)
                        if sequenceType == None:
                            sequenceType = self.decideSequenceTypeSpecificConversions(objInteractor.type.name)
                        #psi_MiFormatted_object.add_attribute(ExternalEntityAttribute("sequence","".join(objInteractor.sequence),"type" : sequenceType})
                    # Insert the entry to the database
                    self.biana_access.insert_new_external_entity( externalEntity = psi_MiFormatted_object )
                    dictIdInteractorToIdExternal[objInteractor.id] = psi_MiFormatted_object.get_id()

                if self.verbose:
                    sys.stderr.write("\nInteractions:\n")
                
		# Create external entity relations for interactions
                for objInteraction in dictInteraction.itervalues():
                    if self.verbose:
                        sys.stderr.write("%s\n" %objInteraction.id)

                    # Start new entry relation
                    if objInteraction.negative:
                        typeRelation = "no_interaction"
                    else:
                        typeRelation = "interaction"
                    psi_MiFormatted_object = ExternalEntityRelation( source_database=self.database, relation_type=typeRelation )
                    # Fill xRef
                    if objInteraction.xRef is not None:
                        self.addXRefAttributesToExternalEntityObject( objPsi_MiXRef= objInteraction.xRef, psi_MiFormatted_object=psi_MiFormatted_object, attribute_class=ExternalEntityRelationAttribute )
                    # Fill name
                    if objInteraction.name is not None:
                        self.addNameAttributesToExternalEntityObject(objInteraction.name, psi_MiFormatted_object)
                    # Fill experimentList
                    listObjXRefMethodParticipantIdentification = []
                    for idExperiment in objInteraction.listExperimentId:
                        experiment = dictExperiment[idExperiment]
                        # Fill experiment description - for now ignored --> add_common_attribute(intactExperiment) would return internal id assigned for each exp desription which would then be inserted as an attribute like methodID                        
                        #if experiment.description.name is not None: # description has no type ###self.addNameAttributesToExternalEntityObject(experiment.description, psi_MiFormatted_object, nameAttribute="description", flagIgnoreAlias=True)
                        #    psi_MiFormatted_object.add_attribute(attributeName="description", attributeFields={"value": experiment.description.name})
                        # Fill experiment bibref                        
                        self.addXRefAttributesToExternalEntityObject(experiment.xRefBib, psi_MiFormatted_object, flagIgnoreRefSecondary=True)
                        # Fill experiment xref - secondary references are ignored
                        if experiment.xRef is not None:
                            self.addXRefAttributesToExternalEntityObject(experiment.xRef, psi_MiFormatted_object, flagIgnoreRefSecondary=True)
                        # Fill experiment identification method
                        ###self.addXRefAttributesToExternalEntityObject(experiment.xRefMethodInteraction, psi_MiFormatted_object, flagIgnoreRefSecondary=True)
                        if experiment.xRefMethodInteraction.refPrimary.db == "psi-mi":
                            psi_MiFormatted_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "method_id", 
                                                                                                  value = experiment.xRefMethodInteraction.refPrimary.id[3:] ) )
                        if experiment.xRefMethodInteraction.refPrimary.db == "grid":
                            if DICT_METHOD_CONVERSION_GRID_TO_PSI_MI.has_key(experiment.nameMethodInteraction.label):
                                psi_MiFormatted_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier="method_id", 
                                                                                                      value = DICT_METHOD_CONVERSION_GRID_TO_PSI_MI[experiment.nameMethodInteraction.label][3:] ))
                            else:
                                sys.stderr.write("Method %s not recognized\n" %experiment.nameMethodInteraction.label)
                        ###else:
                        ###    print "Warning interaction type is not provided as psi-mi db reference:", experiment.xRefMethodInteraction.refPrimary.db
                        # Store participant identification method as xref in a list (method is the same for all participants in this interaction)
                        if experiment.xRefMethodParticipant is not None:
                            listObjXRefMethodParticipantIdentification.append(experiment.xRefMethodParticipant)
                    # Fill participantList
                    dictIdExternalToCardinality = {}
                    for participant in objInteraction.listParticipant:
                        try:
                            idExternal = dictIdInteractorToIdExternal[participant.interactorId]
                        except:
                            sys.stderr.write("Warning: Unassigned interactor %s\n" %participant.interactorId)
                            continue
                        flagFirstTime = insertKeyIntoHistogramDictionary(dictIdExternalToCardinality, idExternal)
                        if flagFirstTime: # need not to repeat same participant information                          
                            # Add new participant
                            psi_MiFormatted_object.add_participant( externalEntityID = idExternal )
                            # Fill participant identification methods using above created list 
                            for objXRefMethodIdentification in listObjXRefMethodParticipantIdentification:
                                psi_MiFormatted_object.add_participant_attribute(externalEntityID = idExternal, 
                                                                                  participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "detection_method",
                                                                                                                                                     value = objXRefMethodIdentification.refPrimary.id[3:]))
                        # Fill biological role
                        if participant.nameRoleBiological is not None:
                            nameRoleConverted = self.decideRoleSpecificConversions(participant.nameRoleBiological.label)
                            if nameRoleConverted != "ignore":
                                psi_MiFormatted_object.add_participant_attribute(externalEntityID = idExternal, 
                                                                                  participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", 
                                                                                                                                                     value = nameRoleConverted ))
                        # Fill experimental roles
                        for objNameRoleExperimental in participant.listNameRoleExperimental:
                            nameRoleConverted = self.decideRoleSpecificConversions(objNameRoleExperimental.label)
                            if nameRoleConverted != "ignore":
                                psi_MiFormatted_object.add_participant_attribute(externalEntityID = idExternal, 
                                                                                  participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", 
                                                                                                                                                     value = nameRoleConverted ))
                    for (idExternal, cardinality) in dictIdExternalToCardinality.iteritems():
                        psi_MiFormatted_object.add_participant_attribute(externalEntityID = idExternal, 
                                                                          participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "cardinality", 
                                                                                                                                             value = cardinality ))
                    # Fill interactionType - physical interaction for each - ignored for now
                    ###self.addXRefAttributesToExternalEntityObject(objInteraction.type, psi_MiFormatted_object, flagIgnoreRefSecondary=True)
                    # Insert the entry to the database
                    self.biana_access.insert_new_external_entity( externalEntity = psi_MiFormatted_object ) 
                    
        return
    
    def addNameAttributesToExternalEntityObject(self, objPsi_MiNames, psi_MiFormatted_object, nameAttribute="name", flagIgnoreAlias=False):
        if objPsi_MiNames.name is not None:
            #nameConverted = objPsi_MiNames.name.replace("&#150;", ' ')
            nameConverted = objPsi_MiNames.name.encode("ascii", "replace")
            if len(nameConverted) > MAX_NAME_LENGTH:
                psi_MiFormatted_object.add_attribute(ExternalEntityAttribute(attribute_identifier="description", value=nameConverted))    
            else:
                psi_MiFormatted_object.add_attribute(ExternalEntityAttribute(attribute_identifier = nameAttribute, value = nameConverted, type="cross-reference"))    
        if objPsi_MiNames.label is not None:
            if len(objPsi_MiNames.label) > MAX_NAME_LENGTH:
                psi_MiFormatted_object.add_attribute(ExternalEntityAttribute("description",objPsi_MiNames.label))
            else:
                psi_MiFormatted_object.add_attribute(ExternalEntityAttribute(attribute_identifier = nameAttribute, value = objPsi_MiNames.label, type="cross-reference"))#, "type": "label"})
        if not flagIgnoreAlias:
            if objPsi_MiNames.listAlias is not None:
                for (type, name) in objPsi_MiNames.listAlias:
                    attribute = self.decideAliasTypeSpecificConversions(type)
                    if attribute != "ignore" and name is not None:
                        psi_MiFormatted_object.add_attribute(ExternalEntityAttribute(attribute_identifier = attribute, value = name, type="alias"))#, "type": "alias"})
        return
    
    def addXRefAttributesToExternalEntityObject(self, objPsi_MiXRef, psi_MiFormatted_object, flagIgnoreRefSecondary=False, attribute_class=ExternalEntityAttribute):
        if objPsi_MiXRef.refPrimary is not None:
            (dbNameConverted, dictFieldName, dictFieldValue, dictFieldType) = self.decideDBReferenceSpecificConversions(objPsi_MiXRef.refPrimary)
            if dbNameConverted != "ignore":
                psi_MiFormatted_object.add_attribute(attribute_class(attribute_identifier = dbNameConverted, 
                                                                     value = dictFieldValue,
                                                                     type = dictFieldType) )
        if not flagIgnoreRefSecondary:                    
            if objPsi_MiXRef.listRefSecondary is not None:
                for objDBReference in objPsi_MiXRef.listRefSecondary:
                    (dbNameConverted, dictFieldName, dictFieldValue, dictFieldType) = self.decideDBReferenceSpecificConversions(objDBReference)
                    if dbNameConverted != "ignore":
                        psi_MiFormatted_object.add_attribute( attribute_class( attribute_identifier = dbNameConverted, 
                                                                               value = dictFieldValue,
                                                                               type = dictFieldType ) )
        return

    def decideInteractorTypeSpecificConversions(self, interactorType):
        interactorTypeConverted = None
        if interactorType == None:
            return interactorTypeConverted
        if interactorType == "protein":
            interactorTypeConverted = "protein"
        elif interactorType == "peptide":
            interactorTypeConverted = "protein"
        elif interactorType == "dna":
            interactorTypeConverted = "DNA"
        elif interactorType == "rna":
            interactorTypeConverted = "RNA"
        elif interactorType.endswith("dna"):
            interactorTypeConverted = "DNA"
        elif interactorType.endswith("rna"):
            interactorTypeConverted = "RNA"
        elif interactorType == "nucleic acid":
            interactorTypeConverted = "DNA"
        #elif interactorType == "mrna":
        #    interactorTypeConverted = "RNA"
        elif interactorType == "small molecule":
            interactorTypeConverted = "compound"
        else:
            sys.stderr.write("Warning: Unkown interactor type: %s\n" %interactorType)
        return interactorTypeConverted
    
    def decideSequenceTypeSpecificConversions(self, sequenceType):
        sequenceTypeConverted = None
        if sequenceType == None:
            return sequenceTypeConverted
        if sequenceType == "protein":
            sequenceTypeConverted = "peptide"
        elif sequenceType == "peptide":
            sequenceTypeConverted = "peptide"
        elif sequenceType == "dna":
            sequenceTypeConverted = "dna"
        elif sequenceType == "rna":
            sequenceTypeConverted = "rna"
        elif sequenceType.endswith("dna"):
            sequenceTypeConverted = "dna"
        elif sequenceType.endswith("rna"):
            sequenceTypeConverted = "rna"
        #elif sequenceType == "mrna":
        #    sequenceTypeConverted = "rna"
        #elif sequenceType == "ds dna":
        #    sequenceTypeConverted = "dna"
        elif sequenceType == "nucleic acid":
            sequenceTypeConverted = "dna"
        else:
            sys.stderr.write("Warning: Unkown sequence type: %s\n" %sequenceType)
        return sequenceTypeConverted

    def decideRoleSpecificConversions(self, nameRole):
        nameRoleConverted = ""
        if nameRole == "bait":
            nameRoleConverted = "bait"
        elif nameRole == "prey":
            nameRoleConverted = "prey"
        elif nameRole == "neutral component":
            nameRoleConverted = "neutral"
        elif nameRole == "unspecified role":
            nameRoleConverted = "ignore"
        elif nameRole == "unspecifiedrole":
            nameRoleConverted = "ignore"
        elif nameRole == "fluorescence acceptor":
            nameRoleConverted = "acceptor" 
        elif nameRole == "fluorescence accept":
            nameRoleConverted = "acceptor" 
        elif nameRole == "fluorescence donor":
            nameRoleConverted = "donor"
        elif nameRole == "self":
            nameRoleConverted = "self"
        elif nameRole == "ancillary":
            nameRoleConverted = "ancillary"
        elif nameRole == "enzyme":
            nameRoleConverted = "enzyme"
        elif nameRole == "enzyme target":
            nameRoleConverted = "enzyme target"
        elif nameRole == "inhibitor":
            nameRoleConverted = "inhibitor"
        elif nameRole == "cofactor":
            nameRoleConverted = "cofactor"
        elif nameRole == "electron acceptor":
            nameRoleConverted = "acceptor" 
        elif nameRole == "electron donor":
            nameRoleConverted = "donor"
        elif nameRole == "stimulator":
            nameRoleConverted = "stimulator"       
        else:
            sys.stderr.write("Warning: decideRoleSpecificConversions - Unknown type identifier: %s\n" %nameRole)
            nameRoleConverted = "ignore" 
        return nameRoleConverted

    def decideAliasTypeSpecificConversions(self, type):
        #nameType = ""
        attributeName = "ignore"
        if type == "gene name" or type == "gene name synonym":
            attributeName = "GeneSymbol"
            #nameType = "alias"
        elif type == "orf name":
            attributeName = "orfName"
            #nameType = "cross-reference"
        elif type == "locus name":
            attributeName = "OrderedLocusName"
            #nameType = "cross-reference"
        elif type == "isoform synonym":
            #attributeName = "isoFormName"
            attributeName = "ignore"
        else:
            sys.stderr.write("Warning: decideAliasTypeSpecificConversions - Unknown type identifier: %s\n" %type)
            
        return attributeName
    
    def getIndexOfFirstOccurenceOfDigit(self, strDBName):
        index = 0
        for char in strDBName:
            if ord(char) >= 48 and ord(char) <= 57:
                return index
            index += 1
        return -1
    
    def decideDBReferenceSpecificConversions(self, objDBReference):
        (db, id, type, secondary) = (objDBReference.db, objDBReference.id, objDBReference.type, objDBReference.secondary) 
        #dbNameConverted = None
        dbNameConverted = "ignore"
        dictFieldName = "value" # the default
        dictFieldValue = id # the default
        dictFieldType = "cross-reference" # the default 
        
        dbUpper = db.upper()
        
        if type == "identity":
            dictFieldType = "unique"
        
        if dbUpper == "UNIPROTKB" or dbUpper == "UNIPROT" or dbUpper == "UNIPROT KNOWLEDGE BASE" or dbUpper == "SWISSPROT" or dbUpper == "TREMBL":
            if id.startswith("unknown"):
                dbNameConverted = "ignore"
            else:
                dbNameConverted = "uniprotaccession"
                index = id.find("-PRO_")
                if index != -1:
                    #dictFieldValue = id[index+5:]
                    dictFieldValue = id[:index]
                    #self.checkOrInsertDBNamePrefix(dbNameConverted, "-PRO_")
                else:
                    index = id.find("NP_")
                    if index != -1:
                        dictFieldValue = id[index+3:]
                        self.checkOrInsertDBNamePrefix(dbNameConverted, "NP_")
                        #self.checkOrInsertDBNamePrefix(dbNameConverted+"2", "NP_")
        elif dbUpper == "INTENZ":
            dbNameConverted = "EC"
        elif dbUpper == "GO":
            dbNameConverted = "GO"
            dictFieldValue = id[3:]
            flagInconsistency = self.checkOrInsertDBNamePrefix(dbNameConverted, id[:3])
            # To correct the cases where the prefix is missing
            if flagInconsistency:
                if ord(id[0]) >= 48 and ord(id[0]) <= 57:
                    dictFieldValue = id
                elif id == "CC":
                    dbNameConverted = "ignore" 
                elif id.startswith("GO ") and ord(id[3]) >= 48 and ord(id[3]) <= 57:
                    dictFieldValue = id[3:]
        elif dbUpper == "INTERPRO":
            dbNameConverted = "interpro"
            dictFieldValue = id[3:]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:3])
        elif dbUpper == "ENSEMBL":
            dbNameConverted = "ensembl"
            # ENS[GBRM] -CGS - 
            #dictFieldValue = id[4:]
            #self.checkOrInsertDBNamePrefix(dbNameConverted, id[:4])
        elif dbUpper == "ENCODE":
            dbNameConverted = "encode"
	        # not always starts with AC
            #dictFieldValue = id[2:] 
            #self.checkOrInsertDBNamePrefix(dbNameConverted, id[:2])
        elif dbUpper == "INTACT":
            index = id.find("MINT-") # handling mint db's crayziness
            if index == -1:
                dbNameConverted = "IntAct"
                dictFieldValue = id[4:]
                self.checkOrInsertDBNamePrefix(dbNameConverted, id[:4])
            else:
                dbNameConverted = "MINT"
                dictFieldValue = id[5:]
                self.checkOrInsertDBNamePrefix(dbNameConverted, id[:5])
        elif dbUpper == "MIPS":
            dbNameConverted = "MIPS"
        elif dbUpper == "MINT":
            dbNameConverted = "MINT"
            dictFieldValue = id[5:]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:5])
        elif dbUpper == "PROTEIN ACCESSION":
            indexDigit = self.getIndexOfFirstOccurenceOfDigit(id)
            if indexDigit == 3:
                dbNameConverted = "AccessionNumber"
            elif indexDigit == 1:
                dbNameConverted = "UniprotAccession"
        elif dbUpper == "PROTEIN GI":
            dbNameConverted = "GI"
        elif dbUpper == "RCSB PDB" or dbUpper == "PDB" or dbUpper == "WWPDB":
            dbNameConverted = "pdb"
        elif dbUpper == "REACTOME COMPLEX" or dbUpper == "REACTOME PROTEIN" or dbUpper == "REACTOME":
            dbNameConverted = "Reactome"
            index = id.rfind('.')
            if index == -1:
                dictFieldValue = id[6:]
            else:
                dictFieldValue = id[6:index]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:6])
        elif dbUpper == "HUGE":
            dbNameConverted = "Huge"
            dictFieldValue = id[4:]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:4])
        elif dbUpper == "DDBJ-EMBL-GENBANK" or dbUpper == "DDBJ/EMBL/GENBANK" or dbUpper == "GENBANK_NUCLEOTIDE_G":
            dbNameConverted = "AccessionNumber"
            if self.sourcedb_name == "mint":
                dictFieldValue = secondary
        elif dbUpper == "GENBANK_PROTEIN_GI":
            if id.lower().startswith("gi:"):
                dictFieldValue = id[3:]
            dbNameConverted = "GI"
        elif dbUpper == "IPI":
            dbNameConverted = "IPI"
            #dictFieldValue = id[3:]
            #self.checkOrInsertDBNamePrefix(dbNameConverted, id[:3])
        elif dbUpper == "DIP":
            dbNameConverted = "DIP"
            dictFieldValue = id[4:]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:4])
        elif dbUpper == "WORMBASE": # wormbase, WormBase
            dbNameConverted = "wormbasegeneid"
            dictFieldValue = id[6:]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:6])
        elif dbUpper == "PUBMED":
    	    if id.startswith("unassigned"):
    	        dbNameConverted = "ignore"
            elif id.startswith("missing"): # "missing_pmid"
                dbNameConverted = "ignore"
    	    else:
                dbNameConverted = "pubmed"
        elif dbUpper == "UNIPARC":
            dbNameConverted = "uniparc"
            dictFieldValue = id[3:]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:3])
        elif dbUpper == "CHEBI":
            dbNameConverted = "chebi"
            dictFieldValue = id[6:]
            self.checkOrInsertDBNamePrefix(dbNameConverted, id[:6])
        elif dbUpper == "REFSEQ":
            dbNameConverted = "refseq"
            index = id.rfind('.')
            if index != -1:
                dictFieldValue = id[:index]
        elif dbUpper == "RGD":
            dbNameConverted = "rgd"
        elif dbUpper == "SGD":
            dbNameConverted = "SGD"
        elif dbUpper == "CYGD":
            dbNameConverted = "cygd"
        elif dbUpper == "FLYBASE":
            dbNameConverted = "FlyBase"
        elif dbUpper == "OMIM" or dbUpper == "MIM":
            dbNameConverted = "MIM"
        elif dbUpper == "INTENZ":
            dbNameConverted = "IntEnz"
        elif dbUpper == "ENTREZGENE":
            dbNameConverted = "geneID"
        elif dbUpper == "ENTREZ GENE/LOCUSLINK":
            dbNameConverted = "geneID"
        elif dbUpper == "HPRD":
            dbNameConverted = "HPRD"
        elif dbUpper == "HGNC":
            dbNameConverted = "HGNC"
        elif dbUpper == "MGI":
            dbNameConverted = "MGI"
        elif dbUpper == "TAIR":
            dbNameConverted = "TAIR"
        elif dbUpper == "RGD":
            dbNameConverted = "rgd"
        elif dbUpper == "RATMAP":
            dbNameConverted = "Ratmap"
        elif dbUpper == "IMGT/GENE-DB":
            dbNameConverted = "IMGT"
        elif dbUpper == "PSI-MI":
            dbNameConverted = "method_id"
        elif dbUpper == "DOI":
            dbNameConverted = "ignore"
        elif dbUpper == "CAMJEDB":
            dbNameConverted = "ignore"
        elif dbUpper == "ecogene":
            dbNameConverted = "ignore"
        elif dbUpper == "NEWT":
            dbNameConverted = "ignore"
        elif dbUpper == "IMEX":
            dbNameConverted = "ignore"
        elif dbUpper == "AFCS":
            dbNameConverted = "ignore"
        elif dbUpper == "PRIDE":
            dbNameConverted = "ignore"
        elif dbUpper == "SO":
            dbNameConverted = "ignore"
        elif dbUpper == "GRID" or dbUpper == "GRID_LEGACY":
            dbNameConverted = "ignore"
        elif dbUpper == "CDNA GI":
            dbNameConverted = "ignore"
        elif dbUpper == "CDNA ACCESSION":
            dbNameConverted = "ignore"
        elif dbUpper == "N/A":
            dbNameConverted = "ignore"
        else:
            if db not in self.not_recognized_cross_refs:
                sys.stderr.write("Warning: decideDBReferenceSpecificConversions - Unknown database identifier: %s" %db)
                self.not_recognized_cross_refs.add(db)
            #dbNameConverted = db.encode("ascii", "strict")
            dbNameConverted = "ignore"
            
            
        return (dbNameConverted, dictFieldName, dictFieldValue, dictFieldType)

    def checkOrInsertDBNamePrefix(self, dbName, prefix):
        flagInconsistency = False
        if Psi_MiFormattedDBParser.dictDBNameToPrefix.has_key(dbName):
            if Psi_MiFormattedDBParser.dictDBNameToPrefix[dbName] != prefix:
                sys.stderr.write("Warning: Database name prefix inconsistency: %s\t%s\n"  %(dbName, prefix))
                flagInconsistency = True
        else:
            Psi_MiFormattedDBParser.dictDBNameToPrefix[dbName] = prefix 
#        if Psi_MiFormattedDBParser.dictPrefixToDBName.has_key(prefix):
#            if Psi_MiFormattedDBParser.dictPrefixToDBName[prefix] != dbName:
#                if self.verbose:
#                    print "Warning: Database name prefix inconsistency", dbName, prefix
#                flagInconsistency = True
#        else:
#            Psi_MiFormattedDBParser.dictPrefixToDBName[prefix] = dbName
        return flagInconsistency
    
def insertKeyIntoHistogramDictionary(dictHistogram, key):
    "Inserts a key to a dictionary whish is designated to be used as histogram, returns True if this is the first occurence of the key"
    if dictHistogram.has_key(key):
        dictHistogram[key] += 1
        return False
    else:
        dictHistogram[key] = 1
        return True
     
 
