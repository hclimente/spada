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
File        : biopaxLevel2Parser.py
Author      : Javier Garcia Garcia
Creation    : June 2008
Contents    : fills up tables in database biana with information from a BIOPAX formatted database
Called from : 

=======================================================================================================
"""

from bianaParser import *
from xml.sax import saxutils, handler, make_parser
from XMLNode import XMLNode
import copy

class BiopaxEntity(object):

    resources = None
    controlled_relations = None
    
    database = None
    dbaccess = None
    _not_recognized = set()
    
    datatype_to_biana_type = { "uniprot": "uniprotaccession",
                               "ncbi_taxonomy": "taxID",
                               "tigr": "TIGR",
                               "tigr cna": "TIGR",
                               "tigr eha": "TIGR",
                               "tigr osa": "TIGR",
                               "reactome": "Reactome",
                               "chebi": "CHEBI",
                               "go": "GO",
                               "pubchem compound" : "PubchemCompound",
                               "glycan" : "keggCode",
                               "compound" : "keggCode",
                               #"entrez": "gi",
                               "entrez": "refseq",  # Previous sept 08 it was gi...
                               "pubmed": "pubmed",
                               "embl": "AccessionNumber",
                               "ensembl": "ensembl",
                               "wormbase": "wormbasesequencename",
                               "sgd": "sgd",
                               "flybase": "flybase" }

    def _identity(a):
        return a
    def _entrez_funct(a):
        # for gi:
        #return str(a).replace("gi|","")
        return str(a).split(".")
    def _reactome_funct(a):
        if a.startswith("REACT_"):
            return a[6:]
        else:
            return a

    datatype_operations = {}
    for x in datatype_to_biana_type:
        datatype_operations[x] = _identity

    datatype_operations["entrez"] =  _entrez_funct
    datatype_operations["reactome"] = _reactome_funct

    def __init__(self, XMLNode):
        
        self.rdf_id = XMLNode.attrs["rdf:ID"]

        self.synonyms = []
        self.comments = []
        self.data_source = None
        self.short_name = None
        self.availability = None
        self.name = None
        self.xrefs = []
        self.synonyms = []

        self.organism = None
        self.set_attributes(XMLNode)

        self.biana_object = None

    def set_attributes(self, XMLNode):
        """
        """
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:NAME":
                self.name = current_child.getValue()
            elif current_child.name == "bp:ORGANISM":
                self.organism = current_child.attrs["rdf:resource"]
            elif current_child.name == "bp:XREF":
                self.xrefs.append(current_child.attrs["rdf:resource"])
            elif current_child.name == "bp:SYNONYMS":
                self.synonyms.append(current_child.getValue())
            elif current_child.name == "bp:SHORT_NAME":
                self.short_name = current_child.getValue()
            elif current_child.name == "bp:COMMENT":
                self.comments.append(current_child.getValue())


    def _get_biana_object(self):
        raise ValueError("%s has not implemented _get_biana_object" %self)


    def toBiana( self ):
        """
        Add general attributes to external entity and adds it to the database
        
        returns the external entity id assigned to it
        """


        externalEntity = self._get_biana_object()

        if externalEntity is None:
            return

        if externalEntity.get_id() is not None:
            return externalEntity.get_id()

        externalEntity = self._get_biana_object()

        if externalEntity is None:
            return
        
        # Add name and short_name

        # PROBLEM: What to do with very large names? There are names with more than 255 characters... As this happens, name is inserted as a description, and short name as name
        if self.name is not None:
            externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value = self.name ) )  # SHOULD BE UNIQUE???! Of course, not as a description
            #externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "defaultID", value = self.name[0:255] ) )  # SHOULD BE UNIQUE???! Of course, not as a description
            
        if self.short_name is not None:
            externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value = self.short_name, type = "synonym" ) )


        # Add entity synonyms 
        # In reactome, synonyms seem to be geneSymbols.... Should they be inserted as gene symbols?
        [ externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value = x, type = "synonym" ) ) for x in self.synonyms ]

        # Add organism if defined
        if self.organism:
            organism_obj = BiopaxEntity.resources[self.organism]
            if organism_obj.tax_ref is not None:
                externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID", value = BiopaxEntity.resources[organism_obj.tax_ref].id, type = "cross-reference" ) )

        # Add comments as a description
        [ externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value = x ) ) for x in self.comments ]

        # Detect FUNCTION, CATALYTIC ACTIVITY, DISEASE, SUBCELLULAR LOCATION, SIMILARITY in comments field (at least for Reactome)
        split_regex = re.compile("(FUNCTION|CATALYTIC ACTIVITY|DISEASE|SUBCELLULAR LOCATION|SIMILARITY|DATABASE)")

        for current_comment in self.comments:
            t = split_regex.split(current_comment)
            for x in xrange(len(t)):
                if current_comment[x] == "FUNCTION":
                    externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "function", value = current_comment[x+1], type="cross-reference" ) )
                    x+=1
                elif current_comment[x] == "DISEASE":
                    externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "disease", value = current_comment[x+1], type="cross-reference" ) )
                    x+=1
                elif current_comment[x] == "SUBCELLULAR LOCATION":
                    externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "SubcellularLocation", value = current_comment[x+1], type="cross-reference" ) )
                    x+=1
                elif current_comment[x] == "CATALYTIC ACTIVITY":
                    # Where insert catalytic activity?
                    x+=1
        

        # Search for EC and MIM codes (maybe not necessary because they are defined as unificationXrefs?)
        ec_regex = re.compile("EC\s+(\d*\.\d*\.\d*.\d*)")
        mim_regex = re.compile("\[MIM\:(\d+)\]")

        if self.name is not None:
            for current_ec in ec_regex.findall(self.name):
                externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "ec", value = current_ec, type = "cross-reference" ) )

        if self.short_name is not None:
            for current_ec in ec_regex.findall(self.short_name):
                externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "ec", value = current_ec, type = "cross-reference" ) )

        for current_comment in self.comments:
            for current_ec in ec_regex.findall(current_comment):
                externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "ec", value = current_ec, type = "cross-reference" ) )
            for current_mim in mim_regex.findall(current_comment):
                externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "mim", value = current_mim, type = "cross-reference" ) )

        # Add all xrefs:
        for current_xref in self.xrefs:
            xref_object = BiopaxEntity.resources[current_xref]
            if xref_object.id is not None and xref_object.db is not None:
                if xref_object.db.lower() not in BiopaxEntity.datatype_to_biana_type:
                    if xref_object.db not in BiopaxEntity._not_recognized:
                        print xref_object.db, " not recognized"
                        BiopaxEntity._not_recognized.add(xref_object.db)
                else:
                    value = BiopaxEntity.datatype_operations[xref_object.db.lower()](xref_object.id)
                    if not isinstance(value,list):
                        externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = BiopaxEntity.datatype_to_biana_type[xref_object.db.lower()], 
                                                                               value = value,
                                                                               type = "cross-reference" ) )
                    else:
                        externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = BiopaxEntity.datatype_to_biana_type[xref_object.db.lower()], 
                                                                               value = value[0],
                                                                               version = value[1],
                                                                               type = "cross-reference" ) )


        # Check if this external entity (should be a relation) is catalyzed or modulated by other entities. If it is not a relation it will give an exception
        if BiopaxEntity.controlled_relations.has_key('#'+self.rdf_id):
            for current_controller in BiopaxEntity.controlled_relations['#'+self.rdf_id]:
                control_obj = BiopaxEntity.resources['#'+current_controller.rdf_id]
                #print "Current controller: %s" %current_controller.rdf_id
                #print "Getting %s for control obj %s" %(control_obj.controller_xref,control_obj.rdf_id)

                if control_obj.controller_xref is None:
                    print "control object %s has no contoller_xref" %control_obj.rdf_id
                    continue

                controller = BiopaxEntity.resources[control_obj.controller_xref]

                participant_eEid = controller.toBiana()
                
                if participant_eEid is None:
                    raise ValueError("In BiopaxEntity. %s" %controller)
                externalEntity.add_participant( externalEntityID =  participant_eEid )
                if control_obj.control_type is None:
                    raise ValueError("How is it possible to not have a controller role?")
                externalEntity.add_participant_attribute( externalEntityID = participant_eEid,
                                                           participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", value = control_obj.control_type ) )

        BiopaxEntity.dbaccess.insert_new_external_entity( externalEntity = externalEntity )

        #print "Going to insert ",self.name," to biana"
        #self.biana_external_entity_id = externalEntity.get_id()

        #return self.biana_external_entity_id

        #if self.rdf_id=="UniProt_P35354_Prostaglandin_G_H_synthase_2_precursor__EC_1_14_99_1___Cyclooxygenase__2___COX_2___Prostaglandin_endoperoxide_synthase_2___Prostaglandin_H2_synthase_2___PGH_synthase_2___PGHS_2___PHS_II_":
        #    print "PROTEIN PTGS2 INSERTED WITH ID ",externalEntity.get_id()

        return externalEntity.get_id()



class BiopaxPhysicalEntityParticipant(BiopaxEntity):
    """
    Biopax definition: any additional special characteristics of a physical entity in the context of an interaction or complex. These currently include stoichiometric coefficient and cellular location, but this list may be expanded in later levels.
    """

    def __init__(self, XMLNode):
        
        self.cellular_location_xref = None
        self.stoichiometric_coefficient = None
        self.physical_entity_xref = None
        self.sequence_features_list = []
        BiopaxEntity.__init__(self, XMLNode)
    
    def set_attributes(self, XMLNode):
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:CELLULAR-LOCATION":
                self.cellular_location_xref = current_child.attrs["rdf:resource"]
            elif current_child.name == "bp:PHYSICAL-ENTITY":
                self.physical_entity_xref = current_child.attrs["rdf:resource"]
            elif current_child.name == "bp:STOICHIOMETRIC-COEFFICIENT":
                self.stoichiometric_coefficient = current_child.getValue()
            elif current_child.name == "bp:SEQUENCE-FEATURE-LIST":
                # TO CHECK WHAT TO DO
                pass
        BiopaxEntity.set_attributes(self,XMLNode)

    def _get_biana_object(self):
        if self.biana_object is None:
            BiopaxEntity.toBiana(BiopaxEntity.resources[self.physical_entity_xref])
            self.biana_object = BiopaxEntity.resources[self.physical_entity_xref]._get_biana_object()
        return self.biana_object
        
    def add_participant_attributes_to_relation(self, eEr):
        """
        """
        if self.stoichiometric_coefficient is not None and self.stoichiometric_coefficient!="1" and self.stoichiometric_coefficient!=1:
            eEr.add_participant_attribute( externalEntityID = self.toBiana(),
                                            participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "cardinality", 
                                                                                                               value = self.stoichiometric_coefficient ) )
        
        if self.cellular_location_xref is not None:

            cellular_location_xref = BiopaxEntity.resources[self.cellular_location_xref]
            xref = BiopaxEntity.resources[cellular_location_xref.xref]
            if xref.db.lower()=="go":
                eEr.add_participant_attribute( externalEntityID = self.toBiana(),
                                                participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "go", value = xref.id ) )
            else:
                print "If cellular location is not in GO... in with controlled vocabulary is? It is %s" %xref.db


class BiopaxPathway(BiopaxEntity):
            
    def __init__(self, XMLNode):
        
        self.evidence = None
        self.pathway_components = []        
        BiopaxEntity.__init__(self, XMLNode)

    def set_attributes(self, XMLNode):
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:PATHWAY-COMPONENTS":
                #print "Adding component ",current_child.attrs["rdf:resource"]
                self.pathway_components.append(current_child.attrs["rdf:resource"])
        BiopaxEntity.set_attributes(self,XMLNode)


    def _get_biana_object(self):
        
        if self.biana_object is None:
            
            eEr = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "pathway" )

            for current_component in self.pathway_components:
                eE_ids_list = BiopaxEntity.resources[current_component].toBiana()
                if eE_ids_list is not None:
                    if isinstance(eE_ids_list,list):
                        for current_id in eE_ids_list:
                            if current_id is not None:
                                eEr.add_participant( externalEntityID = current_id )
                    else:
                        eEr.add_participant( externalEntityID = eE_ids_list )

                BiopaxEntity.resources[current_component].add_participants_to_eEr( eEr = eEr )

            self.biana_object = eEr

        return self.biana_object


class BiopaxPathwayStep(BiopaxEntity):
    """
    A step in a patwhay
    Multiple interactions may occur in a pathway step, each should be listed in the STEP-INTERACTIONS property.
    """

    def __init__(self, XMLNode):

        self.step_interactions_xref = []
        self.next_xref = None
        self.set_attributes(XMLNode)
        BiopaxEntity.__init__(self, XMLNode)

    def set_attributes(self, XMLNode):

        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:NEXT-STEP":  # NOT USED FOR THE MOMENT...
                self.next_xref = current_child.attrs["rdf:resource"]
            elif current_child.name == "bp:STEP-INTERACTIONS":
                self.step_interactions_xref.append(current_child.attrs["rdf:resource"])

    def toBiana(self):
        return [ BiopaxEntity.resources[current_step].toBiana() for current_step in self.step_interactions_xref ]


    def add_participants_to_eEr(self, eEr):
        for current_interaction in self.step_interactions_xref:
            if BiopaxEntity.resources[current_interaction].toBiana() is not None:
                eEr.add_participant( externalEntityID = BiopaxEntity.resources[current_interaction].toBiana() )
            else:
                ## For example, control elements
                pass
                #raise ValueError(current_interaction,BiopaxEntity.resources[current_interaction])

    def _get_biana_object(self):
        return None

    #def _get_biana_object(self):

    #    return [ BiopaxEntity.resources[current_step]._get_biana_object() for current_step in self.step_interactions_xref ]


class BiopaxInteraction(BiopaxEntity):
    
    def __init__(self, XMLNode):
        
        self.participants = None
        self.evidence = None
        BiopaxEntity.__init__(self)
        self.set_attributes(XMLNode)

    def set_attributes(self, XMLNode):
        # TODO
        pass

    def _get_biana_object(self):
        if self.biana_object is None:
            self.biana_object = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "interaction" )
        return self.biana_object


class BiopaxPhysicalInteraction(BiopaxEntity):

    def _get_biana_object(self):
        if self.biana_object is None:
            self.biana_object = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "interaction" )

        return self.biana_object


class BiopaxConversion(BiopaxPhysicalInteraction):

    def _get_biana_object(self):
        print "CONVERSION NOT WELL IMPLEMENTED YET!!!"
        if self.biana_object is None:
            self.biana_object = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "reaction" )
        return self.biana_object


class BiopaxComplexAssembly(BiopaxConversion):

    def toBiana(self):

        raise ValueError("COMPLEX ASSEMBLY NOT IMPLEMENTED YET")

        # TODO

        eEr = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "reaction" )  # COMPLEX AS AN INTERACTION?
        return BiopaxEntity.toBiana(self,eEr)


class BiopaxTransport(BiopaxConversion):
    
    def toBiana(self):

        raise ValueError("TRANSPORT NOT IMPLEMENTED YET")

        # TODO
        eEr = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "reaction" )  # COMPLEX AS AN INTERACTION?
        return BiopaxEntity.toBiana(self,eEr)


class BiopaxBiochemicalReaction(BiopaxConversion):
    
    def __init__(self, XMLNode):

        self.left_xrefs = []
        self.right_xrefs = []
        self.ec_number = []
        BiopaxEntity.__init__(self, XMLNode)
        self.set_attributes(XMLNode)

    def set_attributes(self, XMLNode):
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:LEFT":
                self.left_xrefs.append(current_child.attrs["rdf:resource"])
            elif current_child.name == "bp:RIGHT":
                self.right_xrefs.append(current_child.attrs["rdf:resource"])
            elif current_child.name == "bp:EC-NUMBER":
                self.ec_number.append(current_child.getValue())
        BiopaxEntity.set_attributes(self,XMLNode)

    def _get_biana_object(self):

        if self.biana_object is None:

            eEr = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "reaction" )
            [ eEr.add_attribute( ExternalEntityRelationAttribute( attribute_identifier = "ec", value = current_ec )) for current_ec in self.ec_number ] # ADD EC TO A REACTION? ADD A TYPE AS CROSS-REFERENCE?

            for current_left in self.left_xrefs:

                participant_obj = BiopaxEntity.resources[current_left]

                eEid = participant_obj.toBiana()
                if eEid is None:
                    raise ValueError(participant_obj)
                eEr.add_participant( externalEntityID = eEid )
                eEr.add_participant_attribute( externalEntityID = eEid,
                                                participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", value = "substrate" ) )
                participant_obj.add_participant_attributes_to_relation( eEr = eEr )

            for current_right in self.right_xrefs:

                participant_obj = BiopaxEntity.resources[current_right]
                eEid = participant_obj.toBiana()
                if eEid is None:
                    raise ValueError(participant_obj)
                eEr.add_participant( externalEntityID = eEid )
                eEr.add_participant_attribute( externalEntityID = eEid,
                                                participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", value = "product" ) )
                participant_obj.add_participant_attributes_to_relation( eEr = eEr )

            self.biana_object = eEr

        return self.biana_object


class BiopaxControl(BiopaxPhysicalInteraction):
    """
    An interaction in which one entity regulates, modifies, or otherwise influences another. Two types of control interactions are defined: activation and inhibition
    """

    # CONTROL Instances are not inserted as a relation itself, but forming part of the reactions they control

    def __init__(self, XMLNode):

        self.controller_xref = None
        self.controlled_xref = None
        self.control_type = None
        BiopaxEntity.__init__(self, XMLNode)
        self.set_attributes(XMLNode)

        if self.controlled_xref is not None:
            BiopaxEntity.controlled_relations.setdefault(self.controlled_xref,[]).append(self)

    def set_attributes(self, XMLNode):
        
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:CONTROLLER":
                self.controller_xref = current_child.attrs["rdf:resource"]
            elif current_child.name == "bp:CONTROLLED":
                self.controlled_xref = current_child.attrs["rdf:resource"]
            elif current_child.name == "bp:CONTROL-TYPE":
                control_type = current_child.getValue()
                if control_type == "ACTIVATION":
                    self.control_type = "activates"
                elif control_type == "INHIBITION":
                    self.control_type = "inhibits"
                elif control_type == "INHIBITION-ALLOSTERIC":
                    self.control_type = "allosteric_inhibition"
                elif control_type == "INHIBITION-COMPETITIVE":
                    self.control_type = "competitive_inhibition"
                elif control_type == "INHIBITION-IRREVERSIBLE":
                    self.control_type = "irreversible_inhibition"
                elif control_type == "INHIBITION-NONCOMPETITIVE":
                    self.control_type = "non_competitive_inhibition"
                elif control_type == "INHIBITION-OTHER":
                    self.control_type = "inhibits"
                elif control_type == "INHIBITION-UNCOMPETITIVE":
                    self.control_type = "uncompetitive_inhibition"
                elif control_type == "ACTIVATION-NONALLOSTERIC":
                    self.control_type = "nonallosteric_activation"
                elif control_type == "ACTIVATION-ALLOSTERIC":
                    self.control_type = "allosteric_activation"
                else:
                    raise ValueError("Control type %s not recognized" %control_type)
                    
        BiopaxPhysicalInteraction.set_attributes(self,XMLNode)

    def _get_biana_object(self):
        return None


class BiopaxCatalysis(BiopaxControl):
    """
    A control interaction in which a physical entity (a catalyst) increases the rate of a conversion interaction by lowering its activation energy. Instances of this class describe a pairing between a catalyzing entity and a catalyzed conversion
    """

    def __init__(self,XMLNode):
        self.cofactor = None
        self.direction = None
        BiopaxControl.__init__(self,XMLNode)
        self.set_attributes(XMLNode)

    def set_attributes(self, XMLNode):
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:COFACTOR":
                self.cofactor = current_child.attrs["rdf:resource"]
                # NOT IMPLEMENTED YET!!!
            elif current_child.name == "bp:DIRECTION":
                self.direction = current_child.getValue()
                # NOT IMPLEMENTED YET!!!
        BiopaxControl.set_attributes(self,XMLNode)


class BiopaxModulation(BiopaxControl):

    pass

        
class BiopaxPhysicalEntity(BiopaxEntity):

    def _get_biana_object(self):
        if self.biana_object is None:
            self.biana_object = ExternalEntity( source_database = BiopaxEntity.database, type = "compound" )
        return self.biana_object


class BiopaxComplex(BiopaxPhysicalEntity):

    # BIOPAX Complex element is inserted as a complex relation

    def __init__(self, XMLNode):
        self.components = []
        BiopaxPhysicalEntity.__init__(self,XMLNode)

    def set_attributes(self, XMLNode):

        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:COMPONENTS":
                self.components.append( current_child.attrs["rdf:resource"] )
            BiopaxPhysicalEntity.set_attributes(self,XMLNode)

    def _get_biana_object(self):

        if self.biana_object is None:

            eEr = ExternalEntityRelation( source_database = BiopaxEntity.database, relation_type = "complex" )
            for current_child in self.components:
                participant_obj = BiopaxEntity.resources[current_child]
                eEid = participant_obj.toBiana()
                eEr.add_participant( eEid )
                participant_obj.add_participant_attributes_to_relation( eEr = eEr )

            self.biana_object = eEr

        return self.biana_object


class BiopaxProtein(BiopaxPhysicalEntity):
        
    def __init__(self, XMLNode):
        BiopaxPhysicalEntity.__init__(self, XMLNode)

    def _get_biana_object(self):
        if self.biana_object is None:
            self.biana_object = ExternalEntity( source_database = BiopaxEntity.database, type = "protein" )
        return self.biana_object

class BiopaxDNA(BiopaxPhysicalEntity):
        
    def _get_biana_object(self):
        if self.biana_object is None:
            self.biana_object = ExternalEntity( source_database = BiopaxEntity.database, type = "dna" )
        return self.biana_object

class BiopaxRNA(BiopaxPhysicalEntity):
        
    def _get_biana_object(self):
        if self.biana_object is None:
            self.biana_object = ExternalEntity( source_database = BiopaxEntity.database, type = "rna" )
        return self.biana_object

class BiopaxSmallMolecule(BiopaxPhysicalEntity):

    def _get_biana_object(self):
        if self.biana_object is None:
            self.biana_object = ExternalEntity( source_database = BiopaxEntity.database, type = "compound" )
        return self.biana_object


class BiopaxXREF(BiopaxEntity):

    def __init__(self, XMLNode):
        self.db = None
        self.id = None
        BiopaxEntity.__init__(self, XMLNode)
        
    def set_attributes(self, XMLNode):
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:DB":
                self.db = current_child.getValue()
            elif current_child.name == "bp:ID":
                self.id = current_child.getValue()

    def _get_biana_object(self):
        return None


class BiopaxBioSource(BiopaxEntity):

    def __init__(self, XMLNode):
        self.celltype = None
        self.tissue = None
        self.name = None
        self.tax_ref = None
        BiopaxEntity.__init__(self, XMLNode)

    def set_attributes(self, XMLNode):
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:NAME":
                self.name = current_child.getValue()
            elif current_child.name == "bp:TAXON-XREF":
                self.tax_ref = current_child.attrs["rdf:resource"]
            else:
                print current_child.name," not recognized"

    def _get_biana_object(self):
        return None
    
class BiopaxDataSource(object):

    
    pass


class BiopaxOpenControlledVocabulary(BiopaxEntity):

    def __init__(self, XMLNode):
        self.term = None
        self.xref = None
        BiopaxEntity.__init__(self, XMLNode)

    def set_attributes(self, XMLNode):
        for current_child in XMLNode.getChilds():
            if current_child.name == "bp:TERM":
                self.term = current_child.getValue()
            elif current_child.name == "bp:XREF":
                self.xref = current_child.attrs["rdf:resource"]

    def _get_biana_object(self):
        return None


class BiopaxLevel2Parser(BianaParser):
    """

    """

    name = "biopax_level_2"
    description = "This file implements a program that fills up tables in database biana with information of a BIOPAX Level 2 formatted database"
    external_entity_definition = ""
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "Biopax formatted database",
                             default_script_name = "biopaxLevel2Parser.py",
                             default_script_description = BiopaxLevel2Parser.description,
                             additional_compulsory_arguments = [("default-attribute=",None,"Name of the default identifier that this database gives (such as reactome)")])


    class BiopaxLevel2Handler(handler.ContentHandler):
        """
        Class to handle content in Biopax Level2 XML files
        """

        def _identity(a):
            return a
        def _entrez_funct(a):
            return str(a).replace("gi|","")

        datatype_operations = { "uniprot": _identity,
                                "ncbi_taxonomy": _identity,
                                "tigr": _identity,
                                "reactome": _identity,
                                "chebi": _identity,
                                "go": _identity,
                                "pubchem compound" : _identity,
                                "glycan" : _identity,
                                "compound" : _identity,
                                "entrez": _entrez_funct }
        
        biopax_objects_dict = { "bp:unificationxref": BiopaxXREF,
                                "bp:relationshipxref": BiopaxXREF,
                                "bp:publicationxref": BiopaxXREF,
                                "bp:opencontrolledvocabulary": BiopaxOpenControlledVocabulary,
                                "bp:biosource": BiopaxBioSource,
                                "bp:protein": BiopaxProtein,
                                "bp:complex": BiopaxComplex,
                                "bp:dna": BiopaxDNA,
                                "bp:rna": BiopaxRNA,
                                "bp:smallmolecule": BiopaxSmallMolecule,
				"bp:physicalentity": BiopaxSmallMolecule,
                                "bp:pathway": BiopaxPathway,
                                "bp:interaction": BiopaxInteraction,
                                "bp:physcialinteraction": BiopaxPhysicalInteraction,
                                "bp:conversion": BiopaxConversion,
                                "bp:control": BiopaxControl,
                                "bp:biochemicalreaction": BiopaxBiochemicalReaction,
                                "bp:complexassembly": BiopaxComplexAssembly,
                                "bp:transport": BiopaxTransport,
                                "bp:catalysis": BiopaxCatalysis,
                                "bp:modulation": BiopaxModulation,
                                "bp:sequenceparticipant": BiopaxPhysicalEntityParticipant,
                                "bp:physicalentityparticipant": BiopaxPhysicalEntityParticipant,
                                "bp:pathwaystep": BiopaxPathwayStep }
        #"bp:transportwithbiochemicalreaction":,
         #                       "bp:physicalentityparticipant":,
          #                      "bp:proteinparticipant":,
           #                     "bp:complexparticipant":,
            #                    "bp:rnaparticipant":,
             #                   "bp:dnaparticipant":,
              #                  "bp:smallmoleculeparticipant": }
        
        def get_biana_data_type(self, type):
            return BiopaxLevel2Parser.BiopaxLevel2Handler.datatype_to_biana_type[type.lower()]

        def __init__(self):

            print "initalizing BiopaxLevel2Handler"

            self.current_XMLNode = None
            self.step = 0
            self.xmlnode_hierarchylist = []
            self.biopaxElements = {}
            
            BiopaxEntity.resources = self.biopaxElements
            BiopaxEntity.controlled_relations = {}

            handler.ContentHandler.__init__(self)

        def _get_cross_ref(self, xref_id):

            if xref_id[0] == '#':
                xref_id = xref_id[1:]
                
            try:
                return self.unification_xrefs[xref_id]
            except:
                return self._get_cross_ref( xref_id = self.recursive_xref[xref_id] )

        def _get_link(self, xref_id):
            
            if xref_id[0] == '#':
                xref_id = xref_id[1:]

            return self.links[ xref_id ]


        # ContentHandler methods     
        def startDocument(self):
            return

        def endDocument(self):
            self.step += 1

        def startElement(self, name, attrs):
            if self.current_XMLNode is None:
                self.current_XMLNode = XMLNode(name = name, attrs = attrs)
            else:
                t = XMLNode(name = name, attrs = attrs)
                self.current_XMLNode.addChild(t)
                self.xmlnode_hierarchylist.append(self.current_XMLNode)
                self.current_XMLNode = t

        def endElement(self, name):

            name = name.lower()

            if name!=self.current_XMLNode.name.lower():
                raise ValueError("ERROR IN XML FILE")
            
            if BiopaxLevel2Parser.BiopaxLevel2Handler.biopax_objects_dict.has_key(name):
                self.biopaxElements['#'+self.current_XMLNode.attrs["rdf:ID"]] = BiopaxLevel2Parser.BiopaxLevel2Handler.biopax_objects_dict[name](self.current_XMLNode)

            # Sets the current object the parent
            if len(self.xmlnode_hierarchylist)>0:
                self.current_XMLNode = self.xmlnode_hierarchylist.pop()
            else:
                self.current_XMLNode = None

            return

            
        def characters(self, text):
            #print text
            self.current_XMLNode.addValue(text.replace('.&lt;','<').replace('br&gt;','>').encode("ascii","ignore"))
            #self.current_XMLNode.addValue(text.replace('.&lt;','<').replace('br&gt;','>').decode("ascii","ignore"))


        def toBiana(self):
            for current_element in self.biopaxElements.values():
                current_element.toBiana()
                


    class BiopaxLevel2XMLParser(object):
        """
        Class for parsing individual XML files obeying BIOPAX Level 2 standards
        """

        def __init__(self, flagVerbose=False): #, fileName=None, listEntry=None):
            self.fileName = None
            self.file = None
            self.listEntry = []
            self.handler = BiopaxLevel2Parser.BiopaxLevel2Handler()
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
            self.handler.toBiana()
            if not self.file.closed:
                self.file.close()    
            return


    def parse_database(self):
        """
        """
        BiopaxEntity.database = self.database
        BiopaxEntity.dbaccess = self.biana_access


        # Speficy that this database has relations hierarchies
        self.biana_access.store_relations_hierarchy = True

        parser = self.BiopaxLevel2XMLParser(self.verbose)

        if os.path.isdir(self.input_file):
            files_list = os.listdir(self.input_file)
            if not self.input_file.endswith(os.sep):
                self.input_file += os.sep
            files_list = [ self.input_file+x for x in files_list ]
        else:
            files_list = [self.input_file]

        for current_file in files_list:
            if current_file.endswith(".owl"):
                print "Parsing file %s" %current_file
                it = time.time()
                parser.parseFile(current_file)
                if self.time_control:
                    print "Done in %s seconds" %(time.time()-it)




        
