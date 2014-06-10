
"""
File        : uniprot2piana.py
Author      : Ramon Aragues, Javier Garcia Garcia
Creation    : 16.3.2004
Modified    : Javier Garcia Garcia December 2007
Contents    : fills up tables in database piana with information from uniprot
Called from : 

=======================================================================================================

This file implements a program that fills up tables in database piana with information of uniprot databases

This parser uses biopython libraries and methods

Command line option '--help' describes usage of this program

For more details on how to use it, read piana/README.populate_piana_db
"""


from bianaParser import *
import re, sys, sets

class GenericParser(BianaParser):
    """
    Generic Parser Class
    """

    name = "generic"
    description = "This file implements a program that fills up tables in BIANA database from data in a tabulated file"
    external_entity_definition = ""
    external_entity_relations = ""

    mandatory_columns = sets.Set(["id", "type"])
    mandatory_relation_columns = sets.Set(["id", "interactor_id_list", "type"])

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "Generic Tabulated parser",
                             default_script_name = "GenericParser.py",
                             default_script_description = GenericParser.description,
                             additional_compulsory_arguments = [("default-attribute=",None,"Name of the default identifier that this database gives (such as uniprotentry)")])
        
    def parse_database(self):
        """
        Method that implements the specific operations of a general tabulated file
        """
        
        value_separator = "|"
	participant_re = re.compile("(.+):\s*(.+)") #("(\w[(\s\w)]*):\s()")

        # Speficy that this database has relations hierarchies
        self.biana_access.store_relations_hierarchy = True

        self.initialize_input_file_descriptor()

        self.in_external_entities = False
        self.external_entity_fields = None
        self.external_entity_ids_dict = {}

        self.in_external_entity_relations = False
        self.external_entity_relation_fields = None
        
        for line in self.input_file_fd:

            line = line.strip()
            
            if line=="":
                continue
            
            if line.startswith("@EXTERNAL_ENTITY_DATA"):
                self.in_external_entities = True
                self.in_external_entity_relations = False
                continue
            elif line.startswith("@EXTERNAL_ENTITY_RELATION_DATA"):
                self.in_external_entities = False
                self.in_external_entity_relations = True
                continue
                
            # Parse external entities block
            if self.in_external_entities:
                if self.external_entity_fields is None:
                    values = re.split("\t+",line.strip())
		    column_to_index = dict([ (i.lower(),j) for i,j in zip(values, range(len(values))) ])
		    for x in self.mandatory_columns:
			if not column_to_index.has_key(x):
			    raise Exception("External Entity %s column not found" % x)
                    self.external_entity_fields = column_to_index 
                else:
                    values = re.split("\t+",line.strip())
                    #print values

                    new_external_entity = ExternalEntity( source_database = self.database,
                                                          type = values[self.external_entity_fields["type"]].strip() )
                    #for x in xrange(len(values)):
                    for x,i in self.external_entity_fields.iteritems():
			if x in self.mandatory_columns:
			    continue
                        if values[i].strip()!="-":
                            for current_value in values[i].split(value_separator):
				current_value = current_value.strip()
                                attribute_identifier = x 
                                if attribute_identifier.lower()=="proteinsequence":
                                    current_value = ProteinSequence(current_value)
                                new_external_entity.add_attribute( ExternalEntityAttribute( attribute_identifier= attribute_identifier, 
                                                                                            value=current_value, 
                                                                                            type="cross-reference") )
                                #print attribute_identifier, current_value
                          
                    self.external_entity_ids_dict[values[self.external_entity_fields["id"]]] = self.biana_access.insert_new_external_entity( externalEntity = new_external_entity )


            # Parse external entity relations block
            elif self.in_external_entity_relations:
                if self.external_entity_relation_fields is None:
                    values = re.split("\t+",line.strip())
		    column_to_index = dict([ (i.lower(),j) for i,j in zip(values, range(len(values))) ])
		    for x in self.mandatory_relation_columns:
			if not column_to_index.has_key(x):
			    raise Exception("External Entity Relation %s column not found" % x)
                    self.external_entity_relation_fields = column_to_index 
                else:
                    values = re.split("\t+",line.strip())
                    #print values
                    new_external_entity_relation = ExternalEntityRelation( source_database = self.database,
                                                                           relation_type = values[self.external_entity_relation_fields["type"]].strip() )
                    
		    for id in values[self.external_entity_relation_fields["interactor_id_list"]].split(value_separator):
			id = id.strip()
			new_external_entity_relation.add_participant( externalEntityID =  self.external_entity_ids_dict[id] )
                    
                    for current_attribute,index in self.external_entity_relation_fields.iteritems():
			if current_attribute in self.mandatory_relation_columns:
			    continue

                        v = values[index]
                        if v.strip()!="-":
			    if current_attribute.startswith("participants:"):
				current_attribute = current_attribute.replace("participants:", '')
				for current_value in v.split(value_separator):
				    s = participant_re.search(current_value.strip())
				    if s:
					participant_id = s.group(1)
					attribute_value = s.group(2)
				    else:
					sys.stderr.write("Format error, check file format!\n")

				    #print participant_id, current_attribute, attribute_value
				    new_external_entity_relation.add_participant_attribute( externalEntityID = self.external_entity_ids_dict[participant_id],
                                                                                            participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = current_attribute,
                                                                                                                                                               value = attribute_value ) )
			    else:
				for current_value in v.split(value_separator):
				    current_value = current_value.strip()
				    new_external_entity_relation.add_attribute( ExternalEntityRelationAttribute( attribute_identifier = current_attribute,
                                                                                                             value = current_value ) )
				    #print current_attribute, current_value
                    
                    self.external_entity_ids_dict[values[self.external_entity_relation_fields["id"]]] = self.biana_access.insert_new_external_entity( externalEntity = new_external_entity_relation )
            else:
                sys.stderr.write("Format error, check file format!\n")
                                

