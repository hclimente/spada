from database import Database,TableDB,FieldDB
import biana.BianaObjects


class BianaDatabase(Database):

    externalEntityID_col =  "externalEntityID"
    external_entity_relation_id_col = "externalEntityRelationID"
    externalEntityID_col_type = "integer(4) unsigned"
    externalDatabaseID_col = "externalDatabaseID"
    external_entity_relation_participant_id_col = "externalEntityRelationParticipantID"
    externalDatabaseID_col_type = "integer(2) unsigned"

    def __init__(self):

        Database.__init__(self)


        # DEFINE ALL TABLES = None

        self.DATABASE_VERSION_TABLE = None 
        self.BIANA_DATABASE_TABLE = None
        self.EXTERNAL_ENTITY_TABLE = None 
        self.EXTERNAL_DATABASE_AVAILABLE_eE_TYPES_TABLE = None 
        self.EXTERNAL_DATABASE_AVAILABLE_eEr_ATTRIBUTE_TABLE = None 
        self.EXTERNAL_DATABASE_AVAILABLE_eEr_TYPES_TABLE = None 
        self.EXTERNAL_DATABASE_AVAILABLE_eE_ATTRIBUTE_TABLE = None 
        self.EXTERNAL_DATABASE_TABLE = None 
        self.EXTERNAL_DATABASE_ATTRIBUTE_TRANSFER_TABLE = None 
        
        self.EXTERNAL_ENTITY_RELATION_TABLE = None 
        self.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE = None 
        self.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE = None 
        
        self.USER_ENTITY_TABLE = None 
        self.USER_ENTITY_PROTOCOL_TABLE = None 
        self.USER_ENTITY_PROTOCOL_ATOMS_TABLE = None 
        self.USER_ENTITY_PROTOCOL_ATOM_ATTRIBUTES_TABLE = None 
        
        self.ONTOLOGY_IS_A_TABLE = None 
        self.EXTENDED_ONTOLOGY_HIERARCHY_TABLE = None 
        self.ONTOLOGY_IS_PART_OF_TABLE = None 
        self.ONTOLOGY_INFO_TABLE = None 
        
        self.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES = None 
        self.EXTERNAL_ATTRIBUTE_DESCRIPTIONS_DICT = None 
        self.CD_HIT_CLUSTERING_TABLE = None 
        self.PROTEIN_BLAST_RESULTS_TABLE = None 

        self.TRANSFER_TABLES_LIST = []
        
        # DICTIONARIES
        self.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT = {}
        self.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT = {}
        
        # FIELDS
        self.externalEntityIdentifierType_field = None 

        self.externalEntityID_field = FieldDB(field_name = self.externalEntityID_col,
                                              data_type = "integer(4) unsigned",
                                              null = False)

        self.external_entity_relation_participant_id_field = FieldDB( field_name = self.external_entity_relation_participant_id_col,
                                                                      data_type = "integer(4) unsigned",
                                                                      null = False)
        self.externalDatabaseID_field = FieldDB(field_name = "externalDatabaseID",
                                                data_type = "integer(2) unsigned",
                                                null = False)


        # DATABASE SPECIFICIC DATA TYPES
        self.VALID_EXTERNAL_ENTITY_TYPES_DICT = {}
        self.VALID_EXTERNAL_ENTITY_RELATION_TYPES_DICT = {}

        self.VALID_IDENTIFIER_REFERENCE_TYPES_SET = set()

        self.VALID_EXTERNAL_ENTITY_ATTRIBUTE_TYPES_DICT = {}
        self.VALID_EXTERNAL_ENTITY_ATTRIBUTE_DATA_TYPES = {}
        self.VALID_EXTERNAL_ENTITY_GENERAL_ATTRIBUTE_TYPES_SET = set()
        self.VALID_EXTERNAL_ENTITY_IDENTIFIER_ATTRIBUTE_TYPES_SET = set()                    # Attributes of the type identifier that can be searched by value or by external entity 
        self.VALID_EXTERNAL_ENTITY_VERSIONABLE_IDENTIFIER_ATTRIBUTE_TYPES_SET =  set()       # Attributes of the type identifier with associated version that can be searched by value or by external entity and version
        self.VALID_EXTERNAL_ENTITY_DESCRIPTIVE_SEARCHABLE_ATTRIBUTE_TYPES_SET =  set()       # Descriptive attributes that can be searched by keywords
        self.VALID_EXTERNAL_ENTITY_DESCRIPTIVE_ATTRIBUTE_TYPES_SET = set()                   # Attributes searchable only by external entity
        self.VALID_EXTERNAL_ENTITY_NUMERIC_ATTRIBUTE_TYPES_SET = set()                       # Numeric attributes searchable by numeric value (equal, greater or lower than,...)
        self.VALID_EXTERNAL_ENTITY_SPECIAL_ATTRIBUTE_TYPES_SET = set()                       # Special attributes

        self.VALID_EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TYPES_DICT = {}
        self.VALID_EXTERNAL_ENTITY_RELATION_PARTICIPANT_NUMERIC_ATTRIBUTE_TYPES_SET = set()
        self.VALID_EXTERNAL_ENTITY_RELATION_PARTICIPANT_DESCRIPTIVE_ATTRIBUTE_TYPES_SET = set()
        self.VALID_EXTERNAL_ENTITY_RELATION_PARTICIPANT_DESCRIPTIVE_SEARCHABLE_ATTRIBUTE_TYPES_SET =  set()

        self._create_initial_tables()



    def add_valid_external_entity_relation_participant_attribute_type(self, current_attribute, data_type, category, additional_fields_tuple_list):

        new_table = None
        
        category = category.lower()

        if category == "eerp attribute":
            new_table = TableDB( table_name = "externalEntityRelationParticipant"+current_attribute,
                                 table_fields = [self.external_entity_relation_participant_id_field,
                                                 FieldDB( field_name = "value",
                                                          data_type = data_type,
                                                          null = False )],
                                 indices = [self.external_entity_relation_participant_id_col] )
        elif category == "eerp numeric attribute":
            pass
        elif category == "eerp descriptive attribute":
            pass
        elif category == "eerp descriptive searchable attribute":
            pass
        else:
            raise ValueError("Invalid category when adding valid relation participant attribute type: %s" %category)

        self.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT[current_attribute.lower()] = new_table


    def add_valid_external_entity_attribute_type(self, current_attribute, data_type, category, additional_fields_tuple_list=[]):
        
        category = category.lower()

        new_table = None

        if category == "ee attribute":
            self.VALID_EXTERNAL_ENTITY_GENERAL_ATTRIBUTE_TYPES_SET.add(current_attribute.lower())
            additional_fields_set = set("type")
            new_table =  TableDB( table_name = "externalEntity"+current_attribute,
                                  table_fields = [ FieldDB( field_name = "value",
                                                            data_type = data_type,
                                                            null = False ),
                                                   self.externalEntityID_field,
                                                   self.externalEntityIdentifierType_field ],
                                  primary_key = (self.externalEntityID_col, "value"),  # changed in order to optimize unification and expansions
                                  indices = [("value")] )
  

        elif category == "ee identifier attribute":
            self.VALID_EXTERNAL_ENTITY_IDENTIFIER_ATTRIBUTE_TYPES_SET.add(current_attribute.lower())
            additional_fields_set = set(["type"])
            new_table =  TableDB( table_name = "externalEntity"+current_attribute,
                                  table_fields = [ FieldDB( field_name = "value",
                                                            data_type = data_type,
                                                            null = False ),
                                                   self.externalEntityID_field,
                                                   self.externalEntityIdentifierType_field ],
                                  primary_key = (self.externalEntityID_col, "value"),  # changed in order to optimize unification and expansions
                                  indices = [("value")] )
            
        elif category == "ee versionable identifier attribute":
            self.VALID_EXTERNAL_ENTITY_VERSIONABLE_IDENTIFIER_ATTRIBUTE_TYPES_SET.add(current_attribute.lower())
            additional_fields_set = set(["type","version"])
            new_table =  TableDB( table_name = "externalEntity"+current_attribute,
                                  table_fields = [ FieldDB( field_name = "value",
                                                            data_type = data_type,
                                                            null = False ),
                                                   FieldDB( field_name = "version",
                                                            data_type = "integer(1) unsigned",
                                                            null = False,
                                                            default_value = 0 ),
                                                   self.externalEntityID_field,
                                                   self.externalEntityIdentifierType_field ],
                                  primary_key = (self.externalEntityID_col,"value","version"),
                                  indices = [ ("value",self.externalEntityID_col),
                                              ("value","version")] )


        elif category == "ee descriptive attribute":
            self.VALID_EXTERNAL_ENTITY_DESCRIPTIVE_ATTRIBUTE_TYPES_SET.add(current_attribute.lower())
            additional_fields_set = set()
            new_table = TableDB( table_name = "externalEntity"+current_attribute,
                                 table_fields = [ FieldDB( 	field_name = "value",
                                                                data_type = data_type,
                                                                null = False ),
                                                  self.externalEntityID_field,
                                                  self.externalEntityIdentifierType_field ],
                                 indices = [ (self.externalEntityID_col) ] )
            
            
        elif category == "ee descriptive searchable attribute":
            self.VALID_EXTERNAL_ENTITY_DESCRIPTIVE_SEARCHABLE_ATTRIBUTE_TYPES_SET.add(current_attribute.lower())
            additional_fields_set = set()
            new_table = TableDB( table_name = "externalEntity"+current_attribute,
                                 table_fields = [ FieldDB( field_name = "value",
                                                           data_type = data_type,
                                                           null = False ),
                                                  self.externalEntityID_field,
                                                  self.externalEntityIdentifierType_field ],
                                 indices = [ (self.externalEntityID_col) ],
                                 fulltext_indices = [ ("value") ] )
            
        elif category == "ee numeric attribute":
            self.VALID_EXTERNAL_ENTITY_NUMERIC_ATTRIBUTE_TYPES_SET.add(current_attribute.lower())
            additional_fields_set = set()
            new_table = TableDB( table_name = "externalEntity"+current_attribute,
                                 table_fields = [ FieldDB( field_name = "value",
                                                           data_type = data_type,
                                                           null = False ),
                                                  self.externalEntityID_field],
                                 indices = [ (self.externalEntityID_col,), ("value",) ] )
            
        elif category == "ee special attribute":
            self.VALID_EXTERNAL_ENTITY_SPECIAL_ATTRIBUTE_TYPES_SET.add(current_attribute.lower())
            additional_fields_set = set()
            fields = [ FieldDB( field_name = "value",
                                data_type = data_type,
                                null = False ),
                       self.externalEntityID_field,
                       self.externalEntityIdentifierType_field ]
            
            for current_field_name, current_data_type, current_null in additional_fields_tuple_list:
                if current_null == 0: current_null = False
                else:                 current_null = True
                fields.append( FieldDB( field_name = current_field_name,
                                        data_type = current_data_type,
                                        null = current_null ) )
                additional_fields_set.add(current_field_name)
                
            new_table = TableDB( table_name = "externalEntity"+current_attribute,
                                          table_fields = fields,
                                          indices = [ (self.externalEntityID_col,), ("value",) ] )

        else:
            raise ValueError("Invalid category when adding valid attribute type: %s" %category)
        
        self.VALID_EXTERNAL_ENTITY_ATTRIBUTE_TYPES_DICT[current_attribute.lower()] = {"name": current_attribute, "fields": additional_fields_set }
        self.VALID_EXTERNAL_ENTITY_ATTRIBUTE_DATA_TYPES[current_attribute.lower()] = data_type

        if new_table:
            self.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute.lower()] = new_table

    def is_valid_external_entity_relation_participant_attribute_type(self, type):       return type.lower() in self.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT
    def is_valid_external_entity_attribute_type(self, type):                            return type.lower() in self.VALID_EXTERNAL_ENTITY_ATTRIBUTE_TYPES_DICT
    def add_valid_external_entity_type(self, type):                                     self.VALID_EXTERNAL_ENTITY_TYPES_DICT[type.lower()] = type
    def get_valid_external_entity_types(self):                                          return self.VALID_EXTERNAL_ENTITY_TYPES_DICT.keys()
    def is_valid_external_entity_type(self, type):                                      return type.lower() in self.VALID_EXTERNAL_ENTITY_TYPES_DICT
    def add_valid_external_entity_relation_type(self, type):                            self.VALID_EXTERNAL_ENTITY_RELATION_TYPES_DICT[type.lower()] = type
    def get_valid_external_entity_relation_types(self):                                 return self.VALID_EXTERNAL_ENTITY_RELATION_TYPES_DICT.keys()
    def is_valid_external_entity_relation_type(self, type):                             return type.lower() in self.VALID_EXTERNAL_ENTITY_RELATION_TYPES_DICT
    def is_valid_identifier_reference_type(self, reference_type):                       return reference_type.lower() in self.VALID_IDENTIFIER_REFERENCE_TYPES_SET

    def get_attribute_data_type(self, attribute_identifier):                            return self.VALID_EXTERNAL_ENTITY_ATTRIBUTE_DATA_TYPES[attribute_identifier.lower()]

    def add_valid_identifier_reference_type(self, reference_type):

        self.VALID_IDENTIFIER_REFERENCE_TYPES_SET.add(reference_type.lower())
        valid_identifier_reference_types_list = list(self.VALID_IDENTIFIER_REFERENCE_TYPES_SET)
        valid_identifier_reference_types_list.sort()
        self.externalEntityIdentifierType_field = FieldDB( field_name = "type",
                                                           data_type = "ENUM(\"%s\")" %"\",\"".join(valid_identifier_reference_types_list),
                                                           user_friendly_name = "type" )



    def get_sql_query( self, ignore_primary_keys = False ):
        
        self._load_database_tables()
        return Database.get_sql_query(self, ignore_primary_keys = ignore_primary_keys )


    def get_tables(self):
        self._load_database_tables()
        return Database.get_tables(self)




    def _load_database_tables(self):
        """
        Returns a database object with all the content
        """

        self.remove_tables()
        self._create_initial_tables()
        self.create_specific_database_tables()

        self.add_table(self.DATABASE_VERSION_TABLE)
        self.add_table(self.BIANA_DATABASE_TABLE)
        self.add_table(self.TYPES_AND_ATTRIBUTES_TABLE)
        self.add_table(self.SPECIAL_ATTRIBUTES_TABLE)
        self.add_table(self.EXTERNAL_ENTITY_TABLE)
        self.add_table(self.EXTERNAL_DATABASE_AVAILABLE_eE_TYPES_TABLE)
        self.add_table(self.EXTERNAL_DATABASE_AVAILABLE_eEr_ATTRIBUTE_TABLE)
        self.add_table(self.EXTERNAL_DATABASE_AVAILABLE_eEr_TYPES_TABLE)
        self.add_table(self.EXTERNAL_DATABASE_AVAILABLE_eE_ATTRIBUTE_TABLE)
        self.add_table(self.EXTERNAL_DATABASE_TABLE)
        self.add_table(self.EXTERNAL_DATABASE_ATTRIBUTE_TRANSFER_TABLE)
        
        self.add_table(self.EXTERNAL_ENTITY_RELATION_TABLE)
        self.add_table(self.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE)
        self.add_table(self.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE)
        
        self.add_table(self.USER_ENTITY_TABLE)
        self.add_table(self.USER_ENTITY_PROTOCOL_TABLE)
        self.add_table(self.USER_ENTITY_PROTOCOL_ATOMS_TABLE)
        self.add_table(self.USER_ENTITY_PROTOCOL_ATOM_ATTRIBUTES_TABLE)
        
        self.add_table(self.ONTOLOGY_IS_A_TABLE)
        self.add_table(self.EXTENDED_ONTOLOGY_HIERARCHY_TABLE)
        self.add_table(self.ONTOLOGY_IS_PART_OF_TABLE)
        self.add_table(self.ONTOLOGY_INFO_TABLE) 
        
        for current_table in self.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES.values():
            self.add_table(current_table)

        self.add_table(self.CD_HIT_CLUSTERING_TABLE)
        self.add_table(self.PROTEIN_BLAST_RESULTS_TABLE)


        for current_table in self.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT.values():
            self.add_table(current_table)

        for current_table in self.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT.values():
            self.add_table(current_table)

        
    def _create_initial_tables(self):
        """
        Creates initial tables, not dependant on specific biana databases
        """        

        # First of all, create initial tables that are equal in ALL BIANA Databases
        self.BIANA_DATABASE_TABLE = TableDB( table_name = "BianaDatabase",
                                             table_fields = [ FieldDB( field_name="description",
                                                                       data_type = "varchar(255)" ),
                                                              FieldDB( field_name="last_externalEntityID",
                                                                       data_type = "integer(4) unsigned default 0",
                                                                       null = False),
                                                              FieldDB( field_name="last_externalEntityRelationParticipantID",
                                                                       data_type = "integer(4) unsigned default 0",
                                                                       null = False),
                                                              FieldDB( field_name="last_proteinSequenceID",
                                                                       data_type = "integer(4) unsigned default 0",
                                                                       null = False),
                                                              FieldDB( field_name="last_nucleotideSequenceID",
                                                                       data_type = "integer(4) unsigned default 0",
                                                                       null = False),
                                                              FieldDB( field_name="last_keyID",
                                                                       data_type = "integer(1) unsigned default 0",
                                                                       null = False),
                                                              FieldDB( field_name="source_code_version",
                                                                       data_type = "varchar(20)"),
                                                              FieldDB( field_name="optimized_for",
                                                                       data_type = "ENUM(\"parsing\",\"running\")" ) ] )
        
        self.TYPES_AND_ATTRIBUTES_TABLE = TableDB( table_name = "BianaDatabaseTypesAndAttributes",
                                                   table_fields = [ FieldDB( field_name="name",
                                                                             data_type = "varchar(255)",
                                                                             null = False ),
                                                                    FieldDB( field_name="data_type",
                                                                        data_type = "text" ), # not varchar because it only allows 255 characters
                                                                    FieldDB( field_name="category",
                                                                             data_type = "ENUM(\"ee type\",\"eer type\",\"ee attribute\",\"ee identifier attribute\",\"ee versionable identifier attribute\",\"ee descriptive attribute\",\"ee descriptive searchable attribute\",\"ee numeric attribute\",\"ee special attribute\",\"eerp attribute\",\"eerp descriptive attribute\",\"eerp descriptive searchable attribute\",\"identifier reference type\")",
                                                                             null = False ) ] )
        
        # Table for special attributes: those having more than a field
        self.SPECIAL_ATTRIBUTES_TABLE = TableDB( table_name = "BianaDatabaseSpecialAttributes",
                                                 table_fields = [ FieldDB( field_name="attribute_name",
                                                                           data_type = "varchar(255)",
                                                                           null = False ),
                                                                  FieldDB( field_name="field_name",
                                                                           data_type = "varchar(255)" ), # not varchar because it only allows 255 characters                                                             
                                                                  FieldDB( field_name="data_type",
                                                                           data_type = "text",
                                                                           null = False ),
                                                                  FieldDB( field_name="canbenull",
                                                                           data_type = "tinyint(1)",
                                                                           null = False ) ] )

        self.DATABASE_VERSION_TABLE = TableDB( table_name="database_version",
                                               table_fields = [ FieldDB(field_name = "versionID",
                                                                        data_type = "smallint unsigned auto_increment",
                                                                        null = False),
                                                                FieldDB(field_name = "dbControlID",
                                                                        data_type = "varchar(40)",
                                                                        null = False),
                                                                FieldDB(field_name = "date",
                                                                        data_type = "char(10)",
                                                                        null = False),
                                                                FieldDB( field_name = "stable_externalEntityID",
                                                                         data_type = "integer(4) unsigned"),
                                                                FieldDB( field_name = "stable_externalEntityRelationParticipantID",
                                                                         data_type = "integer(4) unsigned") ],
                                               primary_key = ("dbControlID"),
                                               indices = [("versionID")] )




    def create_specific_database_tables(self):
        """
        Adds to the database all the specific tables for this database
        """

        # GENERAL FIELDS
        
        eE_types_list = self.VALID_EXTERNAL_ENTITY_TYPES_DICT.keys()
        eE_types_list.sort()
        self.enum_eE_types_col = "ENUM(\"%s\")" %"\",\"".join(eE_types_list)

        eEr_types_list = self.VALID_EXTERNAL_ENTITY_RELATION_TYPES_DICT.keys()
        eEr_types_list.sort()
        self.enum_eEr_types_col = "ENUM(\"%s\")" %"\",\"".join(eEr_types_list)
    
        eE_attr_list = self.VALID_EXTERNAL_ENTITY_ATTRIBUTE_TYPES_DICT.keys()
        eE_attr_list.sort()
        self.enum_eEattr_col = "ENUM(\"%s\")" %"\",\"".join(eE_attr_list)


        self.EXTERNAL_ENTITY_TABLE = TableDB(table_name = "externalEntity",
                                             table_fields = [ FieldDB(field_name = self.externalEntityID_col,
                                                                      data_type = "integer(4) unsigned",
                                                                      null = False),
                                                              FieldDB(field_name = self.externalDatabaseID_col,
                                                                      data_type = self.externalDatabaseID_col_type,
                                                                      null = False),
                                                              FieldDB(field_name = "type",
                                                                      data_type = self.enum_eE_types_col,
                                                                      null = False)],
                                             primary_key = (self.externalDatabaseID_col,self.externalEntityID_col), #! if the order is reversed, the eEId index may not be needed
                                             indices = [ (self.externalEntityID_col) ] )

        self.EXTERNAL_DATABASE_AVAILABLE_eEr_TYPES_TABLE = TableDB( table_name = "externalDatabaseExternalEntityRelationType",
                                                                    table_fields = [ self.externalDatabaseID_field,
                                                                                     FieldDB( field_name = "eErType",
                                                                                              data_type = self.enum_eEr_types_col,
                                                                                              null = False) ],
                                                                    primary_key = (self.externalDatabaseID_col,"eErType") )

        self.EXTERNAL_DATABASE_AVAILABLE_eE_TYPES_TABLE = TableDB( table_name = "externalDatabaseExternalEntityType",
                                                                   table_fields = [ self.externalDatabaseID_field,
                                                                                    FieldDB( field_name = "eEType",
                                                                                             data_type = self.enum_eE_types_col,
                                                                                             null = False) ],
                                                                   primary_key = (self.externalDatabaseID_col,"eEType") )

        self.EXTERNAL_DATABASE_AVAILABLE_eEr_ATTRIBUTE_TABLE = TableDB( table_name = "externalDatabaseExternalEntityRelationAttribute",
                                                                          table_fields = [ self.externalDatabaseID_field,
                                                                                           FieldDB( field_name = "attributeType",
                                                                                                    data_type = self.enum_eEattr_col,
                                                                                                    null = False) ],
                                                                          primary_key = (self.externalDatabaseID_col,"attributeType") )

        self.EXTERNAL_DATABASE_AVAILABLE_eE_ATTRIBUTE_TABLE = TableDB( table_name = "externalDatabaseExternalEntityAttribute",
                                                                       table_fields = [ self.externalDatabaseID_field,
                                                                                        FieldDB( field_name = "attributeType",
                                                                                                 data_type = self.enum_eEattr_col,
                                                                                                 null = False) ],
                                                                       primary_key = (self.externalDatabaseID_col,"attributeType") )

        self.EXTERNAL_DATABASE_TABLE = TableDB( table_name = "externalDatabase",
                                                table_fields = [ FieldDB( field_name = self.externalDatabaseID_col,
                                                                          data_type = "smallint unsigned auto_increment",
                                                                          null = False),
                                                                 FieldDB( field_name = "databaseName",
                                                                          data_type = "varchar(255)",
                                                                          null = False),
                                                                 FieldDB( field_name = "databaseVersion",
                                                                          data_type = "varchar(255)",
                                                                          null = False),
                                                                 FieldDB( field_name = "parsedFile",
                                                                          data_type = "varchar(255)"),
                                                                 FieldDB( field_name = "parsedDate",
                                                                          data_type = "char(10)",
                                                                          null = False),
                                                                 FieldDB( field_name = "databaseDescription",
                                                                          data_type = "varchar(255)",
                                                                          null = False),
                                                                 FieldDB( field_name = "parsingTime",
                                                                          data_type = "integer(3) unsigned"),
                                                                 FieldDB( field_name = "defaultExternalEntityAttribute",
                                                                          data_type = self.enum_eEattr_col ),
                                                                 FieldDB( field_name = "isPromiscuous",
                                                                          data_type = "bool",
                                                                          default_value = 0,
                                                                          null = False) ], # + type_fields,
                                                primary_key = ("databaseName","databaseVersion"),
                                                indices = [self.externalDatabaseID_col])

        self.EXTERNAL_DATABASE_ATTRIBUTE_TRANSFER_TABLE = TableDB( table_name = "externalDatabaseAttributeTransfer",
                                                                   table_fields = [ FieldDB( field_name = "keyID",
                                                                                             data_type = "integer(1) unsigned",
                                                                                             null = False ),
                                                                                    FieldDB( field_name = "externalDatabaseID",
                                                                                             data_type = "integer(2) unsigned",
                                                                                             null = False ),
                                                                                    FieldDB( field_name = "attributeKey",
                                                                                             data_type = self.enum_eEattr_col,
                                                                                             null = False ),
                                                                                    FieldDB( field_name = "transferAttribute",
                                                                                             data_type = self.enum_eEattr_col,
                                                                                             null = False ) ] )
        
        self.ALIGNMENT_ELEMENT_TABLE = TableDB( table_name = "alignmentElement",
                                                table_fields = [ self.externalEntityID_field,                                   # Referring to the external Entity Relation ID of the alignment
                                                                 FieldDB( field_name = "sequenceMD5",
                                                                          data_type = "binary(16)"),                       # Previously it was mandatory, but as not always is mandatory to have previously the sequence, now it is not mandatory
                                                                 FieldDB( field_name = "position",
                                                                          data_type = "integer(3) unsigned",
                                                                          null = False ),
                                                                 FieldDB( field_name = "crossID",           # Used to cross the alignment element with another external entity (as pex an uniprot, where then it is possible to obtain complete sequence, taxonomy,...)
                                                                          data_type = "varchar(255)" ),                                                            
                                                                 FieldDB( field_name = "crossID_identifier_type",
                                                                          data_type = self.enum_eEattr_col ),
                                                                 FieldDB( field_name = "alignmentBlocks",
                                                                          data_type = "varchar(255)"),
                                                                 FieldDB( field_name = "alignedSequence",                  # Compressed aligned sequence (used for alignments already done)
                                                                          data_type = "text"),
                                                                 FieldDB( field_name = "identity",
                                                                          data_type = "integer(1) unsigned"),
                                                                 FieldDB( field_name = "weighted_sim",
                                                                          data_type = "integer(1) unsigned"),
                                                                 FieldDB( field_name = "sequence_covered",
                                                                          data_type = "integer(1) unsigned")],
                                                primary_key = (self.externalEntityID_col,"sequenceMD5","position"),
                                                indices = [("sequenceMD5")] )
        
        
        self.EXTERNAL_ENTITY_RELATION_TABLE = TableDB( table_name = "externalEntityRelation",
                                                       table_fields = [ FieldDB( field_name = self.external_entity_relation_id_col,
                                                                                 data_type = self.externalEntityID_col_type,
                                                                                 null = False ),
                                                                        FieldDB(field_name = "type",
                                                                                data_type = self.enum_eEr_types_col,
                                                                                null = False)],
                                                       primary_key = (self.external_entity_relation_id_col))


        self.external_entity_relation_id_field = FieldDB( field_name = self.external_entity_relation_id_col,
                                                          data_type = "integer(4) unsigned",
                                                          null = False,
                                                          foreign_key = (self.EXTERNAL_ENTITY_RELATION_TABLE,
                                                                         self.EXTERNAL_ENTITY_RELATION_TABLE.get_field(field_name="externalEntityRelationID" )))


        self.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE = TableDB( table_name = "externalEntityRelationParticipant",
                                                                   table_fields = [ FieldDB( field_name = self.external_entity_relation_participant_id_col,
                                                                                             data_type = "integer(4) unsigned",
                                                                                             null = False ),
                                                                                    self.external_entity_relation_id_field,
                                                                                    self.externalEntityID_field],
                                                                   primary_key = self.external_entity_relation_participant_id_col,
                                                                   indices = [(self.external_entity_relation_id_col, self.externalEntityID_col),
                                                                              (self.externalEntityID_col, self.external_entity_relation_id_col)])
        
        self.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE = TableDB( table_name = "extendedExternalEntityRelationParticipant",
                                                                            table_fields = [ FieldDB( field_name = self.external_entity_relation_participant_id_col,
                                                                                                      data_type = "integer(4) unsigned",
                                                                                                      null = False ),
                                                                                             self.external_entity_relation_id_field,
                                                                                             self.externalEntityID_field],
                                                                            primary_key = self.external_entity_relation_participant_id_col,
                                                                            indices = [(self.external_entity_relation_id_col, self.externalEntityID_col),
                                                                                       (self.externalEntityID_col,self.external_entity_relation_id_col)])
        



        ##################################################
        ###         USER ENTITIES RELATED              ###
        ##################################################

        # A table to store the information about how the unification has been done
        # A table to store the unified entities

        self.USER_ENTITY_TABLE = TableDB( table_name = "userEntityUnification_protocol_",
                                          table_fields = [ FieldDB( field_name = "userEntityID",
                                                                    data_type = "integer(4) unsigned",
                                                                    null = False ),
                                                           FieldDB( field_name = "externalEntityID",
                                                                    data_type = "integer(4) unsigned",
                                                                    null = False) ],
                                          primary_key = ("userEntityID","externalEntityID"),
                                          indices = [("externalEntityID","userEntityID")])
        
        self.USER_ENTITY_PROTOCOL_TABLE = TableDB( table_name = "userEntityUnificationProtocol",
                                                     table_fields = [ FieldDB( field_name = "unificationProtocolID",
                                                                               data_type = "smallint unsigned auto_increment",
                                                                               null = False),
                                                                      FieldDB( field_name = "description",
                                                                               data_type = "varchar(255)",
                                                                               null = False),
                                                                      FieldDB( field_name = "databaseVersion",
                                                                               data_type = "varchar(40)",
                                                                               null = False )],
                                                     # deleted foreign key for the moment... (because it failed when inserting default protocol)
                                                     #foreign_key = (DATABASE_VERSION_TABLE,
                                                     #               DATABASE_VERSION_TABLE.get_field(field_name = "dbControlID")))],
                                                     primary_key = "unificationProtocolID" )


        self.USER_ENTITY_PROTOCOL_ATOMS_TABLE = TableDB( table_name = "userEntityUnificationProtocolAtoms",
                                                           table_fields = [ FieldDB( field_name = "unificationAtomID",
                                                                                     data_type = "integer(2) unsigned auto_increment",
                                                                                     null = False),
                                                                            FieldDB( field_name = "unificationProtocolID",
                                                                                     data_type = "smallint unsigned",
                                                                                     null = False,
                                                                                     foreign_key = (self.USER_ENTITY_PROTOCOL_TABLE,
                                                                                                    self.USER_ENTITY_PROTOCOL_TABLE.get_field(field_name = "unificationProtocolID") )),
                                                                            FieldDB( field_name = "externalDatabaseID_A",
                                                                                     data_type = "integer(2) unsigned",
                                                                                     null = False,
                                                                                     foreign_key = (self.EXTERNAL_DATABASE_TABLE,
                                                                                                    self.EXTERNAL_DATABASE_TABLE.get_field(field_name = "externalDatabaseID") )),
                                                                            FieldDB( field_name = "externalDatabaseID_B",
                                                                                     data_type = "integer(2) unsigned",
                                                                                     null = False,
                                                                                     foreign_key = (self.EXTERNAL_DATABASE_TABLE,
                                                                                                    self.EXTERNAL_DATABASE_TABLE.get_field(field_name = "externalDatabaseID") )) ],
                                                           indices = ["unificationProtocolID"],
                                                           primary_key = "unificationAtomID" )

        self.USER_ENTITY_PROTOCOL_ATOM_ATTRIBUTES_TABLE = TableDB( table_name = "userEntityUnificationProtocolAtomAttribute",
                                                                     table_fields = [ FieldDB( field_name = "unificationAtomID",
                                                                                               data_type = "integer(2) unsigned",
                                                                                               null = False,
                                                                                               foreign_key = (self.USER_ENTITY_PROTOCOL_ATOMS_TABLE,
                                                                                                              self.USER_ENTITY_PROTOCOL_ATOMS_TABLE.get_field( field_name = "unificationAtomID" ) )),
                                                                                      FieldDB( field_name = "cross_referenced_code",
                                                                                               data_type = self.enum_eEattr_col,
                                                                                               null = False ) ],
                                                                     primary_key = ("unificationAtomID","cross_referenced_code") )


        ##################################################
        ###          ONTOLOGIES RELATED                ###
        ##################################################

        
        self.ONTOLOGY_IS_A_TABLE = TableDB( table_name = "ExternalEntityOntology_isA",
                                       table_fields = [ self.externalEntityID_field,
                                                        FieldDB( field_name = "is_a",
                                                                 data_type = self.externalEntityID_col_type,
                                                                 null = False )],
                                       indices = [(self.externalEntityID_col,),("is_a",)] )

        self.EXTENDED_ONTOLOGY_HIERARCHY_TABLE = TableDB( table_name = "ExtendedOntologyIsA",
                                                     table_fields = [ FieldDB( field_name = "parentExternalEntityID",
                                                                               data_type = self.externalEntityID_col_type,
                                                                               null = False ),
                                                                      FieldDB( field_name = "childExternalEntityID",
                                                                               data_type = self.externalEntityID_col_type,
                                                                               null = False ) ],
                                                     indices = [("parentExternalEntityID","childExternalEntityID")] )

        self.ONTOLOGY_IS_PART_OF_TABLE = TableDB( table_name = "ExternalEntityOntology_isPartOf",
                                             table_fields = [ self.externalEntityID_field,
                                                              FieldDB( field_name = "is_part_of",
                                                                       data_type = self.externalEntityID_col_type,
                                                                       null = False )],
                                             indices = [(self.externalEntityID_col,),("is_part_of",)] )

        self.ONTOLOGY_INFO_TABLE = TableDB( table_name = "ExternalEntityOntology",
                                       table_fields = [ self.externalEntityID_field,
                                                        FieldDB( field_name = "name",
                                                                 data_type = "varchar(255)",
                                                                 null = False ),
                                                        FieldDB( field_name = "linked_attribute",        # Attribute used in the external database as the main key
                                                                 data_type = "varchar(255)",
                                                                 null = False ),
                                                        FieldDB( field_name = "key_id",
                                                                 data_type = "integer(1) unsigned",
                                                                 null = False ),
                                                        FieldDB( field_name = "level_attribute",
                                                                 data_type = "varchar(255)",
                                                                 null = True ),
                                                        FieldDB( field_name = "description_attribute",   # Attribute used to show the tree in a user friendly way
                                                                 data_type = "varchar(255)",
                                                                 null = False )],
                                       primary_key = ["name"] )

        # TO CHECK: Al crear sense primary_keys, deixa de ser primary...    



        #####################################################
        ###   EXTERNAL ATTRIBUTES DESCRIPTION DATABASES   ###
        #####################################################                                                  

        self.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES = { "pdb": TableDB( table_name = "pdb",
                                                                   table_fields = [ FieldDB( field_name = "pdb",
                                                                                             data_type = "char(4)",
                                                                                             null = False ),
                                                                                    FieldDB( field_name = "chain",
                                                                                             data_type = "char(1)"),
                                                                                    FieldDB( field_name = "resolution",
                                                                                             data_type = "float"),
                                                                                    self.externalEntityID_field,
                                                                                    FieldDB( field_name = "total_num_atoms",     # Total number of atoms
                                                                                             data_type = "integer(3) unsigned",
                                                                                             null = False ),
                                                                                    FieldDB( field_name = "residue_num_list",    # Number of the residue (integer list)
                                                                                         data_type = "blob"),
                                                                                FieldDB( field_name = "residue_type_list",   # List of residue type
                                                                                         data_type = "text"),                # List of the types of the atoms (string list)
                                                                                FieldDB( field_name = "residue_atom_num_list",  #Number of atoms/residue (integer list)
                                                                                         data_type = "blob"),
                                                                                FieldDB( field_name = "atom_type_list",      # List of atom type (1character/atom)
                                                                                         data_type = "text"),
                                                                                FieldDB( field_name = "atom_name_list",      # List of the atom names (string list)
                                                                                         data_type = "text"),
                                                                                FieldDB( field_name = "atom_coordinates", # List of coordinates for all atoms (float list: x,y,z)
                                                                                         data_type = "mediumblob"),
                                                                                FieldDB( field_name = "hssp_residue_num_correspondences", # List of residue correspondences between hssp and pdb
                                                                                         data_type = "blob"),
                                                                                FieldDB( field_name = "residue_dssp_results",     # Results of the dssp program
                                                                                         data_type = "text"),
                                                                                FieldDB( field_name = "residue_hssp_entropy",     # Results of the entropy found in hssp
                                                                                         data_type = "blob"),
                                                                                FieldDB( field_name = "residue_hssp_norm_entropy", # Normalized entropy found in hssp
                                                                                         data_type = "blob"),
                                                                                FieldDB( field_name = "residue_hssp_variability", # Results of the variability found in hssp
                                                                                         data_type = "blob"),
                                                                                FieldDB( field_name = "conservation_hssp",      # Conservation found in the hssp
                                                                                         data_type = "blob"),
                                                                                FieldDB( field_name = "solvent_exposure_hssp",  # Solvent exposure found in hssp
                                                                                         data_type = "blob")],
                                                               primary_key = ("pdb","chain") ),

                                               "proteinSequence": TableDB( table_name = "sequenceProtein",
                                                                           table_fields = [ FieldDB( field_name = "proteinSequenceID",
                                                                                                     #data_type = "integer unsigned auto_increment"),
                                                                                                     data_type = "integer(4) unsigned"),
                                                                                                     #null = False),
                                                                                            FieldDB( field_name = "sequenceMD5",
                                                                                                     #data_type = "char(40) collate latin1_general_cs",
                                                                                                     data_type = "binary(16)",
                                                                                                     user_friendly_name = "sequenceMD5",
                                                                                                     optimize_space = 4,
                                                                                                     null = False),
                                                                                            FieldDB( field_name = "sequence",
                                                                                                     user_friendly_name = "sequence",
                                                                                                     data_type = "mediumtext",
                                                                                                     compress = 1,
                                                                                                     null = False),
                                                                                            FieldDB( field_name = "sequenceLength",
                                                                                                     data_type = "integer(2) unsigned",
                                                                                                     user_friendly_name = "length",
                                                                                                     null = False),
                                                                                            FieldDB( field_name = "sequenceMW",
                                                                                                     data_type = "float",
                                                                                                     user_friendly_name = "mw",
                                                                                                     null = False,
                                                                                                     default_value = 0),
                                                                                            FieldDB( field_name = "sequenceIP",
                                                                                                     data_type = "float",
                                                                                                     user_friendly_name = "ip",
                                                                                                     null = False,
                                                                                                     default_value = 0 ) ],
                                                                           primary_key = "sequenceMD5",
                                                                           indices = [("proteinSequenceID")]),

                                               "nucleotideSequence": TableDB( table_name = "sequenceNucleotide",
                                                                              table_fields = [ FieldDB( field_name = "nucleotideSequenceID",
                                                                                                        data_type = "integer(4) unsigned"),
                                                                                                        #ata_type = "integer(4) unsigned auto_increment"),
                                                                                                        #null = False),
                                                                                               FieldDB( field_name = "sequenceMD5",
                                                                                                        #data_type = "char(40) collate latin1_general_cs",
                                                                                                        data_type = "binary(16)",
                                                                                                        user_friendly_name = "sequenceMD5",
                                                                                                        optimize_space = 4,
                                                                                                        null = False),
                                                                                               FieldDB( field_name = "sequenceLength",
                                                                                                        data_type = "integer(2) unsigned",
                                                                                                        user_friendly_name = "length",
                                                                                                        null = False),
                                                                                               FieldDB( field_name = "sequence",
                                                                                                        data_type = "mediumtext",
                                                                                                        compress = 1,
                                                                                                        user_friendly_name = "sequence",
                                                                                                        null = False)],
                                                                              primary_key = "sequenceMD5",
                                                                              indices = [("nucleotideSequenceID")]) }

        self.EXTERNAL_ATTRIBUTE_DESCRIPTIONS_DICT = {}

        for actual_attribute in self.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES.keys():

            temp_fields = self.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES[actual_attribute].get_fields()

            self.EXTERNAL_ATTRIBUTE_DESCRIPTIONS_DICT[actual_attribute.lower()] = {"table": self.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES[actual_attribute],
                                                                              "fields": dict([ (x.get_user_friendly_name(),x) for x in temp_fields ]) }

         ##########################################
         ### SEQUENCE SIMILARITY RELATED TABLES ###
         ##########################################

        # Table to store CD-HIT clusters file
        self.CD_HIT_CLUSTERING_TABLE = TableDB( table_name = "sequenceProteinCD_HIT",
                                                table_fields = [ FieldDB( field_name = "representant_proteinSequenceID",
                                                                          data_type = "integer(4) unsigned",
                                                                          null = False ),
                                                        FieldDB( field_name = "proteinSequenceID",
                                                                 data_type = "integer(4) unsigned",
                                                                 null = False ),
                                                        FieldDB( field_name = "representant_start_position",
                                                                 data_type = "integer(2) unsigned",
                                                                 null = False ),
                                                        FieldDB( field_name = "representant_end_position",
                                                                 data_type = "integer(2) unsigned",
                                                                 null = False ),
                                                        FieldDB( field_name = "start_position",
                                                                 data_type = "integer(2) unsigned",
                                                                 null = False ),
                                                        FieldDB( field_name = "end_position",
                                                                 data_type = "integer(2) unsigned",
                                                                 null = False ),
                                                        FieldDB( field_name = "identity",
                                                                 data_type = "integer(1) unsigned",
                                                                 null = False )],
                                       indices = [("proteinSequenceID"),("representant_proteinSequenceID")] )



        self.PROTEIN_BLAST_RESULTS_TABLE = TableDB( table_name = "sequenceProteinBlastResults",
                                           table_fields = [ FieldDB( field_name = "sequenceID_A",
                                                                     data_type = "integer(4) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "sequenceID_B",
                                                                     data_type = "integer(4) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "evalue",
                                                                     data_type = "float",
                                                                     null = False ),
                                                            FieldDB( field_name = "score",
                                                                     data_type = "smallint(2) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "bit_score",
                                                                     data_type = "float unsigned", #float(4)
                                                                     null = False ),
                                                            FieldDB( field_name = "start_A",
                                                                     data_type = "smallint(2) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "end_A",
                                                                     data_type = "smallint(2) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "coverage_A",
                                                                     data_type = "tinyint(1) unsigned"),
                                                            FieldDB( field_name = "coverage_B",
                                                                     data_type = "tinyint(1) unsigned"),
                                                            FieldDB( field_name = "start_B",
                                                                     data_type = "smallint(2) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "end_B",
                                                                     data_type = "smallint(2) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "identities",
                                                                     data_type = "tinyint(1) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "similarity",
                                                                     data_type = "tinyint(1) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "gaps",
                                                                     data_type = "tinyint(1) unsigned",
                                                                     null = False ),
                                                            FieldDB( field_name = "program",
                                                                     data_type = "ENUM(\"bl2seq\",\"blastall\")"),
                                                            FieldDB( field_name = "filter",
                                                                     data_type = "ENUM(\"T\",\"F\")")],


                                           indices = [("sequenceID_A","identities","coverage_A"),
                                                      ("sequenceID_B","identities","coverage_B")])  # TO CHECK THE INDICES!!!


################
## PROCEDURES ##
################

### NOT USED AS IT IS VERY SLOW ###


REMOVE_SEQUENCE_DUPLICATES_PROCEDURE ="""
DELIMITER //
CREATE PROCEDURE remove_duplicates() 
       
  DETERMINISTIC 
  MODIFIES SQL DATA

BEGIN
 DECLARE sequenceMD5 BINARY(16);
 DECLARE previous_sMD5 BINARY(16) DEFAULT 0;
 DECLARE sequenceID int(10) unsigned;
 DECLARE sequence mediumtext;
 DECLARE sequenceLength int(2) unsigned;
 DECLARE sequenceMW float;
 DECLARE sequenceIP float;

 DECLARE done BOOLEAN DEFAULT FALSE;
 DECLARE cur1 CURSOR FOR SELECT * FROM sequenceProtein ORDER BY sequenceMD5;
 DECLARE CONTINUE HANDLER FOR SQLSTATE '02000' SET done = TRUE;

 OPEN cur1;

 cursor_loop: LOOP
   FETCH cur1 INTO sequenceID,sequenceMD5,sequence,sequenceLength,sequenceMW,sequenceIP;
   IF done THEN LEAVE cursor_loop; END IF;
   IF sequenceMD5!=previous_sMD5 THEN
      SET previous_sMD5 = sequenceMD5;
      INSERT INTO test_nonr_seq VALUES (sequenceID,sequenceMD5,sequence,sequenceLength,sequenceMW,sequenceIP);
   END IF;
 END LOOP cursor_loop;

 DEALLOCATE PREPARE stmt_query;
 CLOSE cur1;

END
//
delimiter ;
"""

