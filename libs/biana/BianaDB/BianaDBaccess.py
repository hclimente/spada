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


# Python library
import time
import sys
import md5
import traceback
import copy
from math import ceil
import re

# General Database
import ConnectorDB
import database
import biana.utilities.int_ascii as int_ascii

# Biana specific
import biana.BianaObjects as BianaObjects
from BianaDatabase import BianaDatabase

# Enable debugging for web php scripts
debug_web = 0 #1; #!


# a hard-coded version id assigned to source to prevent a previously created biana database <-> source code inconsistencies
#BIANA_SOURCE_CODE_VERSION = "Mar_16_09"   #Dynamic attributes biana database specific
BIANA_SOURCE_CODE_VERSION = "July_22_09"   #STRING score subcategory seperation, versionable IPI

class BianaDBaccess(object):
    """
    Class used as an interface with database biana
    """

    def __init__(self, dbname=None, dbhost=None, dbuser=None, dbpassword=None, dbport=None, dbsocket=None, use_buffer=False, lock_tables = False, check_integrity=False ):
        """
        "dbname" is the database name to which you want to connect to (required)
        "dbhost" is the machine with the mysql server that holds the biana database (required)
        "dbuser" is the mysql user (not required in most systems)
        "dbpassword" is the mysql password (not required in most systems)
        "dbport" is the mysql port (not required in most systems)

        The following parameters should not be used by standard users, only for advanced users or developers:
        "use_buffer" must be set True when populating database due to performance issues. Automatically controlled by parsers
        "lock_tables" Allow Connector db to lock tables when using them. It is set to True when populating database due to performance issues
        "check_integrity" determines if integrity of the database must be checked (if there is any parser not finished or if some table definitions have been changed) 
        """

        # opening connection to database biana using class BianaDB
        self.db = ConnectorDB.DB(dbname=dbname, dbhost=dbhost, dbuser=dbuser, dbpassword=dbpassword, dbport=dbport, dbsocket=dbsocket, lock_tables=lock_tables)
        self.dbname = dbname
        self.dbhost = dbhost
        self.dbuser=dbuser
        self.dbpassword=dbpassword
        self.dbport=dbport
        self.dbsocket=dbsocket
        self.lock_tables=lock_tables
        self.use_buffer = use_buffer

        # Check into BIANA database which database sources are available
        self.validSources = None
        self.valid_source_database_ids = {}
        self.available_unification_protocols = None  # Key: Description. Value: ID

        # Ontology related variables
        self.available_ontology_names = None   # Stores the available ontologies
        self.ontology_linked_attributes = None
        self.loaded_ontologies = {}

        self.db_version_modified = None # Variable used to control if it is necessary to update the database version identifier

        self.db_versions_list = None

        self.db.add_autoincrement_columns( table = "externalEntity", attribute = "externalEntityID" )
        self.db.add_autoincrement_columns( table = "externalEntityRelationParticipant", attribute = "externalEntityRelationParticipantID" )
        self.db.add_autoincrement_columns( table = "sequenceProtein", attribute = "proteinSequenceID" )
        self.db.add_autoincrement_columns( table = "sequenceNucleotide", attribute = "nucleotideSequenceID" )
        self.db.add_autoincrement_columns( table = "externalDatabaseAttributeTransfer", attribute = "keyID" )

        self.db_description = None
        self.db_optimized_for = None

        self.biana_database = BianaDatabase()

        self.transferred_attributes = {} # Dictionary to store how transfer attributes between databases. Key: transferred attribute
                                         #                                                                Value: set of (externalDatabaseID, key_attribute)
                                         # This information has to be hardcoded in parsers
                                         # For example, pfam code can be tranferred from pfam database using as key attribute the uniprot accession
                                         #              This means that pfam codes are going to be fetched from pfam database
                                         # Another example is taxonomy name. This can be transferred using as key attribute taxid
        self.key_attribute_ids = {}


	self.versionable_external_entity_identifier_attributes = set() # to decide stripping versioning parts of the values in BianaSessionManager - create set from identifier

        if self.dbname is not None:

            if not self.db.check_consistency_with_given_source_version(BIANA_SOURCE_CODE_VERSION):
                # sys.stderr.write("Bailing.. The database you are trying to use is inconsistent with the current source code!\n")
                raise Exception("Source code - database inconsistency")

            self._load_biana_types_and_attributes()

            # Create new tables if necessary, and modify existing ones if necessary
            self.db.check_database(database = self.biana_database, verbose=True)

            if check_integrity:
                self._check_database_integrity()

            self._load_biana_database_information()  #JAVI: Indented

        # Temporal data is used while parsing for storing temporal data
        self.temporal_data = {}
        self.temporal_data["relations_hierarchy_parents"] = {}  # stores all the parents for each external entity id
        self.store_relations_hierarchy = False

        return


    # ----
    # methods required for using pickle with biana objects
    # ----
    def __getstate__(self):
        odict = self.__dict__.copy() # copy the dict since we change it
        del odict['db']              # remove database entry
        return odict

    def __setstate__(self, dict):

        self.__dict__.update(dict) # update attributes
        self.db = ConnectorDB.DB(dict["dbname"], dict["dbhost"], dict["dbuser"], dict["dbpassword"], dict["dbport"], dict["dbsocket"], dict["lock_tables"])
        try:
            self.db.add_autoincrement_columns( table = "externalEntity", attribute = "externalEntityID" )
            self.db.add_autoincrement_columns( table = "externalEntityRelationParticipant", attribute = "externalEntityRelationParticipantID" )
            self.db.add_autoincrement_columns( table = "sequenceProtein", attribute = "proteinSequenceID" )
            self.db.add_autoincrement_columns( table = "sequenceNucleotide", attribute = "nucleotideSequenceID" )
            self.db.add_autoincrement_columns( table = "externalEntityEquivalencesTransfer", attribute = "keyID" )
        except:
            pass
        return
        
    def __getnewargs__(self):
        return (self.dbname, self.dbhost, self.dbuser, self.dbpassword)

    def __str__(self):
        return "BIANA Connection. Database connection: %s" %self.db

    def isOptimizedForRunning(self):
        return self.db_optimized_for=="running"

    def reconnect(self):
        self.db = ConnectorDB.DB(dbname=self.dbname, dbhost=self.dbhost, dbuser=self.dbuser, dbpassword=self.dbpassword, dbport=self.dbport, dbsocket=self.dbsocket, lock_tables=self.lock_tables)

    def create_database(self, dbname, description="BIANA DATABASE", optimize_for="parsing", ignore_primary_keys=False):
        """
        It creates necessary tables into a database server

        "description" is a tag for label the new database

        "optimize_for" can take two distinct values: "parsing" and "running". By default, it creates a database optimized for parsing (as it will be created empty). 
        """

        self.db.insert_db_content( sql_query = self.biana_database.get_sql_query(ignore_primary_keys=ignore_primary_keys).split(";") )

        self.optimize_database_for( mode=optimize_for )
        
        # Insert initial valid data types and attributes
        import biana.biana_globals as BIANA_GLOBALS

        [ self.add_valid_external_entity_type(current_type) for current_type in BIANA_GLOBALS.EXTERNAL_ENTITY_TYPES ]
        [ self.add_valid_external_entity_relation_type(current_type) for current_type in BIANA_GLOBALS.EXTERNAL_ENTITY_RELATION_TYPES ]
        [ self.add_valid_identifier_reference_types(current_reference_type) for current_reference_type in BIANA_GLOBALS.VALID_IDENTIFIER_REFERENCE_TYPES ]
        [ self.add_valid_external_entity_attribute_type(current_type, data_type, "eE attribute") for (current_type, data_type) in BIANA_GLOBALS.EXTERNAL_ENTITY_GENERAL_ATTRIBUTES ]
        [ self.add_valid_external_entity_attribute_type(current_type, data_type, "eE identifier attribute") for (current_type, data_type) in BIANA_GLOBALS.EXTERNAL_ENTITY_IDENTIFIER_ATTRIBUTES ]
        [ self.add_valid_external_entity_attribute_type(current_type, data_type, "eE versionable identifier attribute") for (current_type, data_type) in BIANA_GLOBALS.EXTERNAL_ENTITY_VERSIONABLE_IDENTIFIER_ATTRIBUTE_TYPES ]
        [ self.add_valid_external_entity_attribute_type(current_type, data_type, "eE descriptive searchable attribute") for (current_type, data_type) in BIANA_GLOBALS.EXTERNAL_ENTITY_DESCRIPTIVE_SEARCHABLE_ATTRIBUTE_TYPES ]
        [ self.add_valid_external_entity_attribute_type(current_type, data_type, "eE descriptive attribute") for (current_type, data_type) in BIANA_GLOBALS.EXTERNAL_ENTITY_DESCRIPTIVE_ATTRIBUTE_TYPES ]
        [ self.add_valid_external_entity_attribute_type(current_type, data_type, "eE numeric attribute") for (current_type, data_type) in BIANA_GLOBALS.EXTERNAL_ENTITY_NUMERIC_ATTRIBUTE_TYPES ]
        [ self.add_valid_external_entity_relation_participant_attribute_type(current_type, data_type) for (current_type, data_type) in BIANA_GLOBALS.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TYPES ]
        [ self.add_valid_external_entity_attribute_type(current_type, data_type_dict, "eE special attribute") for (current_type, data_type_dict) in BIANA_GLOBALS.EXTERNAL_ENTITY_SPECIAL_ATTRIBUTE_TYPES.iteritems()]
        

        self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.BIANA_DATABASE_TABLE,
                                                                              column_values = (("description",description),
                                                                                               ("optimized_for",optimize_for)),
                                                                              use_buffer = False ) )
        self.db.close()

        self.db = ConnectorDB.DB(dbname=dbname, dbhost=self.dbhost, dbuser=self.dbuser, dbpassword=self.dbpassword, dbport=self.dbport, dbsocket=self.dbsocket)
        
        self.dbname = dbname

        self._load_biana_types_and_attributes()

        self.db.check_database(database=self.biana_database, ignore_primary_keys=ignore_primary_keys)

        # Creates a default protocol (used to work without any kind of unification)
        default_protocol = BianaObjects.UnificationProtocol( description = "No unification", BianaDatabaseVersion = "X" )
        
        self.db_description = description
        self.db_optimized_for = optimize_for

        self._create_new_unification_protocol_tables(default_protocol)


    def refresh_database_information(self):

        self._load_biana_types_and_attributes()
        self._load_biana_database_information()
        self.db.check_database(self.biana_database)

          
    def add_valid_external_entity_type(self, type):
        
        if self.biana_database.is_valid_external_entity_type(type):
            sys.stderr.write("Trying to add an existing External Entity Type: %s. Not added again.\n" %type)
            return
        
        self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.TYPES_AND_ATTRIBUTES_TABLE,
                                                                              column_values = [("name",type),
                                                                                               ("category","eE type")],

                                                                              use_buffer = False ) )

        self.biana_database.add_valid_external_entity_type(type)


    def get_valid_external_entity_types(self):

        return self.biana_database.get_valid_external_entity_types()


    def add_valid_external_entity_relation_type(self, type): 

        
        if self.biana_database.is_valid_external_entity_relation_type(type):
            sys.stderr.write("Trying to add an existing External Entity Relation Type: %s. Not added again.\n" %type)
            return

        
        self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.TYPES_AND_ATTRIBUTES_TABLE,
                                                                              column_values = [("name",type),
                                                                                               ("category","eEr type")],
                                                                              use_buffer = False ) )

        self.biana_database.add_valid_external_entity_relation_type(type)


    def get_valid_external_entity_relation_types(self):

        return self.biana_database.get_valid_external_entity_relation_types()

        
    def add_valid_identifier_reference_types(self, current_reference_type):
        
        if self.biana_database.is_valid_identifier_reference_type(current_reference_type):
            sys.stderr.write("Trying to add an existing Identifier reference type: %s. Not added again.\n" %tcurrent_reference_ype)
            return

        self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.TYPES_AND_ATTRIBUTES_TABLE,
                                                                              column_values = [("name",current_reference_type),
                                                                                               ("category","identifier reference type")],
                                                                              use_buffer = False ) )

        self.biana_database.add_valid_identifier_reference_type( current_reference_type )

        

    def add_valid_external_entity_attribute_type(self, name, data_type, category):

        if self.biana_database.is_valid_external_entity_attribute_type(name):
            sys.stderr.write("Trying to add an existing External Entity Attribute Type: %s. Not added again\n" %name)
            return

        additional_fields_tuple_list = []

        if category.lower()=="ee special attribute":
            data_type_dict = data_type
            for current_field_tuple in data_type_dict["fields"]:
                current_data_type = current_field_tuple[1]
                current_field = current_field_tuple[0]
                if current_field_tuple[0] == "value":
                    self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.TYPES_AND_ATTRIBUTES_TABLE,
                                                                                          column_values = [("name",name),
                                                                                                           ("data_type",current_data_type),
                                                                                                           ("category",category)],
                                                                                          use_buffer = False ) )
                else:
                    current_null = current_field_tuple[2]
                    if current_null is False: current_null=0
                    else:                     current_null=1 

                    additional_fields_tuple_list.append((current_field, current_data_type, current_null))

                    self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.SPECIAL_ATTRIBUTES_TABLE,
                                                                                          column_values = [("attribute_name",name),
                                                                                                           ("field_name",current_field),
                                                                                                           ("data_type",current_data_type),
                                                                                                           ("canbenull",current_null)],
                                                                                          use_buffer = False ) )
                
        else:
            self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.TYPES_AND_ATTRIBUTES_TABLE,
                                                                                  column_values = [("name",name),
                                                                                                   ("data_type",data_type),
                                                                                                   ("category",category)],
                                                                                  use_buffer = False ) )

        self.biana_database.add_valid_external_entity_attribute_type(name, data_type, category, additional_fields_tuple_list)


    def add_valid_external_entity_relation_participant_attribute_type(self, name, data_type):

        if self.biana_database.is_valid_external_entity_relation_participant_attribute_type(name):
            sys.stderr.write("Trying to add an existing External Entity Relation Participant Attribute Type: %s. Not added again\n" %name)
            return
            
        self.db.insert_db_content( sql_query = self.db._get_insert_sql_query( table = self.biana_database.TYPES_AND_ATTRIBUTES_TABLE,
                                                                              column_values = [("name",name),
                                                                                               ("data_type",data_type),
                                                                                               ("category","eErP attribute")],
                                                                              use_buffer = False ) )


    def _load_biana_types_and_attributes(self):
        """
        Method to load current biana database types and attributes
        """

        data = self.db.select_db_content( self.db._get_select_sql_query( tables = [ self.biana_database.TYPES_AND_ATTRIBUTES_TABLE ],
                                                                         columns = ["name", "data_type","category"] ),
                                          answer_mode = "raw" )

        for current_data in data:
            if current_data[2].lower() == "identifier reference type":  self.biana_database.add_valid_identifier_reference_type(current_data[0])
            

        for current_data in data:
            current_category = current_data[2]
            current_attribute = current_data[0]
            data_type = current_data[1]
            new_table = None

            if current_category.lower() == "ee type":                  self.biana_database.add_valid_external_entity_type(current_attribute)                
            elif current_category.lower() == "eer type":               self.biana_database.add_valid_external_entity_relation_type(current_attribute)
            elif current_category.lower() == "eerp attribute" or current_category.lower() == "eerp numeric attribute" or current_category.lower() == "eerp descriptive attribute" or  current_category.lower() == "eerp descriptive searchable attribute":
                self.biana_database.add_valid_external_entity_relation_participant_attribute_type( current_attribute, data_type, current_category, [])
            elif current_category.lower() == "ee identifier attribute" or current_category.lower() == "ee versionable identifier attribute" or current_category.lower() == "ee descriptive attribute" or current_category.lower() == "ee descriptive searchable attribute" or current_category.lower() == "ee numeric attribute":          
		self.biana_database.add_valid_external_entity_attribute_type( current_attribute, data_type, current_category, [] )
		if current_category.lower() == "ee versionable identifier attribute":
		    self.versionable_external_entity_identifier_attributes.add(current_attribute.lower())
            elif current_category.lower() == "ee attribute":           self.biana_database.add_valid_external_entity_attribute_type( current_attribute, data_type, current_category, [] )
            elif current_category.lower()== "ee special attribute":
                # Get additional fields
                special_data = self.db.select_db_content( self.db._get_select_sql_query( tables = [ self.biana_database.SPECIAL_ATTRIBUTES_TABLE ],
                                                                                         columns = ["field_name", "data_type", "canbenull"],
                                                                                         fixed_conditions = [("attribute_name","=",current_attribute)] ),
                                                          answer_mode = "raw" )
                self.biana_database.add_valid_external_entity_attribute_type( current_attribute, data_type, current_category, special_data )
            elif current_category.lower()== "identifier reference type": pass
            else:
                raise ValueError("%s category not recognized!" %current_category)

            if new_table:
                self.biana_database.add_table(new_table)

    def get_versionable_external_entity_identifier_attributes(self):
	return self.versionable_external_entity_identifier_attributes

    def _transform_attribute_value_data_type_to_biana_database_attribute_data_type( self, attribute_identifier, value ):

        NUMBER_RE = re.compile("[0-9]+")
        CHAR_RE = re.compile("(char|varchar|text)\((\d+)\)") # (var|)char

        new_value = value

        data_type = self.biana_database.get_attribute_data_type(attribute_identifier)

        # to avoid incorrect string value assignment to a integer field in database
        if isinstance(value, str):
            if value.strip() == "":
                sys.stderr.write("Trying to create an External Entity Attribute with an empty value")
                new_value = None

            if data_type.lower().startswith("int"):
                search = NUMBER_RE.search(value)
                if search:
                    new_value = value[search.start():]
                else:
                    new_value = None
                    sys.stderr.write("Trying to create an External Entity Attribute %s of an incorrect value: %s\n" % (attribute_identifier, value))
            else:
                search = CHAR_RE.match(value)
                char_count = 0

                if search:
                    char_count = int(search.group(2))
                else:
                    char_count = 200000000

                if char_count > 0 and char_count < len(value):
                    new_value = value[:char_count]
                    sys.stderr.write("Trying to create an External Entity Attribute %s of an overfloating value: %s\n" % (attribute_identifier, value))

        return new_value


    def _load_biana_database_information(self):
        """
        Method to load biana database information into this class

        It loads information about the database description, its optimization status, special attributes and special attributes

        It is used when initializing a BianaDBaccess object
        """

        if self.dbname is not None:

            data = self.db.select_db_content( sql_query = self.db._get_select_sql_query( tables = [ self.biana_database.BIANA_DATABASE_TABLE ],
                                                                                         columns = ["description","optimized_for"] ),
                                              answer_mode = "raw" )
            
            self.db_description = data[0][0]
            self.db_optimized_for = data[0][1]


            data = self.db.select_db_content( sql_query = self.db._get_select_sql_query( tables = [ self.biana_database.EXTERNAL_DATABASE_ATTRIBUTE_TRANSFER_TABLE ],
                                                                                         columns = ["keyID","externalDatabaseID", "attributeKey", "transferAttribute"] ),
                                              answer_mode = "raw" )

            for (keyID, externalDatabaseID, key, transfer) in data:
                self.transferred_attributes.setdefault(transfer.lower(),set()).add((externalDatabaseID,key.lower(),keyID))
                self.key_attribute_ids.setdefault((externalDatabaseID,key.lower()),keyID)


    def get_available_ontology_names(self, name=None):
        """
        """        

        if self.available_ontology_names is None:

            data = self.db.select_db_content( sql_query = self.db._get_select_sql_query ( tables = [ self.biana_database.ONTOLOGY_INFO_TABLE ],
                                                                                          columns = ["name","linked_attribute","key_id","level_attribute"] ),
                                              answer_mode = "raw" )

            self.available_ontology_names = dict([ (x[0],x[1]) for x in data ])
            self.ontology_linked_attributes = dict([ (x[1],{"ontology_name":x[0], "key_id":x[2], "level_attribute":x[3]}) for x in data ])

        if name is None:
            return self.available_ontology_names
        else:
            return self.available_ontology_names[name.lower()]

    
    def _is_ontology_linked_attribute(self, attribute_identifier):
        
        if self.ontology_linked_attributes is None:
            self.get_available_ontology_names()

        if( attribute_identifier.lower() in self.ontology_linked_attributes ):
            return True

        return False


    def _add_key_attribute(self, external_database_id, key_attribute):
        """
        Adds a key attribute, used in transfer attributes (key attribute) and ontologies (linked attribute)

        Creates the necessary tables
        """
        key_attribute = key_attribute.lower()
        
        if self.key_attribute_ids.has_key( (external_database_id, key_attribute) ):

            return self.key_attribute_ids[ (external_database_id, key_attribute) ]

        else:

            new_id = self._get_new_key_id()

            value_data_type = self.biana_database.get_attribute_data_type(key_attribute)
            
            new_table = database.TableDB( table_name = self._get_key_attribute_table_name( key_id = new_id ),
                                          table_fields = [ database.FieldDB(field_name = "externalEntityID",
                                                                            data_type = "integer(4) unsigned",
                                                                            null = False ),
                                                           database.FieldDB(field_name = "value",
                                                                            data_type =  value_data_type,
                                                                            null = False ) ],
                                          indices = [("value","externalEntityID")] )
            
            self.db.insert_db_content( sql_query = new_table.create_mysql_query(), answer_mode = None )

            self.key_attribute_ids.setdefault((external_database_id,key_attribute.lower()),new_id)
            
            return new_id
        
        


    def _add_transfer_attribute(self, externalDatabaseID, key_attribute, transfer_attribute):

        key_id = self._add_key_attribute( external_database_id = externalDatabaseID,
                                          key_attribute = key_attribute )

        self.db.insert_db_content( sql_query = self.db._get_insert_sql_query(table = self.biana_database.EXTERNAL_DATABASE_ATTRIBUTE_TRANSFER_TABLE.get_table_name(),
                                                                             column_values = (("keyID",key_id),
                                                                                              ("externalDatabaseID", externalDatabaseID),
                                                                                              ("attributeKey", key_attribute),
                                                                                              ("transferAttribute", transfer_attribute)),
                                                                             use_buffer = False ),
                                   answer_mode = None )

        self.transferred_attributes.setdefault(transfer_attribute.lower(),set()).add((externalDatabaseID,key_attribute.lower()))
        
        

        # Create a new variable for this class. This variable is only created in this method when inserting a new transferred attribute
        try:
            self.transferID_key[transfer_attribute.lower()] = (key_id,key_attribute.lower())
        except:
            self.transferID_key = {}
            self.transferID_key[transfer_attribute.lower()] = (key_id,key_attribute.lower())

    
    def _is_transferred_attribute(self, attribute_identifier):
        return self.transferred_attributes.has_key(attribute_identifier.lower())

    def _is_key_attribute(self, attribute_identifier, external_database_id):
        return self.key_attribute_ids.has_key( (external_database_id, attribute_identifier.lower()) )

    def _get_key_attribute_table_name( self, key_id ):
        return "key_attribute_%s" %key_id

    def _get_linked_attribute_ontology_name(self, attribute_identifier):
        """
        Returns the name of the ontology linked to an attribute
        """

        if self._is_ontology_linked_attribute(attribute_identifier):
            return self.ontology_linked_attributes[attribute_identifier.lower()]["name"]
            raise ValueError("Trying to get the name for an unlinked attribute")
        else:
            raise ValueError("Trying to get the name for an unlinked attribute")



    def optimize_database_for(self, mode, optimize=False):
        """
        Performs some modifications in database in order to optimize database access efficiency for parsing or running
        
        "mode" can take two different values: "running" or "parsing"

        "optimize" parameter will be only used if mode is "running"
        """

        if mode == "parsing":
            self.db._disable_indices()

        elif mode == "running":
            self.db._enable_indices()
            self.db.insert_db_content( sql_query = self.biana_database.optimize_database( optimize = optimize, analyze = True ) )
            
            # NOT DONE AS IT IS VERY INEFFICIENT
            # Update the external entity equivalences for attribute transfer
            #print "Updating attribute crossing table..."
            #self._update_external_entity_equivalences_for_attribute_transfer()


            # Precalculate relations hierarchy
            self._update_relations_hierarchy()

            
            

        
        else:
            raise ValueError("optimizing database mode can only be \"parsing\" or \"running\"")

        # Stores the 
        self.db_optimized_for = mode

        self.db.insert_db_content( sql_query = self.db._get_update_sql_query(table = self.biana_database.BIANA_DATABASE_TABLE,
                                                                             update_column_values = [("optimized_for", self.db_optimized_for)]) )
        

    def close(self):
        """
        Close the connection with database

        Although it is only necessary to close the BianaDBaccess object when parsing,
        is highly recommended to use it always as maybe in the future this requirite may change
        """

        #print "Closing BianaDBaccess"
        
        # Check if there is temporal data to be processed in temporal buffer
        # Process relations hierarchy temporal data
        if len( self.temporal_data["relations_hierarchy_parents"] )>0:
            hierarchy_set = set()

            eE_eERid_dict = self.temporal_data["relations_hierarchy_parents"]

            def get_all_parents(eEid):
                parents = set()
                if( eE_eERid_dict.has_key(eEid) ):
                    parents.update(eE_eERid_dict[eEid])
                    for current_eEid in eE_eERid_dict[eEid]:
                        parents.update(get_all_parents(current_eEid))
                return parents

            for eEid in self.temporal_data["relations_hierarchy_parents"].keys():
                for current_parent in get_all_parents(eEid):
                    hierarchy_set.add((current_parent,eEid))


            for current_eERid, current_eEid in hierarchy_set:
                self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE,
                                                                          column_values = (("externalEntityRelationParticipantID", self._get_new_external_entity_relation_participant_id()), 
                                                                                           (self.biana_database.external_entity_relation_id_col, current_eERid),
                                                                                           (self.biana_database.externalEntityID_col, current_eEid )),
                                                                          use_buffer = True ))
        # If database has been modified, add the control id
        if self.db_version_modified:
            self._update_bianaDB_autoincrement_fields()
            self._update_bianaDB_version()

        self.db.close()
        



    ####################################################################################
    #  METHODS USED TO MANAGE AVAILABLE PROTEIN TYPES AND DATABASES IN BIANA DATABASE  #
    ####################################################################################

    # --------------------------------------------------------------------------
    # Methods used to control BIANA database autoincrement fields
    # --------------------------------------------------------------------------

    def _update_bianaDB_autoincrement_fields(self):
        """
        Method used internally to control biana autoincrement fields.
        Called by _update_bianaDB_version since updation of autoincrement fields is required before updating version stable ids.
        """
        column_values = []
        #if self._get_last_external_entity_id() is not None:
        column_values.append(("last_externalEntityID",self._get_last_external_entity_id()))
        #if self._get_last_external_entity_relation_participant_id() is not None:
        column_values.append(("last_externalEntityRelationParticipantID",self._get_last_external_entity_relation_participant_id()))
        #if self.get_last_sequence_protein_id() is not None:
        column_values.append(("last_proteinSequenceID", self._get_last_sequenceProtein_id()))
        #if self.get_last_sequence_nucleotide_id():
        column_values.append(("last_nucleotideSequenceID", self._get_last_sequenceNucleotide_id()))
        #if self.get_last_transfer_id():
        column_values.append(("last_keyID", self._get_last_key_id()))

        #if len(column_values) > 0:
        self.db.insert_db_content( self.db._get_update_sql_query(table = self.biana_database.BIANA_DATABASE_TABLE.get_table_name(),
                                                                 update_column_values = column_values),
                                   answer_mode = None )

    # --------------------------------------------------------------------------
    # Methods used to control BIANA database version
    # --------------------------------------------------------------------------
    def _update_bianaDB_version(self):
        """
        Method used internally to control biana version to be able to control compatibilities between diffent database versions
        
        It also stores the stable values for the database
        """

        date = time.localtime()
        actual_date = "%s-%s-%s" %(date[0],date[1],date[2])
        
        self.db.insert_db_content( self.db._get_insert_sql_query(table = self.biana_database.DATABASE_VERSION_TABLE.get_table_name(),
                                                                 column_values = (("dbControlID", md5.new(str(time.localtime())+str(time.time())).hexdigest()),
                                                                                  ("date", actual_date),
                                                                                  #("stable_externalEntityID",self._get_new_external_entity_id()-1),
                                                                                  ("stable_externalEntityID",self._get_last_external_entity_id()),
                                                                                  #("stable_externalEntityRelationParticipantID",self._get_new_external_entity_relation_participant_id()-1))),
                                                                                  ("stable_externalEntityRelationParticipantID",self._get_last_external_entity_relation_participant_id()))),
                                   answer_mode = None )


    def _check_database_integrity(self):
        """
        Checks the integrity of the database
        """


        sys.stderr.write("checking database integrity\n")

        self.db._check_locked_table(self.biana_database.EXTERNAL_DATABASE_TABLE.get_table_name())
        if len(self.db.select_db_content( sql_query = "SELECT * FROM %s WHERE parsingTime IS NULL" %self.biana_database.EXTERNAL_DATABASE_TABLE, answer_mode="list" )):
             sys.stderr.write("Database does not have integrity. Fixing it...")
             self._rollback()
	return


    def _rollback(self):
        """
        Method used to undo all the changes from an stable data insert
        """

        max_external_entity_id  = self._get_last_stable_external_entity_id()

        max_external_entity_relation_participant_id  = self._get_last_stable_external_entity_relation_participant_id()

        for current_table in self.biana_database.get_tables():
            if current_table.has_field("externalEntityID"):
                self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = current_table.get_table_name(),
                                                                                      fixed_conditions = [("externalEntityID",">",max_external_entity_id, None)] ),
                                           answer_mode = None )

            if current_table.has_field("externalEntityRelationParticipantID"):
                self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = current_table.get_table_name(),
                                                                                      fixed_conditions = [("externalEntityRelationParticipantID",">",max_external_entity_relation_participant_id, None)] ),
                                           answer_mode = None )


        # Delete keyID
        max_key_id = self._get_last_stable_key_id()
        transfer_tables_list = self.db.select_db_content( sql_query = "SHOW TABLES",
                                                          answer_mode = "list" )

        self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = self.biana_database.EXTERNAL_DATABASE_ATTRIBUTE_TRANSFER_TABLE.get_table_name(),
                                                                              fixed_conditions = [("keyID",">",max_key_id, None)] ),
                                   answer_mode = None )

        regex = re.compile("key_attribute_(\d+)")
        for current_table in transfer_tables_list:
            m = regex.match(current_table)
            if m:
                if( int(m.group(1))>int(max_key_id) ):
                    self.db.insert_db_content( sql_query = self.db._get_drop_sql_query( table_list = [current_table] ) )
        
        # Other tables to delete information
        # Ontology specific information
        for current_table in [self.biana_database.ONTOLOGY_IS_A_TABLE.get_table_name(), self.biana_database.ONTOLOGY_IS_PART_OF_TABLE.get_table_name() ]:
            self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = current_table,
                                                                                  fixed_conditions = [("externalEntityID",">",max_external_entity_id, None)] ),
                                       answer_mode = None )
        
        self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = self.biana_database.EXTENDED_ONTOLOGY_HIERARCHY_TABLE.get_table_name(),
                                                                              fixed_conditions = [("parentExternalEntityID",">",max_external_entity_id, None)] ),
                                   answer_mode = None )

        self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = self.biana_database.EXTENDED_ONTOLOGY_HIERARCHY_TABLE.get_table_name(),
                                                                              fixed_conditions = [("childExternalEntityID",">",max_external_entity_id,None)] ),
                                   answer_mode = None )


        # Remove sequences
        max_protein_sequence_id = self._get_last_stable_sequenceProtein_id()
        self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"].get_table_name(),
                                                                              fixed_conditions = [("proteinSequenceID",">",max_protein_sequence_id,None)] ),
                                   answer_mode = None )

        max_nucleotide_sequence_id = self._get_last_stable_sequenceNucleotide_id()
        self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"].get_table_name(),
                                                                              fixed_conditions = [("nucleotideSequenceID",">",max_nucleotide_sequence_id,None)] ),
                                   answer_mode = None )
            
        # DELETE ALL external Databases without parsing time, as it means the parsing has not been finished...
        self.db.insert_db_content( sql_query = self.db._get_delete_sql_query( table = self.biana_database.EXTERNAL_DATABASE_TABLE.get_table_name(),
                                                                              fixed_conditions = [("parsingTime","IS","NULL",None)] ),
                                   answer_mode = None )

        return


    def _get_db_versions_list(self):

        if self.db_versions_list is None:

            self.db_versions_list = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.DATABASE_VERSION_TABLE.get_table_name()],
                                                                                                     columns = ["dbControlID"]),
                                                               answer_mode = "list" )
        return self.db_versions_list
        

    def _get_table_names(self):

        return self.db.select_db_content( self.db.GetTableNames(), answer_mode="list" )



    def get_external_database(self, database_id):

        self._get_valid_source_dbs()

        return self.valid_source_database_ids[database_id]


    def get_external_database_list(self):
        """
        """
        
        self._get_valid_source_dbs()
        
        return self.valid_source_database_ids.values()

    def _get_valid_source_databases_by_id(self):
        
        self._get_valid_source_dbs()
        return self.valid_source_database_ids


    def _get_valid_source_dbs(self):
        """
        Returns a dictionary with current external databases stored in BIANA database

        Database objects are accessed by its name and version
        """

        if self.validSources is None:
   
            
            columns = ["externalDatabaseID","databaseName","databaseVersion","parsedFile","parsedDate","databaseDescription","defaultExternalEntityAttribute","isPromiscuous"]
            #columns.extend( list(available_attributes) )

            try:
                current_sources = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_DATABASE_TABLE.get_table_name()],
                                                                                            columns = columns ),
                                                             answer_mode="raw", remove_duplicates="no" )
                return_dict = {}
                for actual_source in current_sources:
                    databaseObject = BianaObjects.ExternalDatabase( databaseName = actual_source[1],
                                                                    databaseVersion = actual_source[2],
                                                                    databaseFile = actual_source[3],
                                                                    databaseDescription = actual_source[5],
                                                                    defaultExternalEntityAttribute = actual_source[6],
                                                                    databaseDate = actual_source[4],
                                                                    externalDatabaseID = actual_source[0],
								    isPromiscuous = actual_source[7])

                    self.valid_source_database_ids[ databaseObject.get_id() ] = databaseObject

                    available_eE_attributes = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eE_ATTRIBUTE_TABLE],
                                                                                                        columns = ["attributeType"],
                                                                                                        fixed_conditions = [("externalDatabaseID","=",databaseObject.get_id())] ),
                                                                      answer_mode = "list", remove_duplicates = "no" )
                    
                    for actual_attribute in available_eE_attributes:
                        databaseObject.add_valid_external_entity_attribute_type(actual_attribute)

                    
                    available_eEr_attributes = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eEr_ATTRIBUTE_TABLE],
                                                                                                         columns = ["attributeType"],
                                                                                                         fixed_conditions = [("externalDatabaseID","=",databaseObject.get_id())] ),
                                                                          answer_mode = "list", remove_duplicates = "no" )
                    
                    for actual_attribute in available_eEr_attributes:
                        databaseObject.add_valid_external_entity_relation_attribute_type(actual_attribute)

                    
                    available_eE_types = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eE_TYPES_TABLE],
                                                                                                   columns = ["eEType"],
                                                                                                   fixed_conditions = [("externalDatabaseID","=",databaseObject.get_id())] ),
                                                                    answer_mode = "list", remove_duplicates = "no" )
                    
                    for actual_attribute in available_eE_types:
                        databaseObject.add_valid_external_entity_type(actual_attribute)

                    
                    
	            available_eEr_types = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eEr_TYPES_TABLE],
                                                                                                   columns = ["eErType"],
                                                                                                   fixed_conditions = [("externalDatabaseID","=",databaseObject.get_id())] ),
                                                                    answer_mode = "list", remove_duplicates = "no" )
                    
                    for actual_attribute in available_eEr_types:
                        databaseObject.add_valid_external_entity_relation_type(actual_attribute)

                    

                        
                    if return_dict.has_key(actual_source[1]):
                        return_dict[actual_source[1]][actual_source[2]] = databaseObject
                    else:
                        return_dict[actual_source[1]] = {actual_source[2]: databaseObject}

            except:
            	raise
            
            self.validSources = return_dict

        return self.validSources

    
    ####################################################################################
    #  DATABASE MODIFICATION METHODS (INSERT EXTERNAL INFORMATION TO THE DATABASE      #
    ####################################################################################

    def insert_new_external_database( self, externalDatabase ):
        """
        Inserts into database the information of a new external database that is being integrated into BIANA database
        """
        # First: check if this databaseName is already in the database

        if self._get_valid_source_dbs().has_key(externalDatabase.get_name()) and self._get_valid_source_dbs()[externalDatabase.get_name()].has_key(externalDatabase.get_version()):
            sys.stderr.write("\n\nDatabase Name and Version already existed in the database.\n\n")
            sys.exit(-1)
        else:
            #content_type_list = externalDatabase.get_content_types()
            #for parent, child_list in self.biana_database.EXTERNAL_DATA_TYPE_HIERARCHY_DICTIONARY.iteritems():

            is_promiscuous = 0
            if externalDatabase.get_promiscuity() == True:
                is_promiscuous = 1
            
            new_database_identifier = self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_DATABASE_TABLE.get_table_name(),
                                                                                                column_values = ( ("databaseName", externalDatabase.get_name().strip("\"")),
                                                                                                                  ("databaseVersion", externalDatabase.get_version().strip("\"")),
                                                                                                                  ("parsedFile",externalDatabase.get_parsed_file().strip("\"") ),
                                                                                                                  ("parsedDate", externalDatabase.get_parsing_date().strip("\"") ),
                                                                                                                  ("defaultExternalEntityAttribute", externalDatabase.get_default_eE_attribute().strip("\"") ),
                                                                                                                  ("databaseDescription", externalDatabase.get_description().strip("\"") ),
                                                                                                                  ("isPromiscuous", is_promiscuous) ),
                                                                                                use_buffer=False ),
                                                                 answer_mode="last_id" )

            externalDatabase.set_id(externalDatabaseID = new_database_identifier)
                                                                 
            # Add this external database in the current external databases on memory
            if self._get_valid_source_dbs().has_key(externalDatabase.get_name()):
                self._get_valid_source_dbs()[externalDatabase.get_name()][externalDatabase.get_version()] = externalDatabase
            else:
                self._get_valid_source_dbs()[externalDatabase.get_name()] = {externalDatabase.get_version(): externalDatabase }
                
            self.valid_source_database_ids[new_database_identifier] = externalDatabase

            # Mark the database as modified
            self.db_version_modified = 1
            

    def update_external_database_external_entity_attributes( self, externalDatabase ):
        
    	self.db.insert_db_content( self.db._get_update_sql_query( table = self.biana_database.EXTERNAL_DATABASE_TABLE,
                                                                  update_column_values = (("parsingTime",externalDatabase.get_parsing_time()),),
                                                                  fixed_conditions = (("externalDatabaseID","=",externalDatabase.get_id()),) ),
                                   answer_mode = None )
    																
        

        for current in externalDatabase.get_valid_external_entity_attribute_type():
            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eE_ATTRIBUTE_TABLE,
                                                                      column_values = ( ("externalDatabaseID",externalDatabase.get_id()),
                                                                                        ("attributeType",current) )),
                                       answer_mode = None )


        for current in externalDatabase.get_valid_external_entity_type():
            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eE_TYPES_TABLE,
                                                                      column_values = ( ("externalDatabaseID",externalDatabase.get_id()),
                                                                                        ("eEType",current) )),
                                       answer_mode = None )      

        for current in externalDatabase.get_valid_external_entity_relation_attribute_type():
            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eEr_ATTRIBUTE_TABLE,
                                                                      column_values = ( ("externalDatabaseID",externalDatabase.get_id()),
                                                                                        ("attributeType",current) )),
                                       answer_mode = None )

        for current in externalDatabase.get_valid_external_entity_relation_type():
            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_DATABASE_AVAILABLE_eEr_TYPES_TABLE,
                                                                      column_values = ( ("externalDatabaseID",externalDatabase.get_id()),
                                                                                        ("eErType",current) )),
                                       answer_mode = None )

        # add this database to the default unification protocol...
        # TODO


    def _get_externalDatabaseID_int(self, sourceDBName, sourceDBVersion):
        """
        Returns the decimal value assigned to the sourceDB (externalDatabaseID will be converted to lower_case)

        If "sourceDBVersion" is None, it returns a list with the all the identifiers assigned to the sourceDBName
        """

        try:
            if sourceDBVersion is not None:
                return self._get_valid_source_dbs()[sourceDBName.lower()][sourceDBVersion.lower()].get_id()
            else:
                return [ self._get_valid_source_dbs()[sourceDBName.lower()][x].get_id() for x in self._get_valid_source_dbs()[sourceDBName.lower()] ]
        except:
            sys.stderr.write("Trying to get the identifier for an unexisting external database: %s %s\n" %(sourceDBName, sourceDBVersion))
            raise "Trying to get the identifier for an unexisting external database: %s %s\n" %(sourceDBName, sourceDBVersion)


    ####################################################################################
    #                    EXTERNAL ENTITIES INSERTION METHODS                           #
    ####################################################################################


    # AUTOINCREMENT VALUES

    def _get_new_external_entity_id(self):
        return self.db.get_next_autoincrement(table="externalEntity", attribute="externalEntityID" )

    def _get_last_external_entity_id(self):
        return self.db._get_current_autoincrement(table="externalEntity", attribute="externalEntityID" )

    def _get_last_stable_external_entity_id(self):
        return self.db._get_last_stable_autoincrement(table="externalEntity", attribute="externalEntityID" )

    def _get_new_external_entity_relation_participant_id(self):
        return self.db.get_next_autoincrement(table="externalEntityRelationParticipant", attribute="externalEntityRelationParticipantID" )

    def _get_last_external_entity_relation_participant_id(self):
        return self.db._get_current_autoincrement(table="externalEntityRelationParticipant", attribute="externalEntityRelationParticipantID" )

    def _get_last_stable_external_entity_relation_participant_id(self):
        return self.db._get_last_stable_autoincrement(table="externalEntityRelationParticipant", attribute="externalEntityRelationParticipantID" )

    def _get_new_sequenceProtein_id(self):
        return self.db.get_next_autoincrement(table="sequenceProtein",attribute="proteinSequenceID")

    def _get_last_sequenceProtein_id(self):
        return self.db._get_current_autoincrement(table="sequenceProtein",attribute="proteinSequenceID")

    def _get_last_stable_sequenceProtein_id(self):
        return self.db._get_last_stable_autoincrement(table="sequenceProtein",attribute="proteinSequenceID")

    def _get_new_sequenceNucleotide_id(self):
        return self.db.get_next_autoincrement(table="sequenceNucleotide", attribute="nucleotideSequenceID")

    def _get_last_sequenceNucleotide_id(self):
        return self.db._get_current_autoincrement(table="sequenceNucleotide", attribute="nucleotideSequenceID")

    def _get_last_stable_sequenceNucleotide_id(self):
        return self.db._get_last_stable_autoincrement(table="sequenceNucleotide", attribute="nucleotideSequenceID")

    def _get_new_key_id(self):
        return self.db.get_next_autoincrement(table="externalDatabaseAttributeTransfer", attribute="keyID")

    def _get_last_key_id(self):
        return self.db._get_current_autoincrement(table="externalDatabaseAttributeTransfer", attribute="keyID")

    def _get_last_stable_key_id(self):
        return self.db._get_last_stable_autoincrement(table="externalDatabaseAttributeTransfer", attribute="keyID")
    

    def insert_new_external_entity( self, externalEntity):
        """
        Inserts a new external entity to the biana database

        If the externalEntity has been previously inserted in the database, it will show a message advertising about this. It won't insert the attributes, as an external entity can be inserted only once!

        It is required to use this method once by each externalEntity.

        """

        if externalEntity.get_id() is not None:
            raise ValueError("Trying to insert an inserted external entity...")

        # Mark the database as modified
        self.db_version_modified = 1

        if externalEntity.get_id() is None:

            # Insert the externalEntity and get the new externalEntity ID

            new_external_entity_id = self._get_new_external_entity_id()

            self.db.insert_db_content(self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name(),
                                                                     column_values = ( ("externalEntityID", new_external_entity_id),
                                                                                       ("externalDatabaseID",externalEntity.get_source_database().get_id()),
                                                                                       ("type",externalEntity.get_type())),
                                                                     use_buffer=self.use_buffer ))

            # Insert this externalEntity in the default unification protocol

            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.USER_ENTITY_TABLE,
                                                                      column_values = [("userEntityID",new_external_entity_id),
                                                                                       ("externalEntityID",new_external_entity_id)],
                                                                      use_buffer=self.use_buffer ) )

            externalEntity.set_id(id_value = new_external_entity_id)
            

            externalEntity_attributes = externalEntity.get_attributes_dict()

            # Updates the information that this database is adding to the biana database
            self.valid_source_database_ids[externalEntity.get_source_database().get_id()].add_valid_external_entity_type( eE_type = externalEntity.get_type())

            if externalEntity.get_type()== "ontology":
                self.insert_ontology(externalEntity)

            for current_attribute_identifier in externalEntity_attributes:
                if externalEntity.get_type()== "relation":
                     self.valid_source_database_ids[externalEntity.get_source_database().get_id()].add_valid_external_entity_relation_attribute_type( attribute_identifier = current_attribute_identifier )
                     # Dirty trick to make "psimi_name" to appear in the graphical interface
                     if current_attribute_identifier.lower() == "methodID":
                         self.valid_source_database_ids[externalEntity.get_source_database().get_id()].add_valid_external_entity_relation_attribute_type( attribute_identifier = "psimi_name" ) 
                else:
                     self.valid_source_database_ids[externalEntity.get_source_database().get_id()].add_valid_external_entity_attribute_type( attribute_identifier = current_attribute_identifier )
                    
                [ self._insert_external_entity_attribute( externalEntityID = new_external_entity_id,
                                                          externalEntityAttribute = x ) for x in externalEntity_attributes[current_attribute_identifier] ]


                ### IF THIS ATTRIBUTE IS A TRANSFER ATTRIBUTE, SAVE THE KEY IN THE CORRESPONDING TABLE. THE PARSER HAS TO DEFINE IT PREVIOUS TO ADD ANY THING IN THE DATABASE
                if self._is_key_attribute( attribute_identifier = current_attribute_identifier, external_database_id = externalEntity.get_source_database().get_id() ):
                    for current_key_attribute_object in externalEntity_attributes[current_attribute_identifier]:
                        self.db.insert_db_content( self.db._get_insert_sql_query( table = self._get_key_attribute_table_name( key_id = self.key_attribute_ids[(externalEntity.get_source_database().get_id(), current_attribute_identifier.lower())] ),
                                                                                  column_values = [("externalEntityID",new_external_entity_id),
                                                                                                   ("value",current_key_attribute_object.value)],
                                                                                  use_buffer = self.use_buffer ),
                                               answer_mode = None )

            # Specific for relations
            if externalEntity.get_type()== "relation":

                self.db.insert_db_content(self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_ENTITY_RELATION_TABLE.get_table_name(),
                                                                         column_values = ( ("externalEntityRelationID", new_external_entity_id),
                                                                                           ("type",externalEntity.get_relation_type() ) ),
                                                                         use_buffer=self.use_buffer ),
                                          answer_mode=None)

                self.valid_source_database_ids[externalEntity.get_source_database().get_id()].add_valid_external_entity_relation_type( eEr_type = externalEntity.get_relation_type() )

                for actual_participant_external_entity_id in externalEntity.get_participant_external_entity_ids_list():

                    new_id = self._get_new_external_entity_relation_participant_id()

                    # Save temporarly this information for storing the extended relations hierarchy
                    if self.store_relations_hierarchy:
                        self.temporal_data["relations_hierarchy_parents"].setdefault(actual_participant_external_entity_id, []).append(new_external_entity_id)
                    else:
                        self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE,
                                                                                  column_values = (("externalEntityRelationParticipantID", new_id), 
                                                                                                   (self.biana_database.external_entity_relation_id_col, new_external_entity_id),
                                                                                                   (self.biana_database.externalEntityID_col, actual_participant_external_entity_id )),
                                                                                  use_buffer = True ) )


                    self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE,
                                                                              column_values = (("externalEntityRelationParticipantID", new_id), 
                                                                                               (self.biana_database.external_entity_relation_id_col, new_external_entity_id),
                                                                                               (self.biana_database.externalEntityID_col, actual_participant_external_entity_id) ),
                                                                              use_buffer=self.use_buffer ))
                    

                    participant_attributes = externalEntity.get_participant_attributes( participantExternalEntityID = actual_participant_external_entity_id )

                    [ self._insert_external_entity_relation_participant_attribute( pParticipantID = new_id,
                                                                                   pExternalEntityRelationParticipantAttribute = current_attribute ) for actual_key in participant_attributes  for current_attribute in participant_attributes[actual_key] ]
                                                                   

        else:
            sys.stderr.write("Trying to insert an external Entity which is already in the database! Any new attribute will be inserted!\n")

        return externalEntity.get_id()

        

    def _insert_external_entity_relation_participant_attribute(self, pParticipantID, pExternalEntityRelationParticipantAttribute ):
        """

        """

        tableObject = self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT[pExternalEntityRelationParticipantAttribute.attribute_identifier.lower()]


        column_values = [ ("externalEntityRelationParticipantID",pParticipantID),
                          ("value",pExternalEntityRelationParticipantAttribute.value) ]
        
        for current_field in pExternalEntityRelationParticipantAttribute.additional_fields:
            column_values.append((current_field[0],current_field[1].replace('\\','\\\\').replace('"','\\"')))
        

        self.db.insert_db_content( self.db._get_insert_sql_query( table = tableObject.get_table_name(),
                                                                  column_values = column_values,
                                                                  use_buffer=self.use_buffer ),
                                   answer_mode=None )

        return


    def _insert_external_entity_attribute(self, externalEntityID, externalEntityAttribute):
        """
        Adds a new attribute to the external entity
        externalEntityID must exist previously
        """
        
        tableObject = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[externalEntityAttribute.attribute_identifier.lower()]

        if externalEntityAttribute.attribute_identifier.lower()=="proteinsequence" or externalEntityAttribute.attribute_identifier.lower()=="nucleotidenequence" :
            if externalEntityAttribute.value.get_sequence_MD5() not in self.temporal_data.setdefault("sequence",set()):
            	self._insert_sequence(externalEntityAttribute.value)
            value = externalEntityAttribute.value.get_sequence_MD5()
        else:
            value = externalEntityAttribute.value

        column_values = [ ("externalEntityID",externalEntityID),
                          ("value",value) ]

        for current_field in externalEntityAttribute.additional_fields:
            column_values.append((current_field,externalEntityAttribute.get_field(current_field)))

        if externalEntityAttribute.type is not None:
            column_values.append(("type",externalEntityAttribute.type))

        if externalEntityAttribute.version is not None:
            column_values.append(("version",externalEntityAttribute.version))
        
        self.db.insert_db_content( self.db._get_insert_sql_query( table = tableObject.get_table_name(),
                                                                  column_values = column_values,
                                                                  use_buffer=self.use_buffer ),
                                   answer_mode=None )
        return


    def _insert_sequence(self, sequence):
        """
        Inserts a sequence object into BIANA database and returns its ID.

        """
        
        #if sequence.sequenceID is not None:
        #    column_values = [("sequenceID", sequence.get_sequenceID())]
        #else:
        #    column_values = [("sequenceID",self._get_new_sequenceProtein_id())]

        special_column_values = [ ("sequence", "COMPRESS(\"%s\")" %sequence.get_sequence()) ]

        if sequence.sequence_type == "peptide" :
            if sequence.get_sequenceID() is not None:
                seq_id = sequence.get_sequenceID()
            else:
                seq_id = self._get_new_sequenceProtein_id()
            column_values = [("proteinSequenceID",seq_id)]
            table = self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"]
            column_values.extend([("sequenceMW",sequence.get_proteinMW()),("sequenceIP",sequence.get_proteinIP()) ])

        elif sequence.sequence_type == "rna" or sequence.sequence_type == "dna":
            if sequence.get_sequenceID() is not None:
                seq_id = sequence.get_sequenceID()
            else:
                seq_id = self._get_new_sequenceNucleotide_id()
            column_values = [("nucleotideSequenceID", seq_id)]
            table = self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"]
        
        column_values.extend([("sequenceMD5",sequence.get_sequence_MD5()),
                              ("sequenceLength", sequence.get_length())])

        self.db.insert_db_content(self.db._get_insert_sql_query( table = table,
                                                                 column_values = column_values,
                                                                 special_column_values = special_column_values,
                                                                 use_buffer = self.use_buffer),
                                  answer_mode = None)

        self.temporal_data.setdefault("sequence",set()).add(sequence.get_sequence_MD5())

        return seq_id


    def _insert_sequence_file(self, input_fd, type, format="fasta", verbose=True):
        """
        Inserts a sequence file to the database
        """
        
        num = 0
        
        if type.lower() == "proteinsequence":
            seqObj = BianaObjects.ProteinSequence
        elif type.lower() == "nucleotidesequence":
            seqObj = BianaObjects.NucleotideSequence


        if format.lower()=="fasta":
            temp_seq = []
            for line in input_fd:
                if line[0]==">":
                    if len(temp_seq)>0:
                        self._insert_sequence( sequence = seqObj(sequence="".join(temp_seq),sequenceID=sequenceID ) )
                    sequenceID=line[1:].strip()
                    if( sequenceID == "NULL" or sequenceID == "None" or sequenceID is None ):
                        sequenceID = self._get_new_sequenceProtein_id()
                else:
                    temp_seq.append(line.strip())
                    num += 1
                    if (num%100000)==0 and verbose:
                        sys.stderr.write("%s sequences inserted\n" %num)
 

        elif format.lower() == "seq":
            for line in input_fd:
                splitted = line.strip().split("\t")
                sequenceID = splitted[0]
                if( sequenceID == "NULL" or sequenceID == "None" or sequenceID is None ):
                        sequenceID = self._get_new_sequenceProtein_id()
                self._insert_sequence( sequence = seqObj(sequence=splitted[1],sequenceID=sequenceID ) ) 
                num += 1	
                if (num%100000)==0 and verbose:
                     sys.stderr.write("%s sequences inserted\n" %num)

        else:
            raise ValueError("Format %s not recognized" %(format))


    def _insert_protein_sequence_cd_hit_cluster(self, cd_hit_cluster):
        """
        
        """
        
        # Problem: clusters with only a sequence are not added. It is necessary to add them?
        # It is not necessary, as they are only composed by the representative one...
        # But then... How to be sure that this sequenceID has been used? For the moment, it is not necessary
        for current_match in cd_hit_cluster.get_matches():
            self.db.insert_db_content(self.db._get_insert_sql_query( table = "sequenceProteinCD_HIT",
                                                                     column_values = [("representant_proteinSequenceID", cd_hit_cluster.get_representant()),
                                                                                      ("proteinSequenceID",current_match.sequenceID),
                                                                                      ("representant_start_position",current_match.repr_start_pos),
                                                                                      ("representant_end_position",current_match.repr_end_pos),
                                                                                      ("start_position",current_match.start_pos),
                                                                                      ("end_position",current_match.end_pos),
                                                                                      ("identity",current_match.identity)],
                                                                     use_buffer = self.use_buffer),
                                      answer_mode = None)


    def insert_cd_hit_clusters_to_biana_database(self, cd_hit_clusters_file):
        """
        Inserts information about CD-HIT clusters into database
        """
        ## EXAMPLE OF Clusters output file
        ## >Cluster 0
        ## 0       36805aa, >9662380a987baf844... *
        ## >Cluster 1
        ## 0       35213aa, >0623d1531507bb004... *
        ## 1       90aa, >07b1a5a6aab138163... at 1:90:11834:11923/100%
        ## 2       247aa, >5e00433d0ae984091... at 1:247:12464:12710/100%
        ## 3       153aa, >d845d402ebfa7c203... at 1:153:11430:11582/100%
        ## 4       100aa, >ea8147fd16954563f... at 1:100:11483:11582/100%
        ## 5       183aa, >fab09a87d7f461e75... at 1:183:10973:11155/100%
        ## >Cluster 2
        ## 0       391aa, >1cacc168fca71118c... at 1:391:10887:11277/99%
        ## 1       34942aa, >d861c67fa2f88e4cd... *

        input_file_fd = file(cd_hit_clusters_file, 'r')
        num_clusters=0
        current_cluster = None

        cluster_re = re.compile(">Cluster \d+")
        sequence_re = re.compile("\s*\d+\s+\d+aa,\s+\>(\w+)\.+\s+at\s+(\d+)\:(\d+)\:(\d+)\:(\d+)\/(\d+)\%")
        representant_re = re.compile("\d+\s+\d+aa,\s+\>(\w+)\.+")

        for line in input_file_fd:

            if cluster_re.match(line):
                num_clusters += 1
                if num_clusters % 100000 == 0:
                    sys.stderr.write("%s clusters inserted\n" %(num_clusters))
                if current_cluster is not None:
                    self._insert_protein_sequence_cd_hit_cluster(current_cluster)
                current_cluster = BianaObjects.CDHITCluster()
                continue

            sequence = sequence_re.match(line)
            if sequence:
                current_cluster.add_match( BianaObjects.CDHITMatch( sequenceID= sequence.group(1),
                                                                    start_pos = int(sequence.group(2)),
                                                                    end_pos = int(sequence.group(3)),
                                                                    repr_start_pos = int(sequence.group(4)),
                                                                    repr_end_pos = int(sequence.group(5)),
                                                                    identity = int(sequence.group(6)) ) )
            else:
                # representant
                if line[-2] == "*":
                    sequence = representant_re.search(line)
                    current_cluster.set_representant(sequence.group(1))
                else:
                    sys.stderr.write("Error parsing cd-hit results in line:\n%s" %line)

        self._insert_protein_sequence_cd_hit_cluster(current_cluster)

        # Make sure to empty the buffer...
        self.db._empty_buffer()
        
        return


    def _get_similar_sequences(self, sequenceID, bit_score=None, identity_percent=None, similarity_percent=None, query_coverage_percent=None, match_coverage_percent=None):
        """
        Gets similar sequences from database

        In order to be able to run this method, is mandatory to fill before the table of blast or psi-blast results

        All specified conditions are evaluated jointly, with AND

        TO CHECK!!! The indices are optimized only for identity and coverage!!!

        returns a list of sequenceIDs from similar sequences according to the requeriments.
        """

        sqlStat = database.SQLSelectStatement()
        
        sqlStat.add_element( tables = [self.biana_database.PROTEIN_BLAST_RESULTS_TABLE] )
        sqlStat.add_element( fixed_conditions = [("sequenceID_A","=",sequenceID)] )
        sqlStat.add_element( columns = ["sequenceID_B"] )
        
        if( bit_score is not None ):
            sqlStat.add_element( fixed_conditions = [("bit_score",">=",bit_score)] )
        if( identity_percent is not None ):
            sqlStat.add_element( fixed_conditions = [("identities",">=",identity_percent)])
        if( similarity_percent is not None ):
            sqlStat.add_element( fixed_conditions = [("similarity",">=",similarity_percent)])
        if( query_coverage_percent is not None ):
            sqlStat.add_element( fixed_conditions = [("coverage_A",">=",query_coverage_percent)])
        if( match_coverage_percent is not None ):
            sqlStat.add_element( fixed_conditions = [("coverage_B",">=",match_coverage_percent)] )
        
        query = self.db._get_select_sql_query( tables = sqlStat.tables,
                                               columns = sqlStat.columns,
                                               fixed_conditions = sqlStat.fixed_conditions,
                                               join_conditions = sqlStat.join_conditions )

        print query # JAVI CHECK
        data = self.db.select_db_content( query, answer_mode = "list" )

        #print data
        return data
                                 

    def _insert_blast_results_file(self, file_fd):
        """
        METHOD TO TEST
        """

        columns = ["sequenceID_A","sequenceID_B","evalue","score","bit_score","start_A","end_A","start_B","end_B",
                   "Identities","similarity","gaps","program","filter","coverage_A","coverage_B"]

        for line in file_fd:
            if line.strip().split("\t")[12] != "bl2seq" and line.strip().split("\t")[12] != "blastall":
                sys.stderr.write(line)
            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.PROTEIN_BLAST_RESULTS_TABLE,
                                                                      column_values = zip(columns,line.strip().split("\t")),
                                                                      use_buffer = self.use_buffer ) )
            
        #self.db.insert_db_content("LOAD DATA LOCAL INFILE '%s' INTO TABLE %s (sequenceID_A,sequenceID_B,evalue,score,bit_score,start_A,end_A,start_B,end_B,identities,similarity,gaps,program,filter,coverage_A,coverage_B)" %(file_path, self.biana_database.PROTEIN_BLAST_RESULTS_TABLE),
        #                          answer_mode = None )
                                  

    def _load_sequences(self, sequenceIdList, type="proteinsequence"):
        """
        Gets multiple sequence object from BIANA database

        "sequenceID" can be a single id or a id list

        "type" can be "protein" or "nucleotide"
        """

        results = {}

        if type == "proteinsequence":
            data = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"]],
                                                                             columns = ["proteinSequenceID","UNCOMPRESS(sequence)"],
                                                                             fixed_conditions = [("proteinSequenceID","IN","(%s)" %(", ".join(map(str,sequenceIdList))),None)] ),
                                              answer_mode = "raw" )

            for current_data in data:
                results[current_data[0]] = BianaObjects.ProteinSequence(sequence=current_data[1],sequenceID=current_data[0])

        elif type == "nucleotidesequence":
            data = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"]],
                                                                             columns = ["nucleotideSequenceID","sequenceType","UNCOMPRESS(sequence)"],
                                                                             fixed_conditions = [("nucleotideSequenceID","IN","(%s)" %(", ".join(map(str,sequenceIdList))))] ),
                                              answer_mode = "raw" )
            for current_data in data:
                if current_data[1] == "dna":
                    results[current_data[0]] = BianaObjects.ProteinSequence(sequence=current_data[2],sequenceID=current_data[0])
                elif current_data[1] == "rna":
                    results[current_data[0]] = BianaObjects.ProteinSequence(sequence=current_data[2],sequenceID=current_data[0])

        else:
            raise ValueError("type must be \"proteinsequence\" or \"nucleotidesequence\"")

        return results


    def insert_ontology(self, ontology):
        """
        Inserts to the database the complete ontology object
        
        The external entities of the ontology must be inserted previously
        """

        if not ontology.linked_attribute.lower() in self.biana_database.VALID_EXTERNAL_ENTITY_ATTRIBUTE_TYPES_DICT:
             raise ValueError("Trying to link this ontology to an invalid attribute type")

        key_id = self._add_key_attribute = self._add_key_attribute( external_database_id = ontology.get_source_database().get_id(), key_attribute = ontology.linked_attribute )

        column_values = [ ("externalEntityID", ontology.get_id()),
                          ("linked_attribute", ontology.linked_attribute),
                          ("key_id", key_id),
                          ("name", ontology.name),
                          ("description_attribute", ontology.description_attribute) ]
        
        # In order to simplify and improve the performance of the expansion queries, we need to duplicate data in a table, as in _add_transfer_attribute.
        # It is equivalent to the idea of the key attribute in "transfer attributes"
        # In order to check
       
                
        if ontology.level_attribute is not None:
            column_values.append(("level_attribute", ontology.level_attribute))

        self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.ONTOLOGY_INFO_TABLE,
                                                                  column_values = column_values,
                                                                  use_buffer = self.use_buffer ) )

        for current_eE_id in ontology.get_all_external_entity_ids():
            for current_parent_id in ontology.get_parents_ids(current_eE_id):
                if current_parent_id != current_eE_id:
                    self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.ONTOLOGY_IS_A_TABLE,
                                                                              column_values = [ ("externalEntityID", current_eE_id),
                                                                                                ("is_a",current_parent_id) ],
                                                                              use_buffer = self.use_buffer ) )

            for current_part_parent_id in ontology.get_part_parents_ids(current_eE_id):
                if current_part_parent_id != current_eE_id:
                    self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.ONTOLOGY_IS_PART_OF_TABLE,
                                                                              column_values = [ ("externalEntityID", current_eE_id),
                                                                                                ("is_part_of",current_part_parent_id) ],
                                                                              use_buffer = self.use_buffer ) )

        # Insert extended ontology
        print "Inserting extended ontology"
        for current_eE_id in ontology.get_all_external_entity_ids():
            for current_child in ontology.get_descendants( ontologyElementID = current_eE_id ):
                self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTENDED_ONTOLOGY_HIERARCHY_TABLE,
                                                                          column_values = [ ("parentExternalEntityID",current_eE_id),
                                                                                            ("childExternalEntityID",current_child) ],
                                                                          use_buffer = self.use_buffer ) )
        return

    def get_ontology(self, ontology_name, root_attribute_values = [], load_external_entities=False):
        """
        Loads ontology object
        """

        #print "Getting ontology %s with roots %s" %(ontologyName, root_attribute_values)
        root_attribute_values.sort()
        key = ontology_name.lower()+str(root_attribute_values)

        if key in self.loaded_ontologies:
            return self.loaded_ontologies[key]

        is_a_table = self.biana_database.ONTOLOGY_IS_A_TABLE
        is_part_of_table = self.biana_database.ONTOLOGY_IS_PART_OF_TABLE
        eE_table = self.biana_database.EXTERNAL_ENTITY_TABLE
        eE_field = self.biana_database.externalEntityID_col

        linked_attr = self.get_available_ontology_names(name=ontology_name)
        attr_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[linked_attr.lower()]

        data  = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.ONTOLOGY_INFO_TABLE,
                                                                                    self.biana_database.EXTERNAL_ENTITY_TABLE],
                                                                          columns = ["externalDatabaseID","linked_attribute","key_id","level_attribute","description_attribute"],
                                                                          join_conditions = [("%s.externalEntityID" %self.biana_database.ONTOLOGY_INFO_TABLE,
                                                                                              "=",
                                                                                              "%s.externalEntityID" %self.biana_database.EXTERNAL_ENTITY_TABLE)],
                                                                          fixed_conditions = [("name","=",ontology_name)]),
                                           answer_mode = "raw" )
                                                 
        if len(data) == 0 or len(data)>1:
            raise ValueError("Trying to get an unexsiting ontology or multiple with the same name...")

        ontology = BianaObjects.Ontology( source_database = data[0][0], linkedAttribute = data[0][1], name= ontology_name, descriptionAttribute = data[0][4], levelAttribute=data[0][3] )


        def add_data(sql_results):

            for current_data in sql_results:
                if current_data[1] is not None:
                    is_a = [ int(x) for x in current_data[1].split(",") ]
                else:
                    is_a = []
                if current_data[2] is not None:
                    is_part_of = [ int(x) for x in current_data[2].split(",") ]
                else:
                    is_part_of = []

                ontology.add_element( ontologyElementID = current_data[0],
                                      isA = is_a,
                                      isPartOf = is_part_of,
                                      linkedAttributeValue = current_data[3] )

        # For loading all ontology
        if len(root_attribute_values)==0:
            
            query = "SELECT %s.%s, GROUP_CONCAT(DISTINCT %s), GROUP_CONCAT(DISTINCT %s), value FROM %s LEFT JOIN %s ON %s.%s=%s.%s LEFT JOIN %s ON %s.%s=%s.%s, %s WHERE externalDatabaseID=%s AND %s.%s = %s.%s GROUP BY %s.%s" %(eE_table, eE_field,
                                                                                                                                                                                                                           "is_a", "is_part_of",
                                                                                                                                                                                                                           eE_table,
                                                                                                                                                                                                                           is_a_table,
                                                                                                                                                                                                                           eE_table,eE_field,
                                                                                                                                                                                                                           is_a_table,eE_field,
                                                                                                                                                                                                                           is_part_of_table,
                                                                                                                                                                                                                           eE_table,eE_field,
                                                                                                                                                                                                                           is_part_of_table,eE_field,
                                                                                                                                                                                                                           attr_table,
                                                                                                                                                                                                                           ontology.get_source_database(),
                                                                                                                                                                                                                           eE_table, eE_field,
                                                                                                                                                                                                                           attr_table, eE_field,
                                                                                                                                                                                                                           eE_table, eE_field)

            #print query

            add_data( self.db.select_db_content( query, answer_mode="raw" ) )

        else:
            
            # SELECT ROOT FIRST
            root_list= self.get_list_external_entities_IDs_by_attribute( attribute_identifier = self.get_available_ontology_names(name = ontology_name),
                                                                         field_values = [("value",x) for x in root_attribute_values],
                                                                         source_databases = [ontology.get_source_database() ],
                                                                         expand_ontology_attributes = False )

            root_eE_dict = self.get_external_entities_dict(externalEntityIdsList=root_list,
                                                           attribute_list=[linked_attr])

            add_data( [(x,None,None,root_eE_dict[x].get_attribute(linked_attr).pop().value) for x in root_list ] )

            current_level = map(str, root_list)

            if( len(current_level) == 0 ):
                raise ValueError("Trying to load an ontology with an unexisting root")

            def get_recursive_query(list_root_ids):

                where_sql = "is_a IN (%s)" %(",".join(list_root_ids))
    
                # is part of should be included? TO DECIDE...

                query = "SELECT %s.%s, GROUP_CONCAT(DISTINCT %s), GROUP_CONCAT(DISTINCT %s), value FROM %s LEFT JOIN %s ON %s.%s=%s.%s LEFT JOIN %s ON %s.%s=%s.%s, %s WHERE externalDatabaseID=%s AND %s.%s=%s.%s AND %s GROUP BY %s.%s" %(eE_table, eE_field,
                                                                                                                                                                                                                                    "is_a", "is_part_of",
                                                                                                                                                                                                                                    eE_table,
                                                                                                                                                                                                                                    is_a_table,
                                                                                                                                                                                                                                    is_a_table,eE_field,
                                                                                                                                                                                                                                    eE_table,eE_field,
                                                                                                                                                                                                                                    is_part_of_table,
                                                                                                                                                                                                                                    is_part_of_table,eE_field,
                                                                                                                                                                                                                                    eE_table,eE_field,
                                                                                                                                                                                                                                    attr_table,
                                                                                                                                                                                                                                    ontology.get_source_database(),
                                                                                                                                                                                                                                    eE_table, eE_field,
                                                                                                                                                                                                                                    attr_table, eE_field,
                                                                                                                                                                                                                                    where_sql,
                                                                                                                                                                                                                                    eE_table,eE_field)
                
                return query


            def t(v):
                return str(v[0])

            while( True ):

                if len(current_level)==0:
                    break
                
                data = self.db.select_db_content( get_recursive_query(current_level), answer_mode="raw" )
                add_data( data )
                current_level = map(t,data)


        
        if load_external_entities is True:
            external_entities_dict = {}
             
            #avoid bigg packets...
            all_list = ontology.get_all_external_entity_ids()

            for x in xrange(0,len(all_list)/50):
                current_list = all_list[50*x:50*x+50]
                external_entities_dict.update( self.get_external_entities_dict( externalEntityIdsList = current_list,
                                                                                attribute_list = [ontology.description_attribute],
                                                                                useTransferAttributes = False ) )
            
            ontology._set_external_entities_dict( self.get_external_entities_dict( externalEntityIdsList = ontology.get_all_external_entity_ids(),
                                                                                   attribute_list = [ontology.description_attribute],
                                                                                   useTransferAttributes = False ) )
            
        self.loaded_ontologies[key]  = ontology

        return ontology
                


    ####################################################################################
    #                    EXTERNAL ENTITIES SELECTION METHODS                           #
    ####################################################################################



    def _convertListSourceDBVersionToSourceDBIdList(self, source_databases):
        externalDatabaseID_list = []
        if source_databases is not None:
            for actual_source in source_databases:
                if isinstance(actual_source,int) or isinstance(actual_source,long):
                    externalDatabaseID_list.append(actual_source)
                    continue
                if actual_source[1] is None:
                    externalDatabaseID_list.extend(self._get_externalDatabaseID_int(sourceDBName=actual_source[0],sourceDBVersion=None))
                else:
                    externalDatabaseID_list.append(self._get_externalDatabaseID_int(sourceDBName=actual_source[0],sourceDBVersion=actual_source[1]))
        return externalDatabaseID_list





    def get_external_entities_dict_by_attribute(self, attribute_identifier, field_values, source_databases=None, attribute_restrictions=None, attribute_list=[], relation_attribute_list=[], participant_attribute_list=[] ):
        """
        """ 

        list_eE = self.get_list_external_entities_IDs_by_attribute( attribute_identifier = attribute_identifier, 
                                                                    field_values = field_values,
                                                                    source_databases = source_databases,
                                                                    attribute_restrictions = attribute_restrictions )

        return self.get_external_entities_dict( externalEntityIdsList = list_eE,
                                                attribute_list=attribute_list, 
                                                relation_attribute_list=relation_attribute_list, 
                                                participant_attribute_list=participant_attribute_list)





    def _get_expand_ontology_attribute_table_query(self, attribute_identifier, values_to_expand_list):

        attribute_identifier = attribute_identifier.lower()

        in_str = "(\"%s\")" %"\",\"".join(map(str,values_to_expand_list))

        if not self._is_ontology_linked_attribute(attribute_identifier):
            raise ValueError("Not possible to expand by using %s attribute. Not ontology linked" %attribute_identifier)

        

        attribute_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier]
        
        tables = [ (self.biana_database.EXTENDED_ONTOLOGY_HIERARCHY_TABLE, "o"),
                   (attribute_table, "attr1"),
                   (attribute_table, "attr2") ]
        
        join_conditions = [ ("attr1.value","IN",in_str),
                            ("o.parentExternalEntityID","=","attr1.externalEntityID"),
                            ("o.childExternalEntityID","=","attr2.externalEntityID") ]
        
        columns = [("attr2.value","childs")]
        
        query = self.db._get_select_sql_query( tables = tables,
                                               columns = columns,
                                               join_conditions = join_conditions )
        
        query_list = [query]
        query_list.extend( [ "SELECT %s" %x for x in values_to_expand_list ] )
        
        return "(%s)" %self.db._get_union_queries( queries = query_list )





    def _get_list_external_entities_IDs_by_attribute_SQLstat(self, attribute_identifier, field_values, source_databases=[], attribute_restrictions=None, expand_ontology_attributes=True ):
        """
        Gets a list of external entities that have an attribute with "attribute_value" value
        "field_value" is the field/S and value/S we want to obtain
        If a value is *, it will obtain all having this attribute
        """

        # TODO!!! IF FIELD VALUES HAVE NOT ONLY VALUES... IT WILL CRASH


        sqlStat = database.SQLSelectStatement()

        attribute_identifier = attribute_identifier.lower()

        if attribute_identifier == "proteinsequenceid":

            seq_table = self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"]
            table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT["proteinsequence"]
            sqlStat.add_element( tables = [seq_table, table] )
            sqlStat.add_element( columns = ["%s.externalEntityID" %(table) ] )
            sqlStat.add_element( join_conditions = [("%s.sequenceMD5" %seq_table,
                                                    "=",
                                                    "%s.value" %table)] )
            
        else:

            table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier.lower()]
            sqlStat.add_element( tables = [table] )
            sqlStat.add_element( columns = ["%s.externalEntityID" %(table) ] )

        
        if len(source_databases)>0:
            externalDatabaseID_list = self._convertListSourceDBVersionToSourceDBIdList(source_databases)
            if externalDatabaseID_list != []:
                sqlStat.add_element( tables = [self.biana_database.EXTERNAL_ENTITY_TABLE] )
                sqlStat.add_element( fixed_conditions = [("externalDatabaseID","IN","(%s)" %(",".join([ str(x) for x in externalDatabaseID_list ])),None)] )

                sqlStat.add_element( join_conditions  = [("%s.externalEntityID" %table,
                                                          "=",
                                                          "%s.externalEntityID" %self.biana_database.EXTERNAL_ENTITY_TABLE)] )

        if attribute_restrictions is not None:
            num = 1
            for current_restriction_attribute, restriction_values in attribute_restrictions:
                current_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_restriction_attribute.lower()]
                sqlStat.add_element( tables = [str(current_table)+" AS Q%s" %num] )

                values_list = restriction_values.split(",")   ### WHY????


                if BianaObjects.ExternalEntityAttribute.isFullTextSearchable(attribute_identifier, self.biana_database):
                    sqlStat.add_element( join_conditions = [ ("MATCH (Q%s.value)" %num,
                                                              "AGAINST",
                                                              "('%s' IN BOOLEAN MODE)" %("\",\"".join(values_list))) ] )
                else:

                    if self._is_ontology_linked_attribute( current_restriction_attribute ):
                        sqlStat.add_element( tables = [(self._get_expand_ontology_attribute_table_query( attribute_identifier = current_restriction_attribute,
                                                                                                         values_to_expand_list = values_list ), "childs")] )
                        sqlStat.add_element( join_conditions = [("Q%s.value" %num, "=", "childs.childs")] )

                    else:
                        sqlStat.add_element( join_conditions = [ ("Q%s.value" %num,
                                                                  "IN",
                                                                  "(\"%s\")" %",".join(map(str, values_list))) ] )

                sqlStat.add_element( join_conditions = [ ("%s.externalEntityID" %table,
                                                          "=",
                                                          "Q%s.externalEntityID" %num) ] )
                num += 1
                

        
        # Detect the fields

        values_list = []
        
        for current_field, current_value in field_values:
            if current_field.lower() == "value":
                values_list.append(str(current_value).strip())
            else:
                raise ValueError("TO IMPLEMENT") # by union???

        if BianaObjects.ExternalEntityAttribute.isVersionableIdentifier(attribute_identifier, self.biana_database):
            for current_value in values_list:
                splitted = current_value.split(".")
                if len(splitted)>1:
                    field_values.append(("value",splitted[0]))
                    field_values.append(("version",splitted[1]))
                    raise ValueError("TO IMPLEMENT")
                else:
                    field_values.append(("value",splitted[0]))

        if "*" not in values_list:
            if attribute_identifier == "proteinsequenceid":
                sqlStat.add_element( join_conditions = [("%s.proteinSequenceID" %seq_table,
                                                         "IN",
                                                         "(\"%s\")" % "\",\"".join(values_list) )] )
            elif BianaObjects.ExternalEntityAttribute.isFullTextSearchable(attribute_identifier, self.biana_database):
                sqlStat.add_element( join_conditions = [("MATCH (%s.value)" %table,
                                                         "AGAINST",
                                                         "('%s' IN BOOLEAN MODE)" % " ".join(values_list) )] )
            else:
                if self._is_ontology_linked_attribute( attribute_identifier ):
                    sqlStat.add_element( tables = [(self._get_expand_ontology_attribute_table_query( attribute_identifier = attribute_identifier,
                                                                                                     values_to_expand_list = values_list ), "childs")] )
                    sqlStat.add_element( join_conditions = [("%s.value" %table, "=", "childs.childs")] )
                else:
                    sqlStat.add_element( join_conditions = [("%s.value" %table, "IN", "(\"%s\")" %"\",\"".join(map(str, values_list)))] )

        return sqlStat




    def transform_expanded_attribute_restrictions(self, attribute_restriction_list):
        """

        """
        
        #restrictions = set(attribute_restriction_list)
        restrictions = set()
        
        attribute_values_to_expand = {}
        
        for current_attribute, current_value in attribute_restriction_list:
            if current_attribute.lower() not in self.transferred_attributes:
                restrictions.add((current_attribute, current_value))
                continue
            else:
                attribute_values_to_expand.setdefault(current_attribute.lower(), []).append(current_value)
                
        for current_attribute, values_list in attribute_values_to_expand.iteritems():
            new_attribute_restrictions = self._get_attribute_value_list_for_transfer_attributes( attribute_identifier = current_attribute,
                                                                                                 field_values = [("value",x) for x in values_list] )
            
            map(restrictions.add, new_attribute_restrictions)

        return list(restrictions)




    def _get_attribute_value_list_for_transfer_attributes(self, attribute_identifier, field_values):
        """
        returns a list of (field, value) of the transferred attribute
        """
        
        attribute_identifier = attribute_identifier.lower()

        attribute_values = []

        for current_transfer in self.transferred_attributes[attribute_identifier]:

            # First, get the key attribute values
            tables = [ self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier],            #Attribute transferred
                       (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_transfer[1]],"key1"),    #Attribute used as a key
                       (self.biana_database.EXTERNAL_ENTITY_TABLE, "e1") ]

            columns = ["DISTINCT(key1.value)"]

            fixed_conditions = [("e1.externalDatabaseID",
                                 "=",
                                 current_transfer[0])]

            join_conditions = [("key1.externalEntityID","=","%s.externalEntityID" %(self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier]))]

            if BianaObjects.ExternalEntityAttribute.isFullTextSearchable(attribute_identifier, self.biana_database):
                join_conditions.append(("MATCH (%s.value)" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier],
                                        "AGAINST",
                                        "('%s' IN BOOLEAN MODE)" % " ".join([ str(actual_restriction[1]) for actual_restriction in field_values ] ) ))
            else:
                join_conditions.append(("%s.value" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier],
                                        "IN",
                                        "(\"%s\")" % "\",\"".join([ str(actual_restriction[1]) for actual_restriction in field_values ] ) ))
                                        #"(\"%s\")" % "\",\"".join([ str(actual_restriction) for actual_restriction in field_values ] ) ))

            attribute_values_query = self.db._get_select_sql_query( tables = tables,
                                                                    columns = columns,
                                                                    fixed_conditions = fixed_conditions,
                                                                    join_conditions = join_conditions)

            #print attribute_values_query
            
            attribute_values_list = self.db.select_db_content( attribute_values_query, answer_mode = "list" )

            #attribute_values.append((current_transfer[1],",".join(map(str,attribute_values_list))))
            attribute_values.extend([(current_transfer[1],str(x)) for x in attribute_values_list ])

        return attribute_values




    def _get_list_user_entities_IDs_by_attribute_transfer_query(self, attribute_identifier, field_values, unification_protocol_name, source_databases=[], attribute_restrictions=None, include_type = False, restrict_to_user_entity_ids_list=[] ):


        # MAYBE TO CHANGE... IT WOULD BE BETTER DOING IT WITH THE NEW METHOD

        if unification_protocol_name is None:
            raise ValueError("You are trying to transfer attributes without using a unification protocol...") 

        attribute_identifier = attribute_identifier.lower()
        
        if attribute_identifier not in self.transferred_attributes:
            raise ValueError("Not transferable attribute")


        user_entities_list = []

        for current_transfer in self.transferred_attributes[attribute_identifier]:

            # First, get the key attribute values
            tables = [ self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier],            #Attribute transferred
                       (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_transfer[1]],"key1"),    #Attribute used as a key
                       (self.biana_database.EXTERNAL_ENTITY_TABLE, "e1") ]

            columns = ["DISTINCT(key1.value)"]

            fixed_conditions = [("e1.externalDatabaseID",
                                 "=",
                                 current_transfer[0])]

            join_conditions = [("key1.externalEntityID","=","%s.externalEntityID" %(self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier]))]

            if BianaObjects.ExternalEntityAttribute.isFullTextSearchable(attribute_identifier, self.biana_database):
                join_conditions.append(("MATCH (%s.value)" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier],
                                        "AGAINST",
                                        "('%s' IN BOOLEAN MODE)" % " ".join([ str(actual_restriction[1]) for actual_restriction in field_values ] ) ))
            else:
                join_conditions.append(("%s.value" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier],
                                        "IN",
                                        "(\"%s\")" % "\",\"".join([ str(actual_restriction[1]) for actual_restriction in field_values ] ) ))

            key_values_query = self.db._get_select_sql_query( tables = tables,
                                                              columns = columns,
                                                              fixed_conditions = fixed_conditions,
                                                              join_conditions = join_conditions)
            
            key_values = self.db.select_db_content( key_values_query, answer_mode = "list" )

            key_values = [("value",x) for x in key_values ]

            new_list = self.get_list_user_entities_IDs_by_attribute(unification_protocol_name = unification_protocol_name,
                                                                    attribute_identifier = current_transfer[1],
                                                                    field_values = key_values,
                                                                    attribute_restrictions = attribute_restrictions,
                                                                    restrict_to_user_entity_ids_list = restrict_to_user_entity_ids_list,
                                                                    include_type = include_type)

            user_entities_list.extend(new_list)

        return user_entities_list




    def get_list_external_entities_IDs_by_attribute(self, attribute_identifier, field_values, source_databases=[], attribute_restrictions=None ):

        # Get those external entities directly linked to those attributes
        sqlStat = self._get_list_external_entities_IDs_by_attribute_SQLstat(attribute_identifier=attribute_identifier,
                                                                            field_values=field_values,
                                                                            source_databases=source_databases,
                                                                            attribute_restrictions=attribute_restrictions)

        query = self.db._get_select_sql_query( tables = sqlStat.tables,
                                               columns = sqlStat.columns,
                                               fixed_conditions = sqlStat.fixed_conditions,
                                               join_conditions = sqlStat.join_conditions )

        direct_linked_list = self.db.select_db_content( query, answer_mode = "list" )


        return direct_linked_list




    def _get_attribute_restrictions_query(self, unification_protocol_name, negative_attribute_restrictions ):

        negative_attribute_restrictions = negative_attribute_restrictions

        queries = []

        self._load_available_unification_protocols()
        unif_table = self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)
        

        num_query = 1
        for current_restriction_attribute, current_restriction_values in negative_attribute_restrictions:
            
            #for current_restriction in restrictions_set:
            #print "APPLYING NEGATIVE RESTRICTION ",current_restriction
            current_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_restriction_attribute.lower()]

            tables = [(unif_table,"nu"),
                      (current_table, "nq")]

            columns = ["nu.userEntityID"]

            group_conditions = columns

            join_conditions = [("nq.externalEntityID","=","nu.externalEntityID")]
            
            values_list = [ x for x in str(current_restriction_values).split(",") ]

            #expand_ontology_attributes:
            #values_list.extend( self.expand_ontology_field_values( values_list = values_list, attribute_identifier = current_restriction[0] ) )
            if self._is_ontology_linked_attribute( current_restriction_attribute ):
                tables.append( ( self._get_expand_ontology_attribute_table_query( attribute_identifier = current_restriction_attribute, values_to_expand_list = values_list), "childs" ) )            
                join_conditions.append( ("nq.value","=","childs.childs") )
                
            if BianaObjects.ExternalEntityAttribute.isFullTextSearchable(current_restriction_attribute, self.biana_database):
                join_conditions.append(("MATCH (nq.value)",
                                        "AGAINST",
                                        "('%s' IN BOOLEAN MODE)" %("\",\"".join(map(str,values_list))) ))
            else:
                join_conditions.append(("nq.value",
                                        "IN",
                                        "(\"%s\")" %("\",\"".join(map(str,values_list))) ))

            num_query += 1

            queries.append(self.db._get_select_sql_query( tables = tables,
                                                          columns = columns,
                                                          join_conditions = join_conditions,
                                                          group_conditions = group_conditions ))

        return self.db._get_union_queries( queries )




    def _apply_negative_restrictions_to_query( self, query, unification_protocol_name, negative_attribute_restrictions, column_name_to_restrict="userEntityID" ):
        """
        """

        if len(negative_attribute_restrictions)>0:

            inner_negative_query = self._get_attribute_restrictions_query( unification_protocol_name = unification_protocol_name, 
                                                                           negative_attribute_restrictions = negative_attribute_restrictions )
            tables = [("(%s)" %query,"query_to_restrict_negatively")]

            join_conditions = [("query_to_restrict_negatively.%s" %(column_name_to_restrict),"NOT IN","(%s)" %inner_negative_query)]

            query = self.db._get_select_sql_query( tables = tables,
                                                   columns = ["query_to_restrict_negatively.*"],
                                                   join_conditions = join_conditions )


        return query




    def _apply_restrictions_to_query( self, query, unification_protocol_name, attribute_restrictions, column_name_to_restrict="userEntityID" ):
        """
        Applies the restrictions for a given query, where some user entites must be restricted
        """
        
        unif_table = self._get_user_entity_table_name(unification_protocol_name = unification_protocol_name)
        
        # Restrictions CANNOT be applied at external entity level. They must be applied at USER ENTITY LEVEL.
        # In order to do this, MySQL queries are very slow doing it in a direct way, so, it is faster to do it in nested way

        if attribute_restrictions is not None:

            attribute_restrictions = attribute_restrictions

            num_nested_query = 1

            for current_restriction_attribute, current_restriction_values  in attribute_restrictions:

                current_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_restriction_attribute.lower()]
                query_name = "subrestrictionquery%s" %num_nested_query

                tables = [("(%s)" %query,query_name),
                          (unif_table,"u"),
                          (current_table, "q")]

                columns = [ "%s.*" %query_name ]

                join_conditions = [("%s.%s" %(query_name,column_name_to_restrict),"=","u.userEntityID"),
                                   ("q.externalEntityID","=","u.externalEntityID")]
                
                values_list = [ x for x in str(current_restriction_values).split(",") ]

                if BianaObjects.ExternalEntityAttribute.isFullTextSearchable(current_restriction_attribute, self.biana_database):
                    join_conditions.append(("MATCH (q.value)",
                                            "AGAINST",
                                            "('%s' IN BOOLEAN MODE)" %("\",\"".join(map(str,values_list))) ))

		if BianaObjects.ExternalEntityAttribute.isNumericAttribute(current_restriction_attribute, self.biana_database) or BianaObjects.ExternalEntityAttribute.isSpecialAttribute(current_restriction_attribute, self.biana_database):
                    regex = re.compile("([><=]*)([\d\.]+)")
                    list_of_non_greater_values = []
                    for current_value in values_list:
                        m = regex.match(current_value)
                        if m:
                            join_conditions.append(("q.value",m.group(1),m.group(2)))
                        else:
                            list_of_non_greater_values.append(current_value)
                    
                    if len(list_of_non_greater_values)>0:
                        if self._is_ontology_linked_attribute( current_restriction_attribute ):
                            tables.append( (self._get_expand_ontology_attribute_table_query( attribute_identifier = current_restriction_attribute,
                                                                                             values_to_expand_list = list_of_non_greater_values ), "childs") )
                            join_conditions.append(("q.value","=","childs.childs"))
                        else:
                            join_conditions.append(("q.value",
                                                    "IN",
                                                    "(\"%s\")" %("\",\"".join(map(list_of_non_greater_values)))))
                else:
                    if self._is_ontology_linked_attribute( current_restriction_attribute ):
                        tables.append( (self._get_expand_ontology_attribute_table_query( attribute_identifier = current_restriction_attribute,
                                                                                         values_to_expand_list = values_list ), "childs") )
                        join_conditions.append(("q.value","=","childs.childs"))
                    else:
                        join_conditions.append(("q.value",
                                                "IN",
                                                "(\"%s\")" %("\",\"".join(map(str,values_list)))))
                num_nested_query += 1

                query = self.db._get_select_sql_query( tables = tables,
                                                       columns = columns,
                                                       join_conditions = join_conditions,
                                                       distinct_columns = True )

        return query



    def _apply_relation_restrictions_to_query( self, query, attribute_restrictions_dict, column_name_to_restrict="externalEntityRelationID"):
        """
        """

        num = 1

        for attribute_name, values in attribute_restrictions_dict.iteritems():

            current_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_name.lower()]

            tables = [ ("(%s)" %query, "query_%s" %num),
                       (current_table.get_table_name(), "EERA%s" %num) ]

            columns = [ "query_%s.*" %num ]

            join_conditions = [("query_%s.%s" %(num,column_name_to_restrict),"=","EERA%s.externalEntityID" %num)]
            
            if BianaObjects.ExternalEntityAttribute.isNumericAttribute(attribute_name, self.biana_database) or BianaObjects.ExternalEntityAttribute.isSpecialAttribute(attribute_name, self.biana_database):
                regex = re.compile("([><=]*)([\d\.]+)")
                for current_value in values:
                    m = regex.match(current_value)
                    if m:
                        join_conditions.append(("EERA%s.value" %num,m.group(1),m.group(2)))
                    else:
                        join_conditions.append(("EERA%s.value" %num,"=","\"%s\"" %current_value))
            else:
                join_conditions.append(("EERA%s.value" %num,"IN","(\"%s\")" %("\",\"".join([ str(x) for x in values]))))

            num += 1

            query = self.db._get_select_sql_query( tables = tables,
                                                   columns = columns,
                                                   join_conditions = join_conditions,
                                                   distinct_columns = True )
            
        return query



    def _get_nested_queries_for_getting_equivalent_uE_by_sharing_attributes( self, unification_protocol_name, attribute_list, user_entity_ids_list, restrict_to_user_entity_ids_list=True, ontology_levels={} ):
        """
        Returns a string with a sql query. The sql query returns 2 columns: userEntityID1 and userEntityID2, because they share the combination of attributes in the attribute_list
        """


        # option 1 working, to check if it is efficient...
        unif_table = self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)

        last_query = None

        current_nested_query = 1

        for (current_attribute, parameter_and_value_list) in attribute_list:

            nested_query_tables = [(unif_table,"u1"),
                                   (unif_table,"u2"),
                                   (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute.lower()].get_table_name(),"attr1"),
                                   (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute.lower()].get_table_name(),"attr2")]

            nested_query_columns = [("u1.userEntityID AS userEntityID1"), 
                                    ("u2.userEntityID AS userEntityID2")]
            

            nested_query_join_conditions = [("u1.userEntityID","!=","u2.userEntityID"),
                                            ("u1.externalEntityID", "=", "attr1.externalEntityID"),
                                            ("u2.externalEntityID", "=", "attr2.externalEntityID")]

            nested_query_fixed_conditions = []

            # Non-parameterizable attribute
            if len(parameter_and_value_list) == 0:
                # Use ontologies to determine if they are the same or not
                if self._is_ontology_linked_attribute(current_attribute) and ontology_levels.has_key( self.ontology_linked_attributes[current_attribute.lower()]["level_attribute"] ):
                    # Is ontology linked

                    # SI SOLUCIONO ESTO, SOLUCIONO LAS EXPANSIONES!!!!

                    nested_query_tables.append( (self.biana_database.EXTENDED_ONTOLOGY_HIERARCHY_TABLE, "o") )
                    nested_query_tables.append( (self.biana_database.EXTENDED_ONTOLOGY_HIERARCHY_TABLE, "o2") )
                    nested_query_tables.append( (self._get_key_attribute_table_name( key_id = self.ontology_linked_attributes[current_attribute.lower()]["key_id"] ), "k") )
                    nested_query_tables.append( (self._get_key_attribute_table_name( key_id = self.ontology_linked_attributes[current_attribute.lower()]["key_id"] ), "k2") )
                    nested_query_tables.append( (self._get_key_attribute_table_name( key_id = self.ontology_linked_attributes[current_attribute.lower()]["level_attribute"] ), "level") )

                    nested_query_join_conditions.extend[ ("attr1.value", "=", "k.value"),
                                                         ("k.externalEntityID", "=", "o.childExternalEntityID"),
                                                         ("o.parentExternalEntityID", "=", "level.externalEntityID"),
                                                         ("level.value","=",ontology_levels[self.ontology_linked_attributes[current_attribute.lower()]]["level_attribute"]),
                                                         ("o.parentExternalEntityID", "=", "o2.parentExternalEntityID"),
                                                         ("o2.externalEntityID", "=","k2.externalEntityID"),
                                                         ("k2.value", "=", "attr2.value") ]

                    # TO CHECK TOMORROW... JAVI

                else:
                    nested_query_join_conditions.append(("attr1.value", "=", "attr2.value"))

            # Parameterizable attribute - for instance; sequenceProtein
            else:
                # To make it generic need to change/add user_friendly_names in self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES
                if current_attribute.lower() == "proteinsequence":
                    nested_query_tables.append( ("%s" %(self.biana_database.EXTERNAL_ATTRIBUTE_DESCRIPTIONS_DICT[current_attribute.lower()]["table"].get_table_name()), "d_1") )
                    nested_query_tables.append( ("%s" %(self.biana_database.EXTERNAL_ATTRIBUTE_DESCRIPTIONS_DICT[current_attribute.lower()]["table"].get_table_name()), "d_2") )
                    nested_query_tables.append( ("%s" %(self.biana_database.PROTEIN_BLAST_RESULTS_TABLE.get_table_name()), "b") )
                    nested_query_join_conditions.append(("attr1.value", "=", "d_1.sequenceMD5"))
                    nested_query_join_conditions.append(("attr2.value", "=", "d_2.sequenceMD5"))
                    nested_query_join_conditions.append(("b.sequenceID_A","=","d_1.proteinSequenceID"))
                    nested_query_join_conditions.append(("b.sequenceID_B","=","d_2.proteinSequenceID"))
                    nested_query_fixed_conditions.extend( [ ("b.%s" % parameter, ">=", value) for (parameter, value) in parameter_and_value_list ] )
                else:
                    notSupported = False
                    for (parameter, value) in parameter_and_value_list:
                        if str(value.strip()) != "":
                            notSupported = True
                            break
                    if notSupported:
                        sys.err.write("Warning: Not supported parameterizable attribute (ignored): %s\n" % current_attribute)
                        continue
                    nested_query_join_conditions.extend( [ ("a%s.%s" % (current_attribute, parameter), "=","b%s.%s" % (current_attribute, parameter) ) for (parameter, value) in parameter_and_value_list ] )

            if current_nested_query == 1:
                nested_query_join_conditions.append(("u1.userEntityID","IN","(%s)" %",".join(map(str,user_entity_ids_list))))
                if restrict_to_user_entity_ids_list:
                    nested_query_join_conditions.append(("u2.userEntityID","IN","(%s)" %",".join(map(str,user_entity_ids_list))))
            else:
                nested_query_join_conditions.append(("u1.userEntityID","=","nested_query%s.userEntityID1" %(current_nested_query-1)))
                nested_query_join_conditions.append(("u2.userEntityID","=","nested_query%s.userEntityID2" %(current_nested_query-1)))
                nested_query_tables.append(("(%s)" %last_query, "nested_query%s" %(current_nested_query-1)))


            last_query = self.db._get_select_sql_query( tables = nested_query_tables,
                                                        columns = nested_query_columns,
							fixed_conditions = nested_query_fixed_conditions,
                                                        join_conditions = nested_query_join_conditions,
                                                        distinct_columns = True )

            current_nested_query += 1


	#print last_query
        return self.db._get_select_sql_query( tables = [("(%s)" %last_query, "nested_queries")],
                                              columns = ["userEntityID1","userEntityID2"],
                                              distinct_columns = True )

    
    def get_list_user_entities_IDs_by_attribute(self, unification_protocol_name, attribute_identifier, field_values, attribute_restrictions=None, negative_attribute_restrictions=None, restrict_to_user_entity_ids_list=[], include_type=False ):
        """
        Returns a list of user entities that match with the attributes specified of type attribute_identifier

        If include_type is set to True, it returns a list of tuples 
        """

        sqlStat = self._get_list_external_entities_IDs_by_attribute_SQLstat( attribute_identifier=attribute_identifier,
                                                                             field_values=field_values )

        self._load_available_unification_protocols()
        
        unif_table = self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)

        if attribute_identifier.lower()=="proteinsequenceid":
            attribute_identifier="proteinsequence"

        attr_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier.lower()].get_table_name()
        sqlStat.add_element( tables = [(unif_table,"u")] )
        sqlStat.reset_columns( columns = ["u.userEntityID"] )
        sqlStat.add_element( join_conditions = [("u.externalEntityID",
                                                 "=",
                                                 "%s.externalEntityID" %attr_table )] )
        if include_type is False:
            sqlStat.add_element( group_conditions = ["u.userEntityID"] )
        else:
            sqlStat.add_element( tables = [(self.biana_database.EXTERNAL_ENTITY_TABLE,"e")] )
            sqlStat.reset_columns( columns = ["u.userEntityID","e.type"] )
            sqlStat.add_element( group_conditions = ["u.userEntityID","e.type"] )
            sqlStat.add_element( join_conditions = [("u.externalEntityID","=","e.externalEntityID")] )

        if len(restrict_to_user_entity_ids_list)>0:
            sqlStat.add_element( fixed_conditions = [("u.userEntityID","IN","(%s)" %",".join(map(str,restrict_to_user_entity_ids_list)),None)] )


        # Restrictions CANNOT be applied at external entity level. They must be applied at USER ENTITY LEVEL.
        # In order to do this, MySQL queries are very slow doing it in a direct way, so, it is faster to do it in nesting way

        general_query = self.db._get_select_sql_query( tables = sqlStat.tables,
                                                       columns = sqlStat.columns,
                                                       fixed_conditions = sqlStat.fixed_conditions,
                                                       join_conditions = sqlStat.join_conditions,
                                                       group_conditions = sqlStat.group_conditions )

        restricted_query = self._apply_restrictions_to_query( query = general_query,
                                                              unification_protocol_name = unification_protocol_name,
                                                              attribute_restrictions = attribute_restrictions,
                                                              column_name_to_restrict="userEntityID" )

        negatively_restricted_query = self._apply_negative_restrictions_to_query( query = restricted_query,
                                                                                  unification_protocol_name = unification_protocol_name,
                                                                                  column_name_to_restrict="userEntityID",
                                                                                  negative_attribute_restrictions = negative_attribute_restrictions )


        if len(sqlStat.columns)>1:
            initial_list = list(self.db.select_db_content( negatively_restricted_query, answer_mode = "raw" ))
        else:
            initial_list = self.db.select_db_content( negatively_restricted_query, answer_mode = "list" )

        return initial_list



    ####################################################################################
    #               SPECIFIC METHODS FOR ATTRIBUTE DESCRIPTION DATABASES               #
    ####################################################################################

    def insert_new_attribute_description(self, attribute_identifier, field_values):
        """
        "attribute_identifier" must be an accepted one

        "field_values" must be a dictionary with the keys (accepted ones) and its values
        
        """

        temp_field_values = [ (x.lower(),field_values[x]) for x in field_values ]

        try:
            temp = self.biana_database.EXTERNAL_ATTRIBUTE_DESCRIPTIONS_DICT[attribute_identifier.lower()]
            tableObject = temp["table"]
            if tableObject.has_optimized_fields():
                table = tableObject.get_temp_table_name()
            else:
                table = tableObject.get_table_name()
            fields = temp["fields"]
        except:
            raise "%s is not a valid attribute identifier" %attribute_identifier

        column_values = []

        for actual_field in temp_field_values:
            if actual_field[1] is not None:
                fieldObject = fields[actual_field[0]]

                #check the field name (optimized or not)
                if fieldObject.get_optimized_space() is not None:
                    fieldName = fieldObject.get_optimized_field_name()
                else:
                    fieldName = fieldObject.get_field_name()

                #check if the data is compressed...
                if fieldObject.is_compressed() is not None:
                    #column_values.append((fieldName,"COMPRESS(\"%s\")" %str(actual_field[1]).replace('"',''),1))
                    column_values.append((fieldName,"COMPRESS(\"%s\")" %str(actual_field[1]).replace('\\','\\\\').replace('"','\\"'),1))
                else:
                    #column_values.append((fieldName,str(actual_field[1]).replace('"','')))
                    column_values.append((fieldName,str(actual_field[1]).replace('\\','\\\\').replace('"','\\"')))      

        self.db.insert_db_content(self.db._get_insert_sql_query( table = table,
                                                                 column_values = column_values),
                                                                        #column_values = zip(columns,values) ),
                                  answer_mode=None)

    def output_sequences(self, outmethod, type="proteinsequence", format="fasta", sequenceIDs=None):
        """
        outmethod: destination where data in fasta format would be written
        type: either "proteinsequence" or "nucleotidesequence" - decides which type of sequence is desired
        format: either "fasta" or "seq" corresponding to fasta format or raw sequence format
        sequenceIDs: list of sequence ids whose sequence would be outputted - if None will output all sequences in the database
        """
        # Set block size not to limit mem usage
        database_block_limit=1000000
        database_block_offset=0
        aa_x_line = 79 #80 # amino acid/nucleotide per line
        enter_to_loop = True

        if sequenceIDs is not None:
            enter_to_loop = False
            if type.lower() == "proteinsequence":
                query = "SELECT proteinSequenceID, UNCOMPRESS(sequence) FROM %s WHERE proteinSequenceID IN (%s)" % (self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"].get_table_name(), ",".join(sequenceIDs))
                self.db._check_locked_table(table=self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"].get_table_name())
            elif type.lower() == "nucleotidesequence":
                query = "SELECT nucleotideSequenceID, UNCOMPRESS(sequence) FROM %s WHERE nucleotideSequenceID IN (%s)" % (self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"].get_table_name(), ",".join(sequenceIDs))
                self.db._check_locked_table(table=self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"].get_table_name())
            else:
                raise ValueError("Sequence type not recognized")
            temp_sequences = self.db.select_db_content( query, answer_mode = "raw", remove_duplicates="no" )
        # Gets all the sequences in db
        else: 
            enter_to_loop = True
            if type.lower() == "proteinsequence":
                query = "SELECT proteinSequenceID, UNCOMPRESS(sequence) FROM %s LIMIT %s" %(self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"].get_table_name(), database_block_limit)
                self.db._check_locked_table(table=self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"].get_table_name())
            elif type.lower() == "nucleotidesequence":
                query = "SELECT nucleotideSequenceID, UNCOMPRESS(sequence) FROM %s LIMIT %s" %(self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"], database_block_limit)
                self.db._check_locked_table(table=self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"].get_table_name())
            else:
                raise ValueError("Sequence type not recognized")
            temp_sequences = self.db.select_db_content( query + " OFFSET %s" % database_block_offset, answer_mode = "raw", remove_duplicates="no" )

        while( len(temp_sequences) > 0 ):
            database_block_offset += database_block_limit
            #print "%s sequences done!" %(database_block_offset)
            for actual_sequence in temp_sequences:
                if format.lower() == "fasta":
                        #outmethod(">%s\n%s\n" %(actual_sequence[0],"\n".join([actual_sequence[1][x*aa_x_line:x*aa_x_line+aa_x_line] for x in xrange((len(actual_sequence[1])/aa_x_line)+1) ])))
                        outmethod(">%s\n%s\n" %(actual_sequence[0],"\n".join([ actual_sequence[1][x*aa_x_line:x*aa_x_line+aa_x_line] for x in xrange(int(ceil(len(actual_sequence[1])/float(aa_x_line)))) ])))
                elif format.lower() == "seq":
                        outmethod("%s\t%s\n" %(actual_sequence[0],actual_sequence[1]))
            if not enter_to_loop:
                break
            temp_sequences = self.db.select_db_content( query + " OFFSET %s" % database_block_offset, answer_mode = "raw", remove_duplicates="no" )


    def _empty_sequences_table(self, type):
        """
	
        """

        if type=="proteinsequence":
             table = "sequenceProtein"
        elif type=="nucleotidesequence":
             table = "sequenceNucleotide"
        else:
             raise ValueError("Type %s not recognized to empy sequences table" %(type) )

        self.db.insert_db_content( self.db._get_delete_sql_query(table = table ) )

        return	

    ####################################################################################
    #                              SOME SPECIFIC METHODS                               #
    ####################################################################################

    def get_taxonomy_names_taxID_dict(self, tax_id_name_type=None):
        """
        name is lower cased
        """
	a = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT["taxid"]
	b = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT["taxid_name"]

	self.db._enable_indices( table_list = [a,b] )

        fixed_conditions = []
        if tax_id_name_type is not None:
            fixed_conditions.append(("taxid_name_type","=","scientific name"))

	result = dict(self.db.select_db_content( self.db._get_select_sql_query( tables = [a,b],
										columns = ["LCASE(%s.value)" %b,"%s.value" %a],
										join_conditions = [("%s.externalEntityID" %a, "=", "%s.externalEntityID" %b )],
                                                                                fixed_conditions = fixed_conditions),
                                                 answer_mode = "raw") ) 
	return result
        

    
    def get_taxonomies_from_species_name(self, species_name):

        list_taxID = self.db.select_db_content( self.db.SelectAttributeDescription( attribute_identifier = "ncbi_taxonomy_names",
                                                                                    fields = ["taxID"],
                                                                                    conditions = {"name": species_name} ),
                                                answer_mode = "list", remove_duplicates = "yes" )                                                                                   
        return list_taxID

    
    ####################################################################################
    #                           SEQUENCE RELATED METHODS                               #
    ####################################################################################

    def get_sequence_from_sequenceID( self, sequenceID ):

        # TO CHECK! not only sequence protein...

        return self.db.select_db_content( "SELECT UNCOMPRESS(sequence) FROM sequenceProtein WHERE sequenceProteinID=%s" %sequenceID,
                                          answer_mode = "single" )

    def get_sequence_taxonomies( self, sequenceID ):

        return self.db.select_db_content( "SELECT externalEntityTaxID.value FROM sequenceProtein,externalEntityTaxID,externalEntityProteinSequence WHERE proteinSequenceID=%s AND externalEntityProteinSequence.value = sequenceMD5 AND externalEntityTaxID.externalEntityID = externalEntityProteinSequence.externalEntityID" %sequenceID,
                                          answer_mode = "list", remove_duplicates=True )


    ####################################################################################
    #                           PDB RELATED METHODS                                    #
    ####################################################################################

    def insert_pdb_object( self, PDBObject, source_database, description = None ):
        
        pdb = PDBObject.get_name()

        # Each chain is stored separatedly
        # To each one, a new external Entity Object is created... (I don't know if doing this in this way...

        for actual_chain in PDBObject.get_chain_names():

            externalEntity = BianaObjects.ExternalEntity( source_database= source_database, type="structure")

            externalEntity.add_attribute(BianaObjects.ExternalEntityAttribute( attribute_identifier = "proteinsequence", 
                                                                               value = BianaObjects.ProteinSequence(sequence =  PDBObject.get_sequence( chain_name = actual_chain ))))
                                                                                                                                                             

            residues_list = PDBObject.get_residues( chain_name = actual_chain )
            temp_residue_num_list = [ int(x.get_residue_num()) for x in residues_list ]

            externalEntity.add_attribute(BianaObjects.ExternalEntityAttribute( attribute_identifier = "pdb",
                                                                               value = pdb,
                                                                               type = "unique",
                                                                               additional_fields = { "chain": actual_chain,
                                                                                                     "pdb_range": "%s-%s" %(min(temp_residue_num_list),
                                                                                                                        max(temp_residue_num_list))}))
            if description is not None:
                externalEntity.add_attribute(BianaObjects.ExternalEntityAttribute(attribute_identifier = "description", value = description))

            self.insert_new_external_entity( externalEntity = externalEntity )

            #print "Inserting to database..."

            total_num_atoms = PDBObject.get_chain_num_atoms( chain_name = actual_chain )

            #residue_num_list = unsigned_int_list_to_ascii(2,[ x.get_residue_num() for x in residues_list ]).replace('\\','\\\\').replace('"','\\"')
            residue_num_list = int_ascii.int_list_to_ascii(2,temp_residue_num_list).replace('\\','\\\\').replace('"','\\"')

            residue_type_list = '\x00'.join([ x.get_residue_type() for x in residues_list ]).replace('\\','\\\\').replace('"','\\"')

            #print [ x for x in residues_list ]
            #print [ x.get_num_atoms() for x in residues_list ]
            residue_atom_num_list = int_ascii.unsigned_int_list_to_ascii(1,[ x.get_num_atoms() for x in residues_list ]).replace('\\','\\\\').replace('"','\\"')
            
            atom_type_list = "".join([y.get_type() for x in residues_list for y in x.get_atoms() ]).replace('\\','\\\\').replace('"','\\"')
            
            atom_name_list = '\x00'.join([y.get_name() for x in residues_list for y in x.get_atoms() ]).replace('\\','\\\\').replace('"','\\"')

            atom_coordinates = []
            [ atom_coordinates.extend(y.get_coordinates()) for x in residues_list for y in x.get_atoms() ]

            atom_coordinates_list = int_ascii.float_list_to_ascii(3,3,atom_coordinates).replace('\\','\\\\').replace('"','\\"')

            resolution = PDBObject.resolution

            column_values = [ ("pdb",pdb),
                              ("chain",actual_chain),
                              ("externalEntityID", externalEntity.get_id()),
                              ("total_num_atoms",total_num_atoms) ]

            if resolution is not None:
            	column_values.append(("resolution",str(resolution)))

            self.db.insert_db_content( self.db._get_insert_sql_query( table = "pdb",
                                                                      column_values =  column_values,
                                                                      special_column_values = [("residue_num_list","COMPRESS(\"%s\")" %residue_num_list),
                                                                                               ("residue_type_list","COMPRESS(\"%s\")" %residue_type_list),
                                                                                               ("residue_atom_num_list","COMPRESS(\"%s\")" %residue_atom_num_list),
                                                                                               ("atom_type_list","COMPRESS(\"%s\")" %atom_type_list),
                                                                                               ("atom_name_list","COMPRESS(\"%s\")" %atom_name_list),
                                                                                               ("atom_coordinates","COMPRESS(\"%s\")" %atom_coordinates_list) ],
                                                                      use_buffer=self.use_buffer),
                                       answer_mode = None)


    def add_hssp_info_to_pdb(self,pdb_code, chain, residue_pdb_number, residue_hssp_entropy, residue_hssp_norm_entropy, residue_hssp_variability, conservation_hssp, solvent_exposure_hssp, dssp):
        """
        Adds HSSP information to pdb files

        If the PDB information does not exist, inserts the data 
        """

        print "inserting %s to chain %s" %(residue_pdb_number,chain)
        # It is necessary to lock table pdb, as it is not inserted in the usual and automathic way
        self.db._check_locked_table(table="pdb")

        accessibility = int_ascii.unsigned_int_list_to_ascii(1,solvent_exposure_hssp).replace('\\','\\\\').replace('"','\\"')
        variability = int_ascii.unsigned_int_list_to_ascii(1,residue_hssp_variability).replace('\\','\\\\').replace('"','\\"')
        entropy = int_ascii.unsigned_float_list_to_ascii(2,3,residue_hssp_entropy).replace('\\','\\\\').replace('"','\\"')
        norm_entropy = int_ascii.unsigned_int_list_to_ascii(1,residue_hssp_norm_entropy).replace('\\','\\\\').replace('"','\\"')
        conservation = int_ascii.unsigned_float_list_to_ascii(1,2,conservation_hssp).replace('\\','\\\\').replace('"','\\"')
        residue_numbers = int_ascii.int_list_to_ascii(2,residue_pdb_number).replace('\\','\\\\').replace('"','\\"')

        query = "UPDATE pdb SET hssp_residue_num_correspondences=COMPRESS(\"%s\"),residue_dssp_results=COMPRESS(\"%s\"),residue_hssp_entropy=COMPRESS(\"%s\"),residue_hssp_norm_entropy=COMPRESS(\"%s\"),residue_hssp_variability=COMPRESS(\"%s\"),conservation_hssp=COMPRESS(\"%s\"),solvent_exposure_hssp=COMPRESS(\"%s\") WHERE pdb=\"%s\" AND chain=\"%s\"" %(residue_numbers,"".join(dssp),entropy,norm_entropy,variability,conservation,accessibility,pdb_code,chain)
        #print query
        self.db.insert_db_content( query, answer_mode = None )


    def load_pdb_object(self, pdb_code, fragments=[], merge_fragments = False, request_information=None):
        """
        "pdb_code" is mandatory

        "fragments" is a list of pdb fragments.

        "merge_fragments" is used to merge all fragments or chains as a single chain or to mantain each chain as a different chain. If it is used, atom numbers and residue numbers are changed

        "request_information" indicates which information has to be loaded. It is a list. It can be:
        It can be: "residue_type","atoms_info","hssp_conservation","hssp_entropy","hssp_exposure","hssp_norm_entropy","hssp_variability","dssp"
        """

        information = ["residue_type","atoms_info"]

        if isinstance(request_information,list):
            information.extend(request_information)
        
        
        pdbObject = BianaObjects.PDB( name=pdb_code )


        #specify the pdb code
        fixed_conditions = [("pdb.pdb","=",pdb_code)]

        # specify the chains to obtain
        requested_chains = [ x.get_chain() for x in fragments ]
        
        if len(requested_chains)>0 and None not in requested_chains:
            fixed_conditions.append(("pdb.chain","IN","(\"%s\")" %("\",\"".join(requested_chains)),None))

        possible_requested_data = {"residue_type": "UNCOMPRESS(residue_type_list)",
                                   "atoms_info": ["UNCOMPRESS(atom_type_list)",
                                                  "UNCOMPRESS(atom_name_list)",
                                                  "UNCOMPRESS(atom_coordinates)",
                                                  "UNCOMPRESS(residue_atom_num_list)"],
                                   "hssp_residue_num_correspondences":"UNCOMPRESS(hssp_residue_num_correspondences)",                                   
                                   "hssp_conservation":"UNCOMPRESS(conservation_hssp)",
                                   "hssp_entropy":"UNCOMPRESS(residue_hssp_entropy)",
                                   "hssp_exposure":"UNCOMPRESS(solvent_exposure_hssp)",
                                   "hssp_norm_entropy":"UNCOMPRESS(residue_hssp_norm_entropy)",
                                   "hssp_variability": "UNCOMPRESS(residue_hssp_variability)",
                                   "dssp": "UNCOMPRESS(residue_dssp_results)"}

        columns = ["pdb","chain","resolution","UNCOMPRESS(residue_num_list)"]
        pos = len(columns)
        requested_data_indexes = {}
        append_hssp_eq = None
        try:
            for actual_information in information:
                if actual_information != "atoms_info":
                    if actual_information != "residue_type":
                        append_hssp_eq = 1
                    columns.append(possible_requested_data[actual_information.lower()])
                    requested_data_indexes[actual_information.lower()] = pos
                    pos += 1
                else:
                    columns.extend(possible_requested_data["atoms_info"])
                    requested_data_indexes[actual_information.lower()] = pos
                    pos += 4
        except:
            raise ValueError("Trying to get unavailable pdb information")

        if append_hssp_eq is not None:
            columns.append(possible_requested_data["hssp_residue_num_correspondences"])
            requested_data_indexes["hssp_residue_num_correspondences"] = pos

        query = self.db._get_select_sql_query( tables= ["pdb"],
                                               columns = columns,
                                               fixed_conditions = fixed_conditions )
                           
        
        data = self.db.select_db_content( query, answer_mode = "raw", remove_duplicates = "no" )

        atom_num = 1
        residue_num_value = 1

        for actual_data in data:
            chain = actual_data[1]
            if merge_fragments:
                new_chain = 'A'
            else:
                new_chain = chain
            residue_num_list = int_ascii.ascii_to_int_list(2,actual_data[3])
            residue_type_list = None
            residue_type_value = None
            atom_type_list = None
            atom_name_list = None
            atom_coordinates_list = None
            residue_atom_num_list = None
            hssp_residue_num_correspondence = None
            hssp_conservation = None
            hssp_entropy = None
            hssp_exposure = None
            hssp_norm_entropy = None
            hssp_variability = None
            hssp_conservation_value = None
            hssp_entropy_value = None
            hssp_exposure_value = None
            hssp_norm_entropy_value = None
            hssp_variability_value = None
            dssp = None
            dssp_value = " "
            
            for actual_requested in requested_data_indexes:

                if actual_requested == "residue_type":
                    residue_type_list = actual_data[requested_data_indexes["residue_type"]].split('\x00')
                    #print residue_type_list
                    continue
                if actual_requested == "atoms_info":
                    #print "Getting atoms info"
                    atom_type_list = actual_data[requested_data_indexes["atoms_info"]]
                    atom_name_list = actual_data[requested_data_indexes["atoms_info"]+1].split('\x00')
                    atom_coordinates_list = int_ascii.ascii_to_float_list(3,3,actual_data[requested_data_indexes["atoms_info"]+2])
                    residue_atom_num_list = int_ascii.ascii_to_unsigned_int_list(1,actual_data[requested_data_indexes["atoms_info"]+3])
                    #print atom_type_list
                    #print atom_name_list
                    #print atom_coordinates_list
                    #print residue_atom_num_list
                    continue
                if actual_requested == "hssp_conservation":
                    temp = actual_data[requested_data_indexes["hssp_conservation"]]
                    if temp is not None:
                        hssp_conservation = int_ascii.ascii_to_unsigned_float_list(1,2,temp)
                    else:
                        print "Trying to get the conservation for pdb %s and chain %s that has no data" %(actual_data[0],actual_data[1])
                    continue
                if actual_requested == "hssp_entropy":
                    temp = actual_data[requested_data_indexes["hssp_entropy"]]
                    if temp is not None:
                        hssp_entropy = int_ascii.ascii_to_unsigned_float_list(2,3,temp)
                    else:
                        print "Trying to get the entropy for pdb %s and chain %s that has no data" %(actual_data[0],actual_data[1])
                    continue
                if actual_requested == "hssp_exposure":
                    temp = actual_data[requested_data_indexes["hssp_exposure"]]
                    if temp is not None:
                        hssp_exposure = int_ascii.ascii_to_unsigned_int_list(1,temp)
                    else:
                        print "Trying to get the exposure for pdb %s and chain %s that has no data" %(actual_data[0],actual_data[1])
                    continue
                if actual_requested == "hssp_norm_entropy":
                    temp = actual_data[requested_data_indexes["hssp_norm_entropy"]]
                    if temp is not None:
                        hssp_norm_entropy = int_ascii.ascii_to_unsigned_int_list(1,temp)
                    else:
                        print "Trying to get the normalized entropy for pdb %s and chain %s that has no data" %(actual_data[0],actual_data[1])
                    continue
                if actual_requested == "hssp_variability":
                    temp = actual_data[requested_data_indexes["hssp_variability"]]
                    if temp is not None:
                        hssp_variability = int_ascii.ascii_to_unsigned_int_list(1,temp)
                    else:
                        print "Trying to get the variability for pdb %s and chain %s that has no data" %(actual_data[0],actual_data[1])
                    continue
                if actual_requested == "dssp":
                    temp = actual_data[requested_data_indexes["dssp"]]
                    if temp is not None:
                        dssp = temp
                    else:
                        print "Trying to get the dssp for pdb %s and chain %s that has no data" %(actual_data[0],actual_data[1])
                    continue
                if actual_requested == "hssp_residue_num_correspondences":
                    temp = actual_data[requested_data_indexes["hssp_residue_num_correspondences"]]
                    if temp is not None:
                        hssp_residue_num_correspondence = int_ascii.ascii_to_int_list(2,temp)
                    else:
                        print "Trying to get the hssp equivalences for pdb %s and chain %s that has no data" %(actual_data[0],actual_data[1])
                    continue


            # Get the equivalences between hssp and pdb
            if hssp_residue_num_correspondence is not None:
                hssp_num = 0
                hssp_pdb_equivalences = []
                #print hssp_residue_num_correspondence
                #print residue_num_list
                for x in xrange(len(residue_num_list)):
                    if hssp_num>=len(hssp_residue_num_correspondence):
                        hssp_pdb_equivalences.append(None)
                        continue
                    if residue_num_list[x] == hssp_residue_num_correspondence[hssp_num]:
                        hssp_pdb_equivalences.append(x)
                        hssp_num+=1
                    else:
                        hssp_pdb_equivalences.append(None)

            # Create all the atoms and add them to the pdbObject
            actual_atom_pos = 0

            for x in xrange(len(residue_num_list)):
                residue_requested = None
                #print residue_num_list[x]

                if len(fragments)>0:
                    if True in [ current_fragment.includes( chain=chain, res_num=residue_num_list[x] ) for current_fragment in fragments ]:
                        residue_requested = True
                else:
                    residue_requested = True

                #TODO: Put the equivalences with the hssp numeration...
                if residue_requested:

                    if hssp_residue_num_correspondence is not None and hssp_pdb_equivalences[x] is not None:
                        if hssp_conservation is not None:
                            hssp_conservation_value = hssp_conservation[hssp_pdb_equivalences[x]]
                        if hssp_entropy is not None:
                            hssp_entropy_value = hssp_entropy[hssp_pdb_equivalences[x]]
                        if hssp_exposure is not None:
                            hssp_exposure_value = hssp_exposure[hssp_pdb_equivalences[x]]
                        if hssp_norm_entropy is not None:
                            hssp_norm_entropy_value = hssp_norm_entropy[hssp_pdb_equivalences[x]]
                        if hssp_variability is not None:
                            hssp_variability_value = hssp_variability[hssp_pdb_equivalences[x]]
                        if dssp is not None:
                            dssp_value = dssp[hssp_pdb_equivalences[x]]
                    if residue_type_list is not None:
                        residue_type_value = residue_type_list[x]
                    if merge_fragments:
                        residue_num_value += 1
                    else:
                        residue_num_value = residue_num_list[x]
                                
                    pdbObject.add_residue( chain_name = new_chain,
                                           residue_num = residue_num_value,
                                           residue_type = residue_type_value,
                                           atoms_initial_list = [],
                                           hssp_conservation = hssp_conservation_value,
                                           hssp_entropy = hssp_entropy_value,
                                           hssp_exposure = hssp_exposure_value,
                                           hssp_norm_entropy = hssp_norm_entropy_value,
                                           hssp_variability = hssp_variability_value,
                                           dssp = dssp_value )
                    
                    if requested_data_indexes.has_key("atoms_info"):
                        #print residue_atom_num_list[x]
                        for y in xrange(residue_atom_num_list[x]):
                            current_atom = BianaObjects.PDBAtom( atom_num = atom_num,
                                                                 atom_type = atom_type_list[actual_atom_pos+y],
                                                                 atom_name = atom_name_list[actual_atom_pos+y],
                                                                 x = atom_coordinates_list[(actual_atom_pos+y)*3],
                                                                 y = atom_coordinates_list[(actual_atom_pos+y)*3+1],
                                                                 z = atom_coordinates_list[(actual_atom_pos+y)*3+2])
                            atom_num+=1
                            #actual_atom_pos +=1
                            pdbObject.add_atom( atomObject= current_atom,
                                                chain_name = new_chain,
                                                residue_num = residue_num_value,
                                                residue_type = residue_type_list[x] )

                actual_atom_pos += residue_atom_num_list[x]
                            
        #print pdbObject.get_in_pdb_format()
        #print pdbObject.get_entropy("A")

        if pdbObject.get_num_residues() == 0:
            raise ValueError("Empty PDB created")
        else:
            print pdbObject.get_num_residues()
                
        return pdbObject

                
            #residue_num_list = int_ascii.ascii_to_int_list(2,actual_data[2])
            #print residue_num_list
        


    ####################################################################################
    #                    ALIGNMENTS RELATED METHODS                                    #
    ####################################################################################

    def insert_alignment( self, alignmentObject, externalEntityID, insert_aligned="no" ):
        """
        "insert_aligned" is used to insert exactly the aligned sequence (used in done alignments, as HSSP, as the sequences does not always is exactly the same
                         In this case, the fragments are not inserted
        """
        
        for x in xrange(len(alignmentObject.get_alignment())):

            column_values = [("externalEntityID",externalEntityID),
                             #("sequenceMD5",str(alignmentObject.sequence_ids_list[x]).replace('\\','\\\\').replace('"','\\"')),
                             #("sequenceMD5",str(alignmentObject.sequence_ids_list[x])),
                             ("position",x)]

            if alignmentObject.sequence_ids_list[x] is not None:
                column_values.append(("sequenceMD5",str(alignmentObject.sequence_ids_list[x])))
            if alignmentObject.cross_ids_list[x] is not None:
                column_values.append(("crossID",str(alignmentObject.cross_ids_list[x][1])))
                column_values.append(("crossID_identifier_type",alignmentObject.cross_ids_list[x][0]))

            special_column_values = []
            
            #print "Inserting eE %s with sequenceID %s to position %s" %(externalEntityID,alignmentObject.sequence_ids_list[x],x)

            if insert_aligned == "no":
                column_values.append(("alignmentBlocks",alignmentObject.get_sequence_fragments_in_ascii(x)))
            else:
                special_column_values.append(("alignedSequence","COMPRESS(\"%s\")" %alignmentObject.get_aligned_sequence(x)))
                
            if alignmentObject.sim[x] is not None:
                column_values.append(("identity",alignmentObject.sim[x]))
            if alignmentObject.wsim[x] is not None:
                column_values.append(("weighted_sim",alignmentObject.wsim[x]))
            if alignmentObject.lali[x] is not None:
                column_values.append(("sequence_covered",alignmentObject.lali[x]))


            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.ALIGNMENT_ELEMENT_TABLE,
                                                                      column_values = column_values,
                                                                      special_column_values = special_column_values,
                                                                      use_buffer=self.use_buffer),
                                       answer_mode = None )


    def get_alignment( self, externalEntityID, get_species=False, range=None ):
        """
        Returns an alignment object identified by the externalEntityID used as parameter
        """

        data = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.ALIGNMENT_ELEMENT_TABLE],
                                                                         columns = ["sequenceMD5","position","alignmentBlocks","UNCOMPRESS(alignedSequence)","crossID","crossID_identifier_type"],
                                                                         fixed_conditions = [(self.biana_database.externalEntityID_col,"=",externalEntityID)]),
                                          answer_mode = "raw" )

        aln = BianaObjects.SequenceAlignment()

        aln.set_number_of_sequences( number_of_sequences = len(data) )

        seq_id_tax_id_dict = {}

        temp_list_seq_ids = [ actual_aligned_seq[0] for actual_aligned_seq in data ]

        def start_list(x):
            seq_id_tax_id_dict[x]=[]

        map(start_list,temp_list_seq_ids)

        #[ seq_id_tax_id_dict[x] = [] for x in temp_list_seq_ids ]


        # FOR DOING IT USING SEQUENCE ID... NOT USED NOW.
        #if get_species == "yes":            
        #    sequenceID_taxID_correspondences = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_ENTITY_ATTRIBUTES_DICT["taxid"]["table"],
        #                                                                                                           self.biana_database.EXTERNAL_ENTITY_ATTRIBUTES_DICT["sequence"]["table"]],
        #                                                                                                 columns = ["sequenceMD5","taxID"],
        #                                                                                                 fixed_conditions = [("sequenceMD5","IN","(%s)" %(",".join([str(x) for x in temp_list_seq_ids])),None)],
        #                                                                                                 join_conditions = [("%s.externalEntityID" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTES_DICT["taxid"]["table"],
        #                                                                                                                     "=",
        #                                                                                                                     "%s.externalEntityID" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTES_DICT["sequence"]["table"])]),
        #                                                                  
        #                                                                  answer_mode = "raw" )
            
         #   for actual_seq_taxID in sequenceID_taxID_correspondences:
         #       seq_id_tax_id_dict[actual_seq_taxID[0]].append(actual_seq_taxID[1])

        tax_ids = []

        # COMMENTED AS HSSP GIVES SPECIE... BUT IN THE FUTURE IT MUST BE UNCOMMENTED
        #if get_species == True:
        #    taxid_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT["taxid"]
        #    for actual_aligned_seq in data:
        #        attr_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[actual_aligned_seq[5].lower()]
        #        query = self.db._get_select_sql_query( tables = [taxid_table,attr_table],
        #                                               columns = ["DISTINCT(%s.value)" %taxid_table],
        #                                               fixed_conditions = [("%s.value" %attr_table,"=",actual_aligned_seq[4])],
        #                                               join_conditions = [("%s.externalEntityID" %attr_table,
        #                                                                   "=",
        #                                                                   "%s.externalEntityID" %taxid_table)] )
        #        print query
        #        tax_ids.append( self.db.select_db_content( sql_query = query, answer_mode = "list" ) )
        #else:
        #    tax_ids = [[] for x in xrange(len(data))]


        for x_actual_aligned_seq in xrange(len(data)):
            actual_aligned_seq = data[x_actual_aligned_seq]
            try:
                tax_id_list = [ actual_aligned_seq[4].split("_")[1] ]
            except:
                tax_id_list = []
            aln.add_sequence_to_position( sequence_id = actual_aligned_seq[0],
                                          #taxID_list = seq_id_tax_id_dict[actual_aligned_seq[0]],
                                          #taxID_list = tax_ids[x_actual_aligned_seq],  # TO UNCOMMENT IF UNCOMMENTED ABOVE
                                          taxID_list = tax_id_list,
                                          sequence_position = actual_aligned_seq[1],
                                          sequence_fragments_ascii = actual_aligned_seq[2],
                                          aligned_sequence = actual_aligned_seq[3],
                                          crossID = (actual_aligned_seq[5],actual_aligned_seq[4]) )

        #aln.print_alignment()

        return aln


    def load_hssp_multiple_chains(self, pdb_id, fragments=[], get_species=False):
        """
        
        """
        final_aln = None

        for current_fragment in fragments:
            new_aln = self.load_hssp_alignment(pdb_id=pdb_id, chain=current_fragment.get_chain(), fragments=[current_fragment], get_species=get_species)
            if final_aln is None:
                final_aln = new_aln
            else:
                final_aln = final_aln.concatenate_by_crossID(new_aln)

        return final_aln

        #by_chains = {}
        #for current_fragment in fragments:
        #    by_chains.setdefault(current_fragment.get_chain(),[]).append(current_fragment)

        #print "Requested %s chains!" %(len(by_chains))
        
        #alns = {}

        #for current_chain in by_chains:
        #    self.load_hssp_alignment(pdb_id=pdb_id, chain=current_chain, fragments=by_chains[current_chain], get_species=get_species)
            
        #if len(by_chains)>2:
        #    raise ValueError("for the moment, only prepared to run with a maximum of two chains")

        
        
        

    def load_hssp_alignment(self, pdb_id, chain, fragments=[], get_species=False):
        """
        Returns an alignment for an structure given in the hssp
        
        For the moment, all the fragments must belong to the same chain
        """

	pdb_attr_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT["pdb"]

        fixed_conditions = [("%s.value" %pdb_attr_table,"=",pdb_id)]
        join_conditions = [("%s.externalEntityID" %pdb_attr_table,"=","%s.externalEntityID" %self.biana_database.ALIGNMENT_ELEMENT_TABLE ),
                           ("%s.pdb" %self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["pdb"],"=","%s.value" %pdb_attr_table)]
    
        if chain is not None:
            fixed_conditions.append(("%s.chain" %pdb_attr_table,"=",chain))
            join_conditions.append(("%s.chain" %pdb_attr_table,"=","%s.chain" %pdb_attr_table))

        query = self.db._get_select_sql_query( tables = [self.biana_database.ALIGNMENT_ELEMENT_TABLE,
                                                         self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["pdb"],
                                                         pdb_attr_table],
                                               columns = ["%s.externalEntityID" %pdb_attr_table,
                                                          "UNCOMPRESS(hssp_residue_num_correspondences)",
                                                          "UNCOMPRESS(residue_num_list)"],
                                               fixed_conditions = fixed_conditions,
                                               join_conditions = join_conditions )

        data = self.db.select_db_content( query, answer_mode = "raw" )
        
        if len(data)==0:
            print "HSSP %s NOT FOUND\n" %pdb_id
            # Search for an homologous chain??
            return None
        else:
            alignmentExternalEntityID = data[0][0]
            hssp_residue_num_correspondence = int_ascii.ascii_to_int_list(2,data[0][1])
            residue_num_list = int_ascii.ascii_to_int_list(2,data[0][2])

            #print residue_num_list
            
            align_fragments_to_select = []
            current_start = 0
            current_end = None

            if len(fragments)>0:
                for actual_fragment in fragments:
                    for x in xrange(len(residue_num_list)):
                        if (actual_fragment.get_start_residue() is None and actual_fragment.get_end_residue() is None) or (actual_fragment.get_start_residue()<=residue_num_list[x] and (actual_fragment.get_end_residue()>=residue_num_list[x] or actual_fragment.get_end_residue() is None)):
                            #print "I take position %s (which corresponds to %s)" %(x,residue_num_list[x])
                            current_end = x
                        else:
                            if current_end is not None:
                                #print "Adding (%s,%s)" %(current_start, current_end)
                                align_fragments_to_select.append((current_start,current_end))
                            current_end = None
                            current_start = x+1

                if current_end is not None:
                    align_fragments_to_select.append((current_start,current_end)) 

            else:
                align_fragments_to_select.append((current_start,len(residue_num_list)-1))

            alignment = self.get_alignment(externalEntityID=alignmentExternalEntityID,get_species=get_species)

            return alignment.get_subalignment( fragments = align_fragments_to_select )

        
    ####################################################################################
    #                    USER ENTITY RELATED METHODS                                   #
    ####################################################################################


    def _load_available_unification_protocols(self):
        """
        Gets the information about available unification protocols
        """

        if( self.available_unification_protocols is None ):
            
            self.available_unification_protocols = {}

            data = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.USER_ENTITY_PROTOCOL_TABLE.get_table_name()],
                                                                             columns = ["description","unificationProtocolID","databaseVersion"]),
                                              answer_mode = "raw" )
            
            for actual_protocol in data:

                uP = BianaObjects.UnificationProtocol(description = actual_protocol[0], BianaDatabaseVersion = actual_protocol[1], id = str(actual_protocol[1]))
                
                # Search all the atoms for this unification protocol
                atoms = self.db.select_db_content( self.db._get_select_sql_query( tables = [(self.biana_database.USER_ENTITY_PROTOCOL_ATOMS_TABLE.get_table_name(),"A"),
                                                                                            (self.biana_database.USER_ENTITY_PROTOCOL_ATOM_ATTRIBUTES_TABLE.get_table_name(),"B")],
                                                                                  columns = ["externalDatabaseID_A","externalDatabaseID_B","GROUP_CONCAT(B.cross_referenced_code)"],
                                                                                  join_conditions = [("A.unificationAtomID","=","B.unificationAtomID")],
										  fixed_conditions = [("A.unificationProtocolID","=",actual_protocol[1])],
                                                                                  group_conditions = ["A.unificationAtomID"] ),
                                                   answer_mode = "raw" )
                for current_atom in atoms:
                    uP.add_unification_atom_elements( BianaObjects.UnificationAtomElement( externalDatabaseID_A = current_atom[0],
                                                                                           externalDatabaseID_B = current_atom[1],
                                                                                           externalAttribute = current_atom[2].split(",") ) )

                self.available_unification_protocols[actual_protocol[0].lower()] = uP


            self.available_unification_protocols["no unification"] = None

        return

    
    def _get_user_entity_table_name(self, unification_protocol_name):
        
        self._load_available_unification_protocols()
        
        unification_protocol_name = str(unification_protocol_name).lower()

        if self.available_unification_protocols[unification_protocol_name] is None:
            return self.biana_database.USER_ENTITY_TABLE.get_table_name()
        
        return self.biana_database.USER_ENTITY_TABLE.get_table_name()+self.available_unification_protocols[unification_protocol_name].get_id()


    def get_available_unification_protocols_list(self):
        
        if( self.dbname is None ):
            return []

        self._load_available_unification_protocols()
        
        return list(self.available_unification_protocols.keys())
    

    def _create_new_unification_protocol_tables(self, protocol):
        """
        """
        
        # Inserting the basic information of the protocol to the database
        protocol_id = self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.USER_ENTITY_PROTOCOL_TABLE,
                                                                                column_values = [("description",protocol.get_description()),
                                                                                                 ("databaseVersion","TEST_DB_VERSION")],
                                                                                use_buffer = False ),
                                                 answer_mode = "last_id" )
        protocol.id = protocol_id

        # Save the unification protocol...
        for actual_unification_atom_element in protocol.get_unification_atom_elements():
            atom_id = self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.USER_ENTITY_PROTOCOL_ATOMS_TABLE,
                                                                                column_values = [ ("unificationProtocolID",protocol_id),
                                                                                                  ("externalDatabaseID_A",actual_unification_atom_element.get_externalDatabaseID_A() ),
                                                                                                  ("externalDatabaseID_B",actual_unification_atom_element.get_externalDatabaseID_B() ) ],
                                                                                use_buffer=False ),
                                                 answer_mode = "last_id" )

            crossed_attributes = actual_unification_atom_element.get_external_attribute_list()

            for actual_crossed_attribute in crossed_attributes:
                self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.USER_ENTITY_PROTOCOL_ATOM_ATTRIBUTES_TABLE,
                                                                          column_values = [ ("unificationAtomID",atom_id),
                                                                                            ("cross_referenced_code",actual_crossed_attribute)],
                                                                          use_buffer = False) )
                               

    def create_new_user_entities(self, protocol):
        """

        "protocol" is the UnificationProtocol object that has to be followed
        """

        # WE COULD CHECK IF A PREVIOUS UNIFICATION WITH THE SAME PARAMETERS HAS BEEN PREVIOUSLY DONE...


        if self.isOptimizedForRunning is False:
            self.optimize_database_for( mode="running", optimize=True )

        self._create_new_unification_protocol_tables(protocol)

        protocol_temp_table = database.TableDB( table_name = "temp_equivalences_protocol_%s" %protocol.id,
                                                table_fields = [ database.FieldDB(field_name = "externalEntity1",
                                                                                  data_type = "integer(4) unsigned"),
                                                                 database.FieldDB(field_name = "externalEntity2",
                                                                                  data_type = "integer(4) unsigned")] )

        self.db.insert_db_content( protocol_temp_table.create_mysql_query(), answer_mode = None )
        
        print "Unifying..."
        
        #import C_functions
        import tempfile

        temp_eq_file = tempfile.NamedTemporaryFile(bufsize=0)

        #Done like this because if not it does not work on WINDOWS os...
        temp_eq_file_name = temp_eq_file.name
        temp_eq_file.close()
        temp_eq_file = open(temp_eq_file_name,'w')

        # Get the list of queries to obtain the equivalences
        for actual_unification_atom_element in protocol.get_unification_atom_elements():
            query = self._get_equivalent_external_entities(actual_unification_atom_element)
	    if query is not None:
		self.db.insert_db_content( "INSERT INTO %s (%s)" %(protocol_temp_table.get_table_name(),query),
					   answer_mode = None )
            
            #data = self.db.select_db_content( query, answer_mode = "raw" )
            #[ temp_eq_file.write("%s\t%s\n" %(equivalence[0],equivalence[1])) for equivalence in data ]



        # Get the information from the temporal table

        l=1000000
        o=0


        temp_equivalences = self.db.select_db_content( "SELECT * FROM %s LIMIT %s OFFSET %s" %(protocol_temp_table.get_table_name(),l,o),
                                                       answer_mode = "raw", remove_duplicates="no" )

        while( len(temp_equivalences) > 0 ):
            o += l
            [ temp_eq_file.write("%s\t%s\n" %(equivalence[0],equivalence[1])) for equivalence in temp_equivalences ]
            temp_equivalences = self.db.select_db_content( "SELECT * FROM %s LIMIT %s OFFSET %s" %(protocol_temp_table.get_table_name(),l,o),
                                                           answer_mode = "raw", remove_duplicates="no" )


        # Delete the content of the temporal table
        self.db.insert_db_content( "DELETE FROM %s" %protocol_temp_table.get_table_name() )

        temp_all_file = tempfile.NamedTemporaryFile(bufsize=0)
        temp_all_file_name = temp_all_file.name
        temp_all_file.close()
        temp_all_file = open(temp_all_file_name, 'w')

        # Get the list of queries to obtain all the external entities of desired databases
        # FIRST, IT IS NEEDED TO OBTAIN A LIST OF ALL DATABASES THAT ARE GOING TO BE USED
        databases = protocol.get_database_ids()
        for actual_database in databases:
	    if not self.get_external_database(database_id = actual_database).get_promiscuity():
		query =  self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_ENTITY_TABLE],
							columns = [self.biana_database.externalEntityID_col],
							fixed_conditions = [("externalDatabaseID","=",actual_database),
									    ("type","!=","relation")] )   #JAVI RECENTLY ADDED. MAY DECREASE UNIFYING EFICIENCY...
            
		#print query

		self.db.insert_db_content( "INSERT INTO %s (externalEntity1) (%s)" %(protocol_temp_table.get_table_name(),query),
					   answer_mode = None )

            #data = self.db.select_db_content( query, answer_mode = "raw" )
            #[ temp_all_file.write("%s\n" %x) for x in data ]

            
        # Get the information from the temporal table
        l=1000000
        o=0
        temp_equivalences = self.db.select_db_content( "SELECT externalEntity1 FROM %s LIMIT %s OFFSET %s" %(protocol_temp_table.get_table_name(),l,o),
                                                       answer_mode = "raw", remove_duplicates="no" )
        while( len(temp_equivalences) > 0 ):
            o += l
            [ temp_all_file.write("%s\n" %equivalence) for equivalence in temp_equivalences ]
            temp_equivalences = self.db.select_db_content( "SELECT externalEntity1 FROM %s LIMIT %s OFFSET %s" %(protocol_temp_table.get_table_name(),l,o),
                                                           answer_mode = "raw", remove_duplicates="no" )


        # Delete the temporal table
        self.db.insert_db_content( protocol_temp_table.get_drop_query() )

        unification_temp_file = tempfile.NamedTemporaryFile(bufsize=0)
        unification_temp_file_name = unification_temp_file.name

        unification_temp_file.close()
        temp_eq_file.close()
        temp_all_file.close()

        from biana import __path__ as biana_path
        import os

        command = "\""+biana_path[0]+os.sep+"BianaDB"+os.sep+"%s\" "+temp_eq_file_name+" "+temp_all_file_name+" "+unification_temp_file_name
    
        if sys.platform.lower().startswith("win"):
            command = command %("win_unify.exe")
        else:
            command = command %("unify")
        
        p = os.system(command)

        if p==1:
            raise ValueError("ERROR in unification. Unify compiled program error")

        temp_eq_file.close()
        temp_all_file.close()


        print "Creating the protocol table"
        protocol_table = copy.deepcopy(self.biana_database.USER_ENTITY_TABLE)  # It is necessary to do a copy because we are going to change its name...
        protocol_table.set_table_name(new_name = "%s%s" %(protocol_table.get_table_name(),protocol.id))

        self.db.insert_db_content( protocol_table.create_mysql_query(ignore_primary_key=True), answer_mode = None )
        
        self.db._disable_indices( table_list = [protocol_table] )



        if False:
            self.db._lock_tables( table_list = [protocol_table.get_table_name()] )
            self.db.insert_db_content( "LOAD DATA LOCAL INFILE '%s' INTO TABLE %s" %(unification_temp_file.name.replace("\\","\\\\"),protocol_table.get_table_name()) )
        else:
            # Done like this because not all mysql versions accept LOAD DATA LOCAL INFILE.
            
            # Read and insert the unification
            unification_file = open(unification_temp_file.name,'r')
            for line in unification_file:
                fields = line.strip().split("\t")
                self.db.insert_db_content( self.db._get_insert_sql_query( table = protocol_table,
                                                                          column_values = [("userEntityID", fields[0]),
                                                                                           ("externalEntityID", fields[1])],
                                                                          use_buffer = "yes" ) )
            unification_file.close()

        import os
        os.unlink(unification_temp_file.name)
        
        # Make sure to empty the buffer...
        self.db._empty_buffer()

        self.db._enable_indices( table_list = [protocol_table] )

        self.db._unlock_tables()

        # Here, we should incorporate into the unification the promiscuous external entities (external entities than can belong to different user entities, as scop or pfam domains, etc...
        self._unify_promiscuous_external_entities( protocol = protocol )

        return


    def _unify_promiscuous_external_entities( self, protocol ):
        """
        
        """

        print "Adding promiscuous external entities"

        #self.db.set_lock_tables( False )

        unification_table = self._get_user_entity_table_name(unification_protocol_name=protocol.description)

        # Get all promiscuous external entities
        for actual_unification_atom_element in protocol.get_unification_atom_elements():
	    db_id_A = actual_unification_atom_element.get_externalDatabaseID_A()
	    db_id_B = actual_unification_atom_element.get_externalDatabaseID_B()

	    # If only one of the external database is promiscuous
	    #xor = lambda x,y: x!=y and (x or y) # using built-in xor: x^y
	    if not (self.get_external_database(database_id = db_id_A).get_promiscuity() ^ self.get_external_database(database_id = db_id_B).get_promiscuity()):
		continue

	    if self.get_external_database(database_id = db_id_A).get_promiscuity():
		promiscuous_db_id, non_promiscuous_db_id = db_id_A, db_id_B
	    else:
		promiscuous_db_id, non_promiscuous_db_id = db_id_B, db_id_A

	    attribute_list = actual_unification_atom_element.get_external_attribute_list()

            # Get the equivalences of the user entities for this promiscuous external entities
            tables = [ (unification_table, "u"),
                       (self.biana_database.EXTERNAL_ENTITY_TABLE, "eE1"),
                       (self.biana_database.EXTERNAL_ENTITY_TABLE, "eE2") ]

            columns = ["userEntityID", "eE2.externalEntityID"]

	    join_conditions = [ ("eE1.externalEntityID", "=", "u.externalEntityID") ]
	    fixed_conditions = [ ("eE2.type","!=","relation"), ("eE1.externalDatabaseID","=",non_promiscuous_db_id), ("eE2.externalDatabaseID","=",promiscuous_db_id) ]

	    for attribute_index in xrange(len(attribute_list)):
		unifying_attribute = attribute_list[attribute_index]
		tables.append( (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[unifying_attribute], "a1_%s" % attribute_index) )
		tables.append( (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[unifying_attribute], "a2_%s" % attribute_index) )
		join_conditions.append( ("a1_%s.value" % attribute_index,"=","a2_%s.value" % attribute_index) )
		join_conditions.append( ("eE1.externalEntityID","=","a1_%s.externalEntityID" % attribute_index) )
		join_conditions.append( ("eE2.externalEntityID","=","a2_%s.externalEntityID" % attribute_index) )
		# do not unify based on previous attributes (emre), may decrease unifying efficiency as well
		fixed_conditions.append( ("a1_%s.type" % attribute_index,"!=","previous") ) 
		fixed_conditions.append( ("a2_%s.type" % attribute_index,"!=","previous") )
		
		if unifying_attribute.lower()=="pdb":
		    join_conditions.append( ("a1_%s.chain"%attribute_index,"=","a2_%s.chain"%attribute_index) )

            get_equivalences_query = self.db._get_select_sql_query( tables = tables, 
                                                                    columns = columns,
                                                                    join_conditions = join_conditions,
                                                                    fixed_conditions = fixed_conditions )
            #print get_equivalences_query

            # Insert them in the user entities table

            # Trying to insert directly, without being in a temporary table.
            # Probably this is slow... If it is:
            #     Insert all information in a temporary table
            #     Disable indices from unification table
            #     Put all information of temporary table into unification table
            #     Enable indices from unification table
            #     Delete temporary table
            
            #self.db.insert_db_content( sql_query = "INSERT INTO %s (userEntityID, externalEntityID) (%s)" %(unification_table,
            #                                                                                                get_equivalences_query ),
            #                           answer_mode = None )

            self.db.insert_db_content( sql_query = self.db._get_nested_insert_sql_query( table = unification_table,
                                                                                         columns = ["userEntityID", "externalEntityID"],
                                                                                         subquery = get_equivalences_query ),
                                       answer_mode = None)
            # Finished


    def drop_unification_protocol( self, unification_protocol_name ):
        """
        Drops from the database all the information of a protocol of unification

        "protocol_description" corresponds to the description of the protocol to drop
        """

        self._load_available_unification_protocols()

        unification_protocol_name = unification_protocol_name.lower()

        if not self.available_unification_protocols.has_key( unification_protocol_name ):
            raise ValueError("ERROR. Trying to drop an unexisting unification protocol")

        if unification_protocol_name == "no unification":
            raise ValueError("Default unification protocol cannot be deleted")

        unificationProtocolObj = self.available_unification_protocols[unification_protocol_name]

        # Delete information about unification protocol atoms
        atom_ids_query = self.db._get_select_sql_query( tables = [self.biana_database.USER_ENTITY_PROTOCOL_ATOMS_TABLE],
                                                        columns = ["unificationAtomID"],
                                                        fixed_conditions = [("unificationProtocolID","=", unificationProtocolObj.get_id())] )

        self.db.insert_db_content( self.db._get_delete_sql_query(table = self.biana_database.USER_ENTITY_PROTOCOL_ATOM_ATTRIBUTES_TABLE,
                                                                 fixed_conditions = [("unificationAtomID","IN","(%s)" %atom_ids_query, None)]) )
        
        self.db.insert_db_content( self.db._get_delete_sql_query(table = self.biana_database.USER_ENTITY_PROTOCOL_ATOMS_TABLE,
                                                                 fixed_conditions = [("unificationProtocolID","=", unificationProtocolObj.get_id())]) )

        # Delete information in unification protocols table
        self.db.insert_db_content( self.db._get_delete_sql_query(table = self.biana_database.USER_ENTITY_PROTOCOL_TABLE,
                                                                 fixed_conditions = [("unificationProtocolID","=", unificationProtocolObj.get_id())]) )

        # Drop unification table
        self.db.insert_db_content( self.db._get_drop_sql_query( [self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)] ) )

        return

    def get_unification_protocol_atoms( self, unification_protocol_name ):
        """
        Fetchs from the database all the atom information of a protocol of unification

        "unification_protocol_name" corresponds to the description of the protocol
        """
        self._load_available_unification_protocols()

        unification_protocol_name = unification_protocol_name.lower()

        if not self.available_unification_protocols.has_key( unification_protocol_name ):
            raise ValueError("ERROR. Trying to get an unexisting unification protocol")

        if unification_protocol_name == "no unification":
            raise ValueError("Default unification protocol does not have any atom information")

        unificationProtocolObj = self.available_unification_protocols[unification_protocol_name]

        unification_protocol_atoms = unificationProtocolObj.get_unification_atom_elements()
        
        #for atom in atoms:
        #    atom.get_externalDatabaseID_A()
        #    atom.get_externalDatabaseID_B()
        #    atom.get_external_attribute_list()

        return unification_protocol_atoms 


    def get_equivalent_user_entities_from_list(self, userEntitiesList, attribute, protocol_id):
        """
        """
        table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute.lower()]["table"].get_table_name()
        protocol_table = "%s%s" %(self.biana_database.USER_ENTITY_TABLE,protocol_id)
        
        data = self.db.select_db_content( self.db._get_select_sql_query( tables = [(protocol_table,"u1"),
                                                                                   (protocol_table,"u2"),
                                                                                   (table,"e1"),
                                                                                   (table,"e2")],
                                                                         columns = ["u1.userEntityID", "u2.userEntityID"],
                                                                         join_conditions = [("e1.%s" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute.lower()]["fields"]["value"].get_field_name(),
                                                                                             "=",
                                                                                             "e2.%s" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute.lower()]["fields"]["value"].get_field_name()),
                                                                                            ("u1.externalEntityID","<","u2.externalEntityID"),
                                                                                            ("u1.externalEntityID","=","e1.externalEntityID"),
                                                                                            ("u2.externalEntityID","=","e2.externalEntityID")],
                                                                         fixed_conditions = [("u1.userEntityID","IN","(%s)" %(",".join([ str(x) for x in userEntitiesList])),None),
                                                                                             ("u2.userEntityID","IN","(%s)" %(",".join([ str(x) for x in userEntitiesList])),None)] ),
                                          answer_mode = "raw", remove_duplicates="no" )

        import biana.ext.biana_networkx as networkx
        
        g = networkx.Graph()
        g.add_nodes_from(userEntitiesList)

        for current_data in data:
            g.add_edge(current_data[0],current_data[1])
            
        return networkx.connected_components(g)

        

    def _get_equivalent_external_entities(self, unification_atom_element):
        """
        returns the query to obtain the list of equivalent
        """
        
        actual_atom_element = unification_atom_element
        attribute_list = actual_atom_element.get_external_attribute_list()

        tables = [(self.biana_database.EXTERNAL_ENTITY_TABLE,"e1"),
                  (self.biana_database.EXTERNAL_ENTITY_TABLE,"e2")]

        fixed_conditions = []
        join_conditions = []

	# Skip (delay processing till all user entities created) external databases providing promiscuous data
	db_id_A = actual_atom_element.get_externalDatabaseID_A()
	db_id_B = actual_atom_element.get_externalDatabaseID_B()
	if (self.get_external_database(database_id = db_id_A).get_promiscuity() or self.get_external_database(database_id = db_id_B).get_promiscuity()):
	    return None

	# If only none of the external database is promiscuous
        for attribute_index in xrange(len(attribute_list)):

            actual_attribute = attribute_list[attribute_index]

            tables.append( (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[actual_attribute.lower()],"a%s" %(attribute_index)) )
            tables.append( (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[actual_attribute.lower()],"b%s" %(attribute_index)) )
            
	    # To check... Not being used (field conditions are like restrictions on attribute values, i.e. { "type": "unique" })
            fixed_conditions.extend( [ ("a%s.%s" %(attribute_index,self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[actual_attribute]["fields"][x[0].lower()]),"=",x[1])
                                       for x in actual_atom_element.get_field_conditions_A() ] )

            fixed_conditions.extend( [ ("b%s.%s" %(attribute_index,self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[actual_attribute]["fields"][x[0].lower()]),"=",x[1])
                                       for x in actual_atom_element.get_field_conditions_B() ] )

	    # do not unify based on previous attributes (emre), may decrease unifying efficiency as well
	    fixed_conditions.append( ("a%s.type" % attribute_index,"!=","previous") ) 
	    fixed_conditions.append( ("b%s.type" % attribute_index,"!=","previous") )

            join_conditions.extend( [ ("a%s.value" %attribute_index,
                                       "=",
                                       "b%s.value" %attribute_index)
                                      for x in actual_atom_element.get_field_cross_references() ] )

            # Special treatment for PDB: in PDB, pdb code and chain must be shared.
            if actual_attribute.lower()=="pdb":
                join_conditions.extend( [ ("a%s.chain" %attribute_index,
                                           "=",
                                           "b%s.chain" %attribute_index)
                                          for x in actual_atom_element.get_field_cross_references() ] )

            join_conditions.append( ("e1.externalEntityID","=","a%s.externalEntityID" %attribute_index) )
            join_conditions.append( ("e2.externalEntityID","=","b%s.externalEntityID" %attribute_index) )

	# Consider equal only external entities have the same type
        join_conditions.append( ("e1.type","=","e2.type") ) 

        fixed_conditions.append( ("e1.externalDatabaseID","=",actual_atom_element.get_externalDatabaseID_A()) )
        fixed_conditions.append( ("e2.externalDatabaseID","=",actual_atom_element.get_externalDatabaseID_B()) )

        fixed_conditions.append( ("e1.type","!=","relation") )        #JAVI RECENTLY ADDED. MAY DECREASE UNIFYING EFICIENCY... ADDED TO AVOID ADDING relations to unification. It can also be filtered in a posterior step
        fixed_conditions.append( ("e2.type","!=","relation") )
        
        if( actual_atom_element.get_externalDatabaseID_A() == actual_atom_element.get_externalDatabaseID_B() ):
            join_conditions.append( ("e1.externalEntityID","<","e2.externalEntityID") )    # JAVI QUESTION: Does this affect the unification? Is it correct???


        return self.db._get_select_sql_query( tables = tables,
                                              columns = ["e1.externalEntityID","e2.externalEntityID"],
                                              fixed_conditions = fixed_conditions,                                                    
                                              join_conditions = join_conditions )

    
    def get_external_entities_dict(self, externalEntityIdsList, attribute_list=[], relation_attribute_list=[], participant_attribute_list=[], useTransferAttributes=True):
        """
        Returns a dict of external Entity Objects with the attributes specified in the "attribute_identifier_list"

        The key in the dictionary corresponds to the external Entity ID

        The external entity can be of any type (included relations)

        "attribute_field_restrictions" is to restrict the attributes by additional fields. Sintaxis: [(attribute_identifier, field, value)]
        """

        if len(externalEntityIdsList)==0:
            return {}

        eE_id_str_list = ", ".join([ str(x) for x in externalEntityIdsList ])

        # Get the basic information for the external entities
        data = self.db.select_db_content(self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name()],
                                                                        columns = ["%s.externalEntityID" %self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name(),
                                                                                   "%s.externalDatabaseID" %self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name(),
                                                                                   "%s.type" %self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name()],
                                                                        fixed_conditions = [("%s.externalEntityID" %self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name(),
                                                                                             "IN",
                                                                                             "(%s)" %eE_id_str_list, None)] ),
                                         answer_mode="raw", remove_duplicates="no")
        eE_dict = {}

        eEr_id_list = []
        
        for current_data in data:

            if current_data[2]=="relation":
                # Get the relation type
                relation_type = self.db.select_db_content( self.db._get_select_sql_query( tables=[self.biana_database.EXTERNAL_ENTITY_RELATION_TABLE],
                                                                                          columns = ["type"],
                                                                                          fixed_conditions = [("externalEntityRelationID","=",current_data[0])] ),
                                                           answer_mode = "single")

                relationObj = BianaObjects.ExternalEntityRelation( id = current_data[0],
                                                                   source_database = current_data[1],
                                                                   relation_type = relation_type )

                eE_dict[current_data[0]] = relationObj

                # Add the rest of participants of this relation
                participants = self.db.select_db_content( self.db._get_select_sql_query( tables=[self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE],
                                                                                         columns = ["externalEntityID"],
                                                                                         fixed_conditions = [("externalEntityRelationID","=",current_data[0])] ),
                                                          answer_mode = "raw" )

                [ relationObj.add_participant( externalEntityID=x[0] ) for x in participants ]

                eEr_id_list.append(str(current_data[0]))

            else:
                eE_dict[current_data[0]] = BianaObjects.ExternalEntity( id = current_data[0],
                                                                        source_database = current_data[1],
                                                                        type = current_data[2] )



        # Prepare restrictions
        #attribute_field_restrictions_dict = dict([(x[0].lower(),(x[1],x[2])) for x in attribute_field_restrictions ])

        # Now, we can merge simple eE with eEr, as they are treated in the same way
        attributes = set([ x.lower() for x in attribute_list ])
        attributes.update([ x.lower() for x in relation_attribute_list ])

        #for current_attribute in set([ x.lower() for x in attribute_list ]):
        for current_attribute in attributes:

            if current_attribute == "proteinsequence":
                # Get the sequence and sequence ID
                tables = [self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute],
                          self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"]]
                columns = ["externalEntityID","UNCOMPRESS(sequence)"]
                data = self.db.select_db_content(self.db._get_select_sql_query( tables = tables,
                                                                                columns=columns,
                                                                                fixed_conditions = [("externalEntityID",
                                                                                                     "IN",
                                                                                                     "(%s)" %eE_id_str_list, None)],
                                                                                join_conditions = [("%s.value" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute],
                                                                                                    "=",
                                                                                                    "%s.sequenceMD5" %self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"])]),
                                                 answer_mode = "raw" )
                for current_sequence in data:
                    eE_dict[current_sequence[0]].add_attribute(BianaObjects.ExternalEntityAttribute(attribute_identifier=current_attribute,
                                                                                                    value=BianaObjects.ProteinSequence(sequence=current_sequence[1])))

                continue

            if current_attribute == "nucleotidesequence":
                # Get the sequence and sequence ID
                tables = [self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute],
                          self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"]]
                columns = ["externalEntityID","sequenceType","UNCOMPRESS(sequence)"]
                data = self.db.select_db_content(self.db._get_select_sql_query( tables = tables,
                                                                                columns=columns,
                                                                                fixed_conditions = [("externalEntityID",
                                                                                                     "IN",
                                                                                                     "(%s)" %eE_id_str_list, None)],
                                                                                join_conditions = [("%s.value" %self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute],
                                                                                                    "=",
                                                                                                    "%s.sequenceMD5" %self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["nucleotideSequence"])]),
                                                 answer_mode = "raw" )
                for current_sequence in data:
                    if( current_sequence[1]=="dna" ):
                        eE_dict[current_data[0]].add_attribute(BianaObjects.ExternalEntityAttribute(attribute_identifier=current_attribute,
                                                                                                    value=BianaObjects.DNASequence(sequence=current_sequence[2])))
                    elif( current_sequence[1]=="rna" ):
                        eE_dict[current_data[0]].add_attribute(BianaObjects.ExternalEntityAttribute(attribute_identifier=current_attribute,
                                                                                                    value=BianaObjects.RNASequence(sequence=current_sequence[2])))
                        
                continue
            
            try:
                table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute]

                columns = ["%s.externalEntityID" %table,
                           "%s.value" %table]
                           #"%s.type" %table]

                fixed_conditions = [("%s.externalEntityID" %table,
                                     "IN",
                                     "(%s)" %eE_id_str_list, None)]
                
                data = self.db.select_db_content(self.db._get_select_sql_query( tables=[table],
                                                                                columns=columns,
                                                                                fixed_conditions = fixed_conditions),
                                                 answer_mode="raw", remove_duplicates="no")

                for current_data in data:
                    eE_dict[current_data[0]].add_attribute(BianaObjects.ExternalEntityAttribute(attribute_identifier=str(current_attribute).replace("\n"," "),   
                                                                                                value=str(current_data[1]).replace("\n"," ")) )  # Changed to avoid inserting new lines


                ### ADD TRANSFERRED ATTRIBUTES ###

                if self._is_transferred_attribute( current_attribute ) and useTransferAttributes is True:

                    for current_transfer in self.transferred_attributes[current_attribute.lower()]:

                        tables = [ (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute],"transferred"),
                                   (self._get_key_attribute_table_name( key_id = self.key_attribute_ids[(current_transfer[0],current_transfer[1])] ),"key_attr"),
                                   (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_transfer[1]],"key_attr2") ]

                        fixed_conditions = [("key_attr2.externalEntityID",
                                             "IN",
                                             "(%s)" %eE_id_str_list, None)]

                        join_conditions = [("key_attr.value","=","key_attr2.value"),
                                           ("transferred.externalEntityID","=","key_attr.externalEntityID")]

                        columns  = ["key_attr2.externalEntityID","transferred.value","transferred.type"]

                                   
                        data = self.db.select_db_content(self.db._get_select_sql_query( tables=tables,
                                                                                        columns=columns,
                                                                                        fixed_conditions = fixed_conditions,
                                                                                        join_conditions = join_conditions),
                                                         answer_mode="raw", remove_duplicates="no")

                        for current_data in data:
                            eE_dict[current_data[0]].add_attribute(BianaObjects.ExternalEntityAttribute(attribute_identifier=current_attribute,
                                                                                                        value=current_data[1],
                                                                                                        type="transferred_%s" %(current_data[2]) ))
                                   


            except:
                traceback.print_exc()
                sys.stderr.write("Attribute %s is not found in available attributes\n" %(current_attribute))


            for current_attribute in set([ x.lower() for x in participant_attribute_list]):

                try:
                    for current_eErID in eEr_id_list:
                        current_eEr_obj = eE_dict[int(current_eErID)]
                        data = self.db.select_db_content( self.db._get_select_sql_query( tables=[self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE,
                                                                                                 self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT[current_attribute]],
                                                                                                 columns = ["externalEntityID","value"],
                                                                                                 join_conditions = [("%s.externalEntityRelationParticipantID" %self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE,
                                                                                                                     "=",
                                                                                                                     "%s.externalEntityRelationParticipantID" %self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT[current_attribute])],
                                                                                                 fixed_conditions = [("externalEntityID","IN","(%s)" %",".join([ str(x) for x in current_eEr_obj.get_participant_external_entity_ids_list()]),None)]),
                                                                  answer_mode = "raw" )

                        for current_data in data:
                            current_eEr_obj.add_participant_attribute( externalEntityID = current_data[0],
                                                                       participantAttribute = BianaObjects.ExternalEntityRelationParticipantAttribute( attribute_identifier = current_attribute,
                                                                                                                                                       value = current_data[1].replace("\n", " ") ))

                            
                except:
                    traceback.print_exc()
                    sys.stderr.write("Attribute %s is not found in available relation participant attributes\n" %(current_attribute))

                    
        return eE_dict


    def get_user_entity_attributes(self, unification_protocol_name, listUserEntityID, attribute_identifier):
        """
        Returns a dictionary with { userEntityID: list of attributes }

        This method is intended to be faster for getting user entity attributes than getting external entity objects and its attribute objects

        Numerical values are stored as strings!!!
        """

        if debug_web == 1:
            import os
            debug_file_fd = file("/usr/local/apache2/htdocs/webiana/sessions/debug_file_biana.txt", "w")
            os.chmod("/usr/local/apache2/htdocs/webiana/sessions/debug_file_biana.txt", 0777)
            debug_file_fd.write("attribute_identifier is %s" %(attribute_identifier))
            debug_file_fd.write("external entity attribute tables dict has keys %s" %(self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT.keys()))

        attr_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier]
        unif_table = self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)

        if attribute_identifier.lower()=="proteinsequence":
            query = self.db._get_select_sql_query( tables = [attr_table,unif_table,
                                                             self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"]],
                                                   columns = ["userEntityID","UNCOMPRESS(sequence) AS seq"],
                                                   join_conditions = [("userEntityID","IN","(%s)" %",".join(map(str,listUserEntityID))),
                                                                      ("%s.externalEntityID" %attr_table,"=","%s.externalEntityID" %unif_table),
                                                                      ("%s.value" %attr_table,"=","%s.sequenceMD5" %self.biana_database.EXTERNAL_ATTRIBUTES_DESCRIPTION_TABLES["proteinSequence"])],
                                                   group_conditions = ["userEntityID","seq"] )
        else:
            query = self.db._get_select_sql_query( tables = [attr_table,unif_table],
                                                   columns = ["userEntityID","value"],
                                                   join_conditions = [("userEntityID","IN","(%s)" %",".join(map(str,listUserEntityID))),
                                                                      ("%s.externalEntityID" %attr_table,"=","%s.externalEntityID" %unif_table)],
                                                   group_conditions = ["userEntityID","value"] )

        data = self.db.select_db_content( query, answer_mode = "raw" )

        return_dict = {}

        for current_data in data:
            return_dict.setdefault(current_data[0],[]).append(str(current_data[1]).replace("\n"," ")) # CHANGED TO REMOVE new lines in attributes
            

        ### ADD TRANSFERRED ATTRIBUTES ###
        if self._is_transferred_attribute( attribute_identifier ):
            
            for current_transfer in self.transferred_attributes[attribute_identifier.lower()]:

                tables = [ unif_table,
                           (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_identifier],"transferred"),
                           (self._get_key_attribute_table_name( key_id = self.key_attribute_ids[(current_transfer[0],current_transfer[1])] ),"key1"),
                           (self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_transfer[1]],"key2") ]

                columns = ["userEntityID","transferred.value"]

                join_conditions = [("userEntityID","IN","(%s)" %",".join(map(str,listUserEntityID))),
                                   ("%s.externalEntityID" %unif_table,"=","key2.externalEntityID"),
                                   ("key1.value","=","key2.value"),
                                   ("transferred.externalEntityID","=","key1.externalEntityID")]

                query = self.db._get_select_sql_query( tables = tables,
                                                       columns = columns,
                                                       join_conditions = join_conditions )

                data = self.db.select_db_content( query, answer_mode = "raw" )

                for current_data in data:
                    return_dict.setdefault(current_data[0],[]).append(str(current_data[1]))

        return return_dict

    def get_user_entity_type(self, unification_protocol_name, user_entity_ids):
        """
        Gets types associated with given user entity ids from database

        Returns a dictionary of user entity id and type
        """
        
        self._load_available_unification_protocols()

        #fixed_conditions = [("u.externalEntityID","IN","(%s)" %",".join(map(str,external_entity_IDs_list)),None)]
	fixed_conditions = [ ("u.userEntityID","IN","(%s)" %",".join(map(str,user_entity_ids)),None) ]
	tables = [ "%s u" % self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name),
		   "%s e" % self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name() ]
	columns = [ "u.userEntityID", "e.type" ]
	join_conditions = [ ("u.externalEntityID", "=", "e.externalEntityID") ]
	group_conditions = [ "u.userEntityID", "e.type" ]

	query = self.db._get_select_sql_query( tables = tables,
					       columns = columns,
					       fixed_conditions = fixed_conditions,
					       join_conditions = join_conditions,
					       group_conditions = group_conditions)

	userEntityId_type_tuples = self.db.select_db_content(query,
						     answer_mode = "raw" )
	# gene has type precedence over protein
	uEId_to_type = {}
	for uEId, type in userEntityId_type_tuples:
	    if uEId_to_type.setdefault(uEId, type) != "gene":
		if type == "gene":
		    uEId_to_type[uEId] = type

	return uEId_to_type



    def _get_list_eE_for_uE(self, unification_protocol_name, userEntityID):

        return self.db.select_db_content( self.db._get_select_sql_query( tables = [self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)],
                                                                         columns = ["externalEntityID"],
                                                                         fixed_conditions = [("userEntityID","=",userEntityID)] ),
                                          answer_mode = "list" )


    def get_user_entity_relations_by_sharing_attributes(self, unification_protocol_name, userEntityID_list, listAttributes, limit_to_userEntityID_list=False, attribute_restrictions=[], negative_attribute_restrictions=[], ontology_expansion_level=None):
        """
        Returns a list of relations between user entities, in which the share attributes listed in listAttributes
        
        "listAttributes" is a list of lists of external entity attributes. Each sublist is a "sharing" restriction

        "expand_ontology_attributes": Boolean to specify if ontology attributes should be expanded to lower levels

        "ontology_expansion_level": Dictionary to specigy the category level to which should be considered equivalent two ontology attributes (for example, family level in scop)
        """

        # TODO: APPLY "ontology_expansion_level" attribute

        if len(userEntityID_list)==0:
            return []

        equivalent_uE_pair_queries = []

        for current_attribute_group in listAttributes:

            equivalent_uE_pair_queries.append(self._get_nested_queries_for_getting_equivalent_uE_by_sharing_attributes( unification_protocol_name = unification_protocol_name,
                                                                                                                        attribute_list = current_attribute_group,
                                                                                                                        user_entity_ids_list = userEntityID_list,
                                                                                                                        restrict_to_user_entity_ids_list = limit_to_userEntityID_list))
        
        union_query = self.db._get_union_queries( equivalent_uE_pair_queries )

        
        # Apply restrictions

        restricted_query = self._apply_restrictions_to_query( query = union_query, 
                                                              unification_protocol_name = unification_protocol_name, 
                                                              attribute_restrictions = attribute_restrictions,
                                                              column_name_to_restrict="userEntityID2" )

        negatively_restricted_query = self._apply_negative_restrictions_to_query( query = restricted_query,
                                                                                  unification_protocol_name = unification_protocol_name, 
                                                                                  column_name_to_restrict="userEntityID2",
                                                                                  negative_attribute_restrictions = negative_attribute_restrictions )


        #print negatively_restricted_query
        

        results = self.db.select_db_content( negatively_restricted_query,
                                             answer_mode = "raw" )

        if results:
	    uE1s, uE2s = zip(*results)
	    uE2_types_dict = self.get_user_entity_type(unification_protocol_name, uE2s)
            #new_results = [ (x[0], x[1][0], x[1][1]) for x in zip(uE1s, tuples) ]
            new_results = []
            for current_result in results:
                new_results.append((current_result[0],current_result[1],uE2_types_dict[current_result[1]]))
        else:
            new_results = []

        return new_results

    def _get_list_external_entities_for_user_entities(self, unification_protocol_name, userEntityID_list):
        
        tables = [ "%s u" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name) ) ]
        columns = ["externalEntityID"]
        join_conditions = [("u.userEntityID","IN","(%s)" %",".join(map(str,userEntityID_list)))]

        return self.db.select_db_content( self.db._get_select_sql_query( tables = tables,
                                                                         columns = columns,
                                                                         join_conditions = join_conditions ),
                                          answer_mode = "list" )

    def _get_list_user_entities_for_external_entities(self, unification_protocol_name, externalEntityID_list):
        
        tables = [ "%s u" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name) ) ]
        columns = ["userEntityID"]
        join_conditions = [("u.externalEntityID","IN","(%s)" %",".join(map(str,externalEntityID_list)))]

        return self.db.select_db_content( self.db._get_select_sql_query( tables = tables,
                                                                         columns = columns,
                                                                         join_conditions = join_conditions ),
                                          answer_mode = "list" )


    def get_default_external_entity_ids(self, externalEntityIDsList):

        # MONDAY: METHOD TO CHANGE !!!!!!!!
        
        if len(externalEntityIDsList)==0:
            return {}

        default_attribute_identifiers = self.db.select_db_content( self.db._get_select_sql_query( tables = [ self.biana_database.EXTERNAL_DATABASE_TABLE, self.biana_database.EXTERNAL_ENTITY_TABLE ],
                                                                                                  columns = [ "%s.externalEntityID" %(self.biana_database.EXTERNAL_ENTITY_TABLE),
                                                                                                              "%s.defaultExternalEntityAttribute" %(self.biana_database.EXTERNAL_DATABASE_TABLE) ],
                                                                                                  
                                                                                                  join_conditions = [("%s.externalDatabaseID" %(self.biana_database.EXTERNAL_DATABASE_TABLE),"=","%s.externalDatabaseID" %(self.biana_database.EXTERNAL_ENTITY_TABLE)),
                                                                                                                     ("%s.externalEntityID" %(self.biana_database.EXTERNAL_ENTITY_TABLE),"IN","(%s)" %(",".join(map(str,externalEntityIDsList))))] ),
                                                                   answer_mode = "raw" )
        
        attribute_identifiers_dict = {}
        for externalEntityID, defaultAttribute in default_attribute_identifiers:
            attribute_identifiers_dict.setdefault(defaultAttribute,[]).append(externalEntityID)

        return_dict = {}

        # Get the id for external entities
        for current_attribute, list_externalEntityIds in attribute_identifiers_dict.iteritems():
            if current_attribute is None:
                for current_eEid in list_externalEntityIds:
                    return_dict[current_eEid] = ""  # What should I do here?
            else:
                tables = [self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[current_attribute.lower()]]
                join_conditions = [("externalEntityID","IN","(%s)" %",".join(map(str,list_externalEntityIds)))]
                columns = ["externalEntityID","value"]

                #print self.db._get_select_sql_query( tables = tables,
                #                                     columns = columns,
                #                                     join_conditions = join_conditions )

                raw_data = self.db.select_db_content( self.db._get_select_sql_query( tables = tables,
                                                                                     columns = columns,
                                                                                     join_conditions = join_conditions,
										     fixed_conditions = [("type","=","unique")],
                                                                                     group_conditions = ["externalEntityID"]),
                                                      answer_mode = "raw" )

                return_dict.update( dict( [ (x[0],"%s: %s" %(current_attribute,x[1])) for x in raw_data ] ) )

                del raw_data

                #return_dict.update(dict(self.db.select_db_content( self.db._get_select_sql_query( tables = tables,
                #                                                                                  columns = columns,
                #                                                                                  join_conditions = join_conditions,
                #                                                                                  group_conditions = ["externalEntityID"]),
                #                                                   answer_mode = "raw" )))

        return return_dict


    def get_relations_hierarchy(self, externalEntityRelationIDs):
        """
        Returns a list of tuples (relationID1, relationID2) in where relationID1 is a child of relationID2
        """

        tables = ["%s p" %self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE]
        
        columns = ["p.externalEntityID","p.externalEntityRelationID"]
        join_conditions = [("p.externalEntityID","IN","(%s)" %",".join(map(str,externalEntityRelationIDs))),
                           ("p.externalEntityRelationID","IN","(%s)" %",".join(map(str,externalEntityRelationIDs)))]
        
        return self.db.select_db_content( self.db._get_select_sql_query( tables = tables,
                                                                         columns = columns,
                                                                         join_conditions = join_conditions ),
                                          answer_mode = "raw" )


    def _update_relations_hierarchy(self):
        """
        Precalculates and stores in database relations hierachy
        """

        return "DEPRECATED"

        hierarchy_set = set()  # Set that contains tuples (externalEntityRelationID, externalEntityID), where externalEntityID is included in externalEntityRelationID. externalEntityID can also be a relation

        # Get all the information about relations and participants
        eE_eERid_list = self.db.select_db_content( self.db._get_select_sql_query( tables = [self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE],
                                                                                  columns = ["externalEntityID","externalEntityRelationID"]),
                                                        answer_mode = "raw" )
        eE_eERid_dict = {}

        for current_pair in eE_eERid_list:
            eE_eERid_dict.setdefault(current_pair[0],[]).append(current_pair[1])


        def get_all_parents(eEid):
            parents = set()
            if( eE_eERid_dict.has_key(eEid) ):
                parents.update(eE_eERid_dict[eEid])
                for current_eEid in eE_eERid_dict[eEid]:
                    parents.update(get_all_parents(current_eEid))
            return parents


        for eEid in eE_eERid_dict.keys():

            for current_parent in get_all_parents(eEid):
                hierarchy_set.add((current_parent,eEid))


        self.db.set_lock_tables(True)

        self.db._disable_indices([self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE])

        for current_eERid, current_eEid in hierarchy_set:
            self.db.insert_db_content( self.db._get_insert_sql_query( table = self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE,
                                                                      column_values = (("externalEntityRelationParticipantID", self._get_new_external_entity_relation_participant_id()), 
                                                                                       (self.biana_database.external_entity_relation_id_col, current_eERid),
                                                                                       (self.biana_database.externalEntityID_col, current_eEid )),
                                                                      use_buffer = True ))
        self.db._empty_buffer()
        self.db._enable_indices([self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE])

        self.db.set_lock_tables(False)
        return


    def OLDget_user_entity_relations(self, unification_protocol_name, userEntityID_list, attribute_restrictions = [], negative_attribute_restrictions = [], listRelationType=[], dictRelationAttributeRestriction={}, use_self_relations=True, limit_to_userEntityID_list=False, use_nested_relations=True):
	"""
        #! Nasty trick (not used anymore): to not to join externalEntityType table in the query in case sth is interacting with itself (see use_self_relations) there is an internal type called "self_type" which should be handled carefully whenever this method is called and type is going to be used
	"""

        if len(userEntityID_list)==0:
            return []

        if use_nested_relations:
            PARTICIPANT_TABLE = self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE.get_table_name()
        else:
            PARTICIPANT_TABLE = self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE.get_table_name()
        
        tables = ["%s u1" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                  "%s p1" %PARTICIPANT_TABLE,
                  "%s u2" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                  "%s p2" %PARTICIPANT_TABLE,
                  "%s r" %(self.biana_database.EXTERNAL_ENTITY_RELATION_TABLE.get_table_name()),
                  "%s e" %(self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name())]
        
        columns = [ ("u1.userEntityID","userEntityID1"),
                    ("u2.userEntityID","userEntityID2"),
                    ("r.externalEntityRelationID", "externalEntityRelationID"),
                    ("r.type", "type"),
                    ("e.type", "etype") ]

        #group_conditions = ["u1.userEntityID","u2.userEntityID","r.externalEntityRelationID, r.type, e.type"]

        fixed_conditions = [("u1.userEntityID","IN","(%s)" %",".join([ str(x) for x in userEntityID_list]),None)]
        join_conditions = [("p1.externalEntityID","=","u1.externalEntityID"),
                           ("p1.externalEntityRelationID","=","p2.externalEntityRelationID"),
                           ("p2.externalEntityID","=","u2.externalEntityID"),
                           ("p1.externalEntityID","!=","p2.externalEntityID"),
                           ("p1.externalEntityRelationID","=","r.externalEntityRelationID"),
                           ("p2.externalEntityID", "=", "e.externalEntityID")]

        if limit_to_userEntityID_list is True:
            join_conditions.append(("u2.userEntityID","IN","(%s)" %",".join([ str(x) for x in userEntityID_list ]),None))

        # Apply relation attribute restrictions
        num = 1
        for attribute_name, values in dictRelationAttributeRestriction.iteritems():

            current_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_name.lower()]
            tables.append( "%s AS EERA%s" %(current_table.get_table_name(),num) )
            if BianaObjects.ExternalEntityAttribute.isNumericAttribute(attribute_name, self.biana_database) or BianaObjects.ExternalEntityAttribute.isSpecialAttribute(attribute_name, self.biana_database):
                regex = re.compile("([><=]*)([\d\.]+)")
                for current_value in values:
                    m = regex.match(current_value)
                    if m:
                        join_conditions.append(("EERA%s.value" %num,m.group(1),m.group(2)))
                    else:
                        join_conditions.append(("EERA%s.value" %num,"=","\"%s\"" %current_value))
            else:
                join_conditions.append(("EERA%s.value" %num,"IN","(\"%s\")" %("\",\"".join([ str(x) for x in values]))))
            join_conditions.append( ("r.externalEntityRelationID","=","EERA%s.externalEntityID" %num ) )
            num += 1


        if len(listRelationType) > 0:
            fixed_conditions.append(("r.type","IN","(\"%s\")" % "\",\"".join([ relationType for relationType in listRelationType]),None ))                
                
        # Get the interacting user entities
        query = self.db._get_select_sql_query( tables = tables,
                                               columns = columns,
                                               fixed_conditions = fixed_conditions, 
                                               join_conditions = join_conditions,
                                               distinct_columns = True )
                                               # group_conditions = group_conditions )


        query = self._apply_restrictions_to_query( query = query,
                                                   unification_protocol_name = unification_protocol_name,
                                                   attribute_restrictions= attribute_restrictions,
                                                   column_name_to_restrict="userEntityID2" )
                                                   #column_name_to_restrict="userEntityID" )

        query = self._apply_negative_restrictions_to_query( query = query,
                                                            unification_protocol_name = unification_protocol_name,
                                                            negative_attribute_restrictions = negative_attribute_restrictions,
                                                            column_name_to_restrict="userEntityID2" )
                                                            #column_name_to_restrict="userEntityID" )


        #print query
        
        interacting_uE = list(self.db.select_db_content( query, answer_mode = "raw" ))

        if( use_self_relations is True ):

            # Get self interactions
            fixed_conditions = [("u1.userEntityID","IN","(%s)" %",".join([ str(x) for x in userEntityID_list]),None)]

            if len(listRelationType) > 0:
                fixed_conditions.append(("r.type","IN","(\"%s\")" % "\",\"".join([ relationType for relationType in listRelationType]),None ))

            tables = ["%s u1" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                      "%s p1" %PARTICIPANT_TABLE,
                      "%s c" %self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT["cardinality"].get_table_name(),
                      "%s r" %(self.biana_database.EXTERNAL_ENTITY_RELATION_TABLE.get_table_name()),
		      "%s e" %(self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name())]

            columns = ["u1.userEntityID","r.externalEntityRelationID", ("r.type", "type"), ("e.type", "etype") ]


            join_conditions = [("p1.externalEntityID","=","u1.externalEntityID"),
                               ("c.externalEntityRelationParticipantID","=","p1.externalEntityRelationParticipantID"),
                               ("c.value",">",1),
                               ("r.externalEntityRelationID","=","p1.externalEntityRelationID"), 
			       ("p1.externalEntityID", "=", "e.externalEntityID")]


            query = self.db._get_select_sql_query( tables = tables,
                                                   columns = columns,
                                                   fixed_conditions =  fixed_conditions,
                                                   join_conditions = join_conditions )
            
            query = self._apply_relation_restrictions_to_query( query = query,
                                                                attribute_restrictions_dict = dictRelationAttributeRestriction, 
                                                                column_name_to_restrict="externalEntityRelationID" )
            

            query = self._apply_restrictions_to_query( query = query,
                                                       unification_protocol_name = unification_protocol_name,
                                                       attribute_restrictions= attribute_restrictions,
                                                       column_name_to_restrict="userEntityID" )
            
            query = self._apply_negative_restrictions_to_query( query = query,
                                                                unification_protocol_name = unification_protocol_name,
                                                                negative_attribute_restrictions = negative_attribute_restrictions,
                                                                column_name_to_restrict="userEntityID" )

            #print query

            #interacting_uE.extend([ (x,x,y,z,"self_type") for x,y,z in self.db.select_db_content( query, answer_mode="raw" )])
            interacting_uE.extend([ (x,x,y,z,t) for x,y,z,t in self.db.select_db_content( query, answer_mode="raw" )])

        #print len(interacting_uE)

        return interacting_uE


    def NEWget_user_entity_relations(self, unification_protocol_name, userEntityID_list, attribute_restrictions = [], negative_attribute_restrictions = [], listRelationType=[], dictRelationAttributeRestriction={}, use_self_relations=True, limit_to_userEntityID_list=False, use_nested_relations=True):
        """
        Returns the list of relations where the userEntity is involved
        use_nested_relations ==> Include relations greater than a nested level
        """

        if len(userEntityID_list)==0:
            return []

        if use_nested_relations:
            PARTICIPANT_TABLE = self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE.get_table_name()
        else:
            PARTICIPANT_TABLE = self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE.get_table_name()
  
        # First, get the relations
        tables = [ "%s u1" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                   "%s p1" %PARTICIPANT_TABLE,
                   "%s r" %(self.biana_database.EXTERNAL_ENTITY_RELATION_TABLE.get_table_name())]

        columns = [ ("u1.userEntityID", "userEntityID"), 
                    ("r.externalEntityRelationID", "externalEntityRelationID"),
                    ("r.type", "type"),
                    ("p1.externalEntityID", "externalEntityID")]

        join_conditions = [("u1.userEntityID","IN","(%s)" %",".join([ str(x) for x in userEntityID_list])),
                           ("p1.externalEntityID","=","u1.externalEntityID"),
                           ("p1.externalEntityRelationID","=","r.externalEntityRelationID")]

        if len(listRelationType) > 0:
            join_conditions.append(("r.type","IN","(\"%s\")" % "\",\"".join([ relationType for relationType in listRelationType])))

       # Apply relation attribute restrictions                                                                                                                                                                  
        num = 1
        for attribute_name, values in dictRelationAttributeRestriction.iteritems():

            current_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_name.lower()]
            tables.append( "%s AS EERA%s" %(current_table.get_table_name(),num) )
            if BianaObjects.ExternalEntityAttribute.isNumericAttribute(attribute_name, self.biana_database) or BianaObjects.ExternalEntityAttribute.isSpecialAttribute(attribute_name, self.biana_database):
                regex = re.compile("([><=]*)([\d\.]+)")
                for current_value in values:
                    m = regex.match(current_value)
                    if m:
                        join_conditions.append(("EERA%s.value" %num,m.group(1),m.group(2)))
                    else:
                        join_conditions.append(("EERA%s.value" %num,"=","\"%s\"" %current_value))
            else:
                join_conditions.append(("EERA%s.value" %num,"IN","(\"%s\")" %("\",\"".join([ str(x) for x in values]))))
            join_conditions.append( ("r.externalEntityRelationID","=","EERA%s.externalEntityID" %num ) )
            num += 1

        query = self.db._get_select_sql_query( tables = tables,
                                               columns = columns,
                                               join_conditions = join_conditions )

        # Then, get the u2 participants
        tables = [ ("(%s)" %query, "query1"),
                   "%s p2" %PARTICIPANT_TABLE,
                   "%s u2" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                   "%s e" %(self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name()) ]

        columns = [ "query1.userEntityID", 
                    ("u2.userEntityID", "userEntityID2"), 
                    "query1.externalEntityRelationID", 
		    "query1.type",
		    ("e.type", "etype") ]

        join_conditions = [ ("query1.externalEntityRelationID","=","p2.externalEntityRelationID"),
                            ("p2.externalEntityID","=","u2.externalEntityID"),
                            ("query1.externalEntityID","!=","p2.externalEntityID"),
                            ("p2.externalEntityID", "=", "e.externalEntityID")]

        if limit_to_userEntityID_list is True:
            join_conditions.append(("u2.userEntityID","IN","(%s)" %",".join([ str(x) for x in userEntityID_list ])))


        #query = self.db._get_select_sql_query( tables = tables,
        #                                       columns = columns,
        #                                       join_conditions = join_conditions )
        

        #if limit_to_userEntityID_list is True:
        #    join_conditions.append(("u2.userEntityID","IN","(%s)" %",".join([ str(x) for x in userEntityID_list ]),None))

        # Get the interacting user entities                                                                                                                                                                                                              
        query = self.db._get_select_sql_query( tables = tables,
                                               columns = columns,
                                               join_conditions = join_conditions,
                                               distinct_columns = True )

        query = self._apply_restrictions_to_query( query = query,
                                                   unification_protocol_name = unification_protocol_name,
                                                   attribute_restrictions= attribute_restrictions,
                                                   column_name_to_restrict="userEntityID2" )

        query = self._apply_negative_restrictions_to_query( query = query,
                                                            unification_protocol_name = unification_protocol_name,
                                                            negative_attribute_restrictions = negative_attribute_restrictions,
                                                            column_name_to_restrict="userEntityID2" )

        interacting_uE = list(self.db.select_db_content( query, answer_mode = "raw" ))

        return interacting_uE


    def get_relations(self, unification_protocol_name, attribute_restrictions = [], negative_attribute_restrictions = [], listRelationType = [], dictRelationAttributeRestriction={}, use_self_relations=True, use_nested_relations=True):
        """
        
        """

        if use_nested_relations:
            PARTICIPANT_TABLE = self.biana_database.EXTENDED_EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE.get_table_name()
        else:
            PARTICIPANT_TABLE = self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_TABLE.get_table_name()
        
        tables = ["%s u1" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                  "%s p1" %PARTICIPANT_TABLE,
                  "%s u2" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                  "%s p2" %PARTICIPANT_TABLE,
                  "%s r" %(self.biana_database.EXTERNAL_ENTITY_RELATION_TABLE.get_table_name()),
                  "%s e" %(self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name())]
        
        columns = [ ("u1.userEntityID","userEntityID1"),
                    ("u2.userEntityID","userEntityID2"),
                    ("r.externalEntityRelationID", "externalEntityRelationID"),
                    ("r.type", "type"),
                    ("e.type", "etype") ]

        join_conditions = [("p1.externalEntityRelationID","=","p2.externalEntityRelationID"),
                           ("p1.externalEntityID","=","u1.externalEntityID"),
                           ("p2.externalEntityID","=","u2.externalEntityID"),
                           ("p1.externalEntityID","!=","p2.externalEntityID"),
                           ("p1.externalEntityRelationID","=","r.externalEntityRelationID"),
                           ("p2.externalEntityID", "=", "e.externalEntityID")]

        fixed_conditions = []

        # Apply relation attribute restrictions
        num = 1
        for attribute_name, values in dictRelationAttributeRestriction.iteritems():

            current_table = self.biana_database.EXTERNAL_ENTITY_ATTRIBUTE_TABLES_DICT[attribute_name.lower()]
            tables.append( "%s AS EERA%s" %(current_table.get_table_name(),num) )
            if BianaObjects.ExternalEntityAttribute.isNumericAttribute(attribute_name, self.biana_database) or BianaObjects.ExternalEntityAttribute.isSpecialAttribute(attribute_name, self.biana_database):
                regex = re.compile("([><=]*)([\d\.]+)")
                for current_value in values:
                    m = regex.match(current_value)
                    if m:
                        join_conditions.append(("EERA%s.value" %num,m.group(1),m.group(2)))
                    else:
                        join_conditions.append(("EERA%s.value" %num,"=","\"%s\"" %current_value))
            else:
                join_conditions.append(("EERA%s.value" %num,"IN","(\"%s\")" %("\",\"".join([ str(x) for x in values]))))
            join_conditions.append( ("r.externalEntityRelationID","=","EERA%s.externalEntityID" %num ) )
            num += 1


        if len(listRelationType) > 0:
            fixed_conditions.append(("r.type","IN","(\"%s\")" % "\",\"".join([ relationType for relationType in listRelationType]),None ))                
                
        # Get the interacting user entities
        query = self.db._get_select_sql_query( tables = tables,
                                               columns = columns,
                                               fixed_conditions = fixed_conditions, 
                                               join_conditions = join_conditions,
                                               distinct_columns = True )

        query = self._apply_restrictions_to_query( query = query,
                                                   unification_protocol_name = unification_protocol_name,
                                                   attribute_restrictions= attribute_restrictions,
                                                   column_name_to_restrict="userEntityID2" )

        query = self._apply_negative_restrictions_to_query( query = query,
                                                            unification_protocol_name = unification_protocol_name,
                                                            negative_attribute_restrictions = negative_attribute_restrictions,
                                                            column_name_to_restrict="userEntityID2" )

        interacting_uE = list(self.db.select_db_content( query, answer_mode = "raw" ))

        if( use_self_relations is True ):

            if len(listRelationType) > 0:
                fixed_conditions.append(("r.type","IN","(\"%s\")" % "\",\"".join([ relationType for relationType in listRelationType]),None ))

            tables = ["%s u1" %(self._get_user_entity_table_name(unification_protocol_name=unification_protocol_name)),
                      "%s p1" %PARTICIPANT_TABLE,
                      "%s c" %self.biana_database.EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TABLES_DICT["cardinality"].get_table_name(),
                      "%s r" %(self.biana_database.EXTERNAL_ENTITY_RELATION_TABLE.get_table_name()),
		      "%s e" %(self.biana_database.EXTERNAL_ENTITY_TABLE.get_table_name())]

            columns = ["u1.userEntityID","r.externalEntityRelationID", ("r.type", "type"), ("e.type", "etype") ]

            join_conditions = [("p1.externalEntityID","=","u1.externalEntityID"),
                               ("c.externalEntityRelationParticipantID","=","p1.externalEntityRelationParticipantID"),
                               ("c.value",">",1),
                               ("r.externalEntityRelationID","=","p1.externalEntityRelationID"), 
			       ("p1.externalEntityID", "=", "e.externalEntityID")]

            query = self.db._get_select_sql_query( tables = tables,
                                                   columns = columns,
                                                   fixed_conditions =  fixed_conditions,
                                                   join_conditions = join_conditions )
            
            query = self._apply_relation_restrictions_to_query( query = query,
                                                                attribute_restrictions_dict = dictRelationAttributeRestriction, 
                                                                column_name_to_restrict="externalEntityRelationID" )
            

            query = self._apply_restrictions_to_query( query = query,
                                                       unification_protocol_name = unification_protocol_name,
                                                       attribute_restrictions= attribute_restrictions,
                                                       column_name_to_restrict="userEntityID" )
            
            query = self._apply_negative_restrictions_to_query( query = query,
                                                                unification_protocol_name = unification_protocol_name,
                                                                negative_attribute_restrictions = negative_attribute_restrictions,
                                                                column_name_to_restrict="userEntityID" )

            interacting_uE.extend([ (x,x,y,z,t) for x,y,z,t in self.db.select_db_content( query, answer_mode="raw" )])

        return interacting_uE



    def get_user_entity_relations(self, unification_protocol_name, userEntityID_list, attribute_restrictions = [], negative_attribute_restrictions = [], listRelationType=[], dictRelationAttributeRestriction={}, use_self_relations=True, limit_to_userEntityID_list=False, use_nested_relations=True):
	#return self.NEWget_user_entity_relations(unification_protocol_name, userEntityID_list, attribute_restrictions, negative_attribute_restrictions, listRelationType, dictRelationAttributeRestriction, use_self_relations, limit_to_userEntityID_list, use_nested_relations)
	return self.OLDget_user_entity_relations(unification_protocol_name, userEntityID_list, attribute_restrictions, negative_attribute_restrictions, listRelationType, dictRelationAttributeRestriction, use_self_relations, limit_to_userEntityID_list, use_nested_relations)

    

    # TO CHECK NEGATIVE RESTRICTIONS
    #def get_expanded_entity_relations(self, unification_protocol_name, userEntityID_list, expansionAttributesList=[], listRelationType=[], use_self_relations=True, limit_to_userEntityID_list=False, expansionLevel=2, attribute_restrictions=[], negative_attribute_restrictions=[]):
    def get_expanded_entity_relations(self, unification_protocol_name, userEntityID_list, expansionAttributesList=[], listRelationType=[], use_self_relations=True, limit_to_userEntityID_list=False, expansionLevel=2, attribute_restrictions=[], negative_attribute_restrictions=[], dictRelationAttributeRestriction = {}):
        """
        unification_protocol_name: name of the unification protocol to be used 
        
        userEntityID_list: user entity ids to be searched for relations based on the shared attributes
        
        expansionAttributesList: list of lists containing (attribute, list_of_parameters_and_values) where list_of_parameter_and_values contains tuples like (parameter, value). Each tuple in the inner list is an expansion criterion (anded with the other expansion criteria in the inner list whereas all inner lists are ored in the query). List contains all the attributes that must be used together for the expansion.
        
        listRelationType: type of relations between attribute sharing entry and its partners to be inferred
        
        use_self_relations: including self relations or not
        
        limit_to_userEntityID_list: use given user entity ids only in inference  # TODO!!! NOT USED NOW.
        
        expansionLevel: number of relations (edges) to look further during inference based on shared attributes
        
        attribute_restrictions: restrictions to be applied on the attributes  # TODO!!!! NOT USED NOW.
        """


        if expansionLevel != 1 and expansionLevel != 2:
            raise ValueError("You must specify level to 1 or 2")


        # First, it is necessary to obtain u2 user entities (those that share attributes with u1) WITHOUT USING RESTRICTIONS
        u1_u2_set = set( self.get_user_entity_relations_by_sharing_attributes( unification_protocol_name = unification_protocol_name,
                                                                               userEntityID_list = userEntityID_list,
                                                                               listAttributes = expansionAttributesList ) )

        # Then, it is necessary to obtain u3 user entities (those having relations with u2). u3 ARE PREDICTIONS!        
        u2_list = [ x[1] for x in u1_u2_set ]


        u2_to_u1 = dict([ (x[1],(x[0], x[2])) for x in u1_u2_set ])

        temporary_predictions_uEr = self.get_user_entity_relations( unification_protocol_name = unification_protocol_name,
                                                                    userEntityID_list = u2_list,
                                                                    listRelationType = listRelationType,
                                                                    dictRelationAttributeRestriction = dictRelationAttributeRestriction,   # Added to add restrictions
                                                                    use_self_relations = True,
                                                                    limit_to_userEntityID_list = False )


        #(80171L, 8118L, 1182323L, 'interaction', 'protein')
        predicted_uEr = []

        # for level 2 predictions 
        u3_externalEntityRelationID = {}
        u3_u1_predictions = {} # Dictionary with u3ids as keys and a list as u1ids values //used in next step
        u3_relation_type = {} # contains type of relation_type between u2-u3
        u1_u2 = {}

	#print temporary_predictions_uEr
        for (uE2, uE3, eErID, eErType, eE3Type) in temporary_predictions_uEr:

            uE1, type = u2_to_u1[uE2]

            #if eE3Type == "self_type": 
            #    eE3Type = type
            
            if use_self_relations is False and uE1==uE3:
                pass
            else:
                predicted_uEr.append( ( uE1,
                                        uE3,
                                        None,
                                        #BianaObjects.UserEntityExpandedRelation( userEntity1Source = uE2,
                                        #                                         userEntity2Source = uE3,
                                        #                                         attribute_identifier = None,   #Not stored
                                        #                                         attribute_values = None,       #Not stored
                                        #                                         externalEntityRelationID = eErID),
                                        "inferred_" + eErType,
                                        eE3Type) )
                
            u3_externalEntityRelationID.setdefault(uE3,set()).add(eErID)
            u3_relation_type[uE3] = eErType
            u3_u1_predictions.setdefault(uE3,[]).append(uE1)
            u1_u2.setdefault(uE1,set()).add(uE2)

        del u2_to_u1


        if expansionLevel==2 and len(temporary_predictions_uEr)>0:

            del temporary_predictions_uEr

            u3set = set(map(str,u3_externalEntityRelationID.keys()))

            done = 0

            itime = time.time()

            for current_uE3 in u3set:
                done += 1

                temp_u3_u4_sharing_attributes = self.get_user_entity_relations_by_sharing_attributes( unification_protocol_name = unification_protocol_name,
                                                                                                      userEntityID_list = [current_uE3],
                                                                                                      listAttributes = expansionAttributesList
                                                                                                      )

                for current_uE3, current_uE4, current_uE4_type in temp_u3_u4_sharing_attributes:
                    for current_eErID in u3_externalEntityRelationID[current_uE3]:                 # TODO: Problem, we don't store the original relation....
                        for current_uE1 in u3_u1_predictions[current_uE3]:
                            if use_self_relations is False and current_uE1==current_uE4:
                                pass
                            else:
                                predicted_uEr.append((current_uE1,
                                                      current_uE4,
                                                      None,
                                                      #BianaObjects.UserEntityExpandedRelation( userEntity1Source = u1_u2[current_uE1],
                                                      #                                         userEntity2Source = current_uE3,
                                                      #                                         attribute_identifier = None,
                                                      #                                         attribute_values = None,
                                                      #                                         externalEntityRelationID = current_eErID ),
                                                      "inferred_" + u3_relation_type[current_uE3],
                                                      current_uE4_type))

                #if done%100 == 0:
                #    print "%s done in %s seconds [%s]" %(done, time.time()-itime,len(predicted_uEr))



            # Filter the predicted uEr with the restrictions, negative_restrictions and limit_to_user_entity_ids

            #print "Filternig the predictions"
            predicted_uE_ids_set = set( [int(x[1]) for x in predicted_uEr] )

            if limit_to_userEntityID_list is True:
                predicted_uE_ids_set = set(map(int,userEntityID_list)).intersection(predicted_uE_ids_set)

            if len(attribute_restrictions)>0:
                all_possible = self.db.select_db_content( self._get_attribute_restrictions_query( unification_protocol_name = unification_protocol_name,
                                                                                                  negative_attribute_restrictions = attribute_restrictions ),
                                                          answer_mode = "list" )
                predicted_uE_ids_set = predicted_uE_ids_set.intersection(set(map(int,all_possible)))

            if len(negative_attribute_restrictions)>0:
                excluded_list = self.db.select_db_content( self._get_attribute_restrictions_query( unification_protocol_name = unification_protocol_name,
                                                                                                   negative_attribute_restrictions = negative_attribute_restrictions ),
                                                           answer_mode = "list" )

                predicted_uE_ids_set = predicted_uE_ids_set.difference(set(map(int,excluded_list)))

        filtered_predicted_uEr = []
        for current_prediction in predicted_uEr:
            if int(current_prediction[1]) in predicted_uE_ids_set:
                filtered_predicted_uEr.append(current_prediction)

	#print filtered_predicted_uEr
        return filtered_predicted_uEr


    ####################################################################################
    #                    EXTERNAL ENTITIES DELETION METHODS                            #
    ####################################################################################

    # At this moment, BIANA does not allow to delete external entities


    ####################################################################################
    #                                   TESTING                                        #
    ####################################################################################



if __name__=="__main__":

    db_access = BianaDBaccess(dbname="biana_v5", dbhost="localhost", dbuser="root", dbsocket="/home/jgarcia/local/mysql/var/mysql.sock")
    print db_access
    db_access.create_database()
    sys.exit()
