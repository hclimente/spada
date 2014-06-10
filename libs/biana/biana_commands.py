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

import traceback
import sys
from OutBianaInterface import OutBianaInterface

# Change prompt
sys.ps1 = "BIANA> "

available_sessions = {}    # Dictionary to store all available sessions

#####################################
## DATABASE ADMINISTRATION METHODS ##
#####################################

class administration:

    def create_biana_database(dbname,dbhost,dbuser,dbpassword,description,dbport=None):
        """
        Creates the BIANA database in the mysql server

        "dbname": Desired database name. It must not exist a database with the same name in the same database server (required)

        "dbhost" is the machine with the mysql server that holds the biana database (required)

        "dbuser" is the mysql user (not required in most systems)

        "dbpassword" is the mysql password (not required in most systems)

        "dbport" is the mysql port (not required in most systems)

        In order to make parsers faster, it creates tables without indices and addicional checkings. Biana automatically enable indices when necessary. If user wants to enable indices manually, it is necessary to use the method "enable_indices"
        """

        OutBianaInterface.send_process_message("Creating new biana database...")

        import BianaDB

        try:
            dbaccess = BianaDB.BianaDBaccess( dbhost = dbhost,
                                              dbuser = dbuser,
                                              dbpassword = dbpassword,
                                              dbport = dbport )
            databases_list = dbaccess.db.select_db_content( sql_query = "SHOW DATABASES", answer_mode = "list" )

            if dbname in databases_list:
                OutBianaInterface.send_error_notification("ERROR: Duplicated database name %s at %s" %(dbname, dbhost), "Database %s exists at %s. If it is a biana database, import it by using the \"Add existing database\" option" %(dbname, dbhost) )

            else:
                dbaccess.db.insert_db_content( sql_query = "CREATE DATABASE IF NOT EXISTS %s DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci" %dbname )
                dbaccess.db.insert_db_content( sql_query = "USE %s" %dbname )
                dbaccess.create_database(description=description,ignore_primary_keys=True,dbname=dbname)
                dbaccess.close()
                administration.check_database(dbname,dbhost,dbuser,dbpassword,dbport)
                OutBianaInterface.send_info_message("BIANA Database correctly created. You can start now populating it using \"Parse external databases\" option")

        except:
            OutBianaInterface.send_error_notification("Error while creating database %s at %s" %(dbname, dbhost), traceback.format_exc())

        OutBianaInterface.send_end_process_message()

    create_biana_database = staticmethod(create_biana_database)


    def delete_biana_database(dbname,dbhost,dbuser,dbpassword,dbport=None):
        """
        Deletes the specified BIANA database.

        ALERT!!! This action is not reversible! All data in BIANA database will be deleted permanently!

        "dbname": BIANA database name

        "dbhost" is the machine with the mysql server that holds the biana database (required)

        "dbuser" is the mysql user (not required in most systems)

        "dbpassword" is the mysql password (not required in most systems)

        "dbport" is the mysql port (not required in most systems)
        """

        import BianaDB
        try:
            dbaccess = BianaDB.BianaDBaccess( dbhost = dbhost,dbuser = dbuser, dbpassword = dbpassword, dbport = dbport )
            dbaccess.db.insert_db_content( sql_query = "DROP DATABASE IF EXISTS %s" %dbname )
            dbaccess.close()
            OutBianaInterface.send_data("<biana_database_deleted dbname=\"%s\" dbhost=\"%s\"/>" %(dbname,dbhost))
        except:
            OutBianaInterface.send_error_notification("Error while deleting database %s at %s" %(dbname,dbhost),"ERROR MESSAGE") # TO CHECK ERROR MESSAGE
        return

    delete_biana_database = staticmethod(delete_biana_database)


    def reset_biana_database(dbname,dbhost,dbuser,dbpassword,dbport=None):
        """
        Deletes all the tables in the specified BIANA database without deleting the database itself.

        ALERT!!! This action is not reversible! All data in BIANA database will be deleted permanently!

        "dbname": Biana database name

        "dbhost" is the machine with the mysql server that holds the biana database (required)

        "dbuser" is the mysql user (not required in most systems)

        "dbpassword" is the mysql password (not required in most systems)

        "dbport" is the mysql port (not required in most systems)
        """

        import BianaDB
        try:
            dbaccess = BianaDB.BianaDBaccess( dbname = dbname, dbhost = dbhost,dbuser = dbuser, dbpassword = dbpassword, dbport = dbport )
            tables = dbaccess.db._get_table_names()
            dbaccess.db.insert_db_content( sql_query = dbaccess.db._get_drop_sql_query(table_list = tables) )
            dbaccess.close()
        except:
            OutBianaInterface.send_error_notification("Error while reseting database %s at %s" %(dbname,dbhost),"ERROR MESSAGE") # TO CHECK ERROR MESSAGE
        return

    reset_biana_database = staticmethod(reset_biana_database)


    def check_database(dbname,dbhost,dbuser,dbpassword,dbport=None):
        """
        Examines the current database and gets its basic information (external databases, unification protocols...)
        """

        print "Checking BIANA database %s at %s..." %(dbname, dbhost)

        OutBianaInterface.send_process_message("Checking BIANA database %s at %s..." %(dbname, dbhost) )

        import BianaDB
        try:
            dbaccess = BianaDB.BianaDBaccess( dbname = dbname,
                                              dbhost = dbhost,
                                              dbuser = dbuser,
                                              dbpassword = dbpassword,
                                              dbport = dbport )


            if( OutBianaInterface.outmethod is not None ):
                tstr = [ "<unification_protocol description=\"%s\"/>" %x for x in dbaccess.get_available_unification_protocols_list() ] 

                source_databases = dbaccess.get_external_database_list()

                eDstr = []
                for current_eD in source_databases:
                    attrStrList = [ "<eDB_external_entity_attribute name=\"%s\"/>" %x for x in current_eD.get_valid_external_entity_attribute_type() ]
                    attrStrList.extend( ["<eDB_external_entity_relation_attribute name=\"%s\"/>" %x for x in current_eD.get_valid_external_entity_relation_attribute_type() ] )
                    attrStrList.extend( ["<eDB_external_entity_type name=\"%s\"/>" %x for x in current_eD.get_valid_external_entity_type() ] )
                    attrStrList.extend( ["<eDB_external_entity_relation_type name=\"%s\"/>" %x for x in current_eD.get_valid_external_entity_relation_type() ] )
                    eDstr.append("<external_database name=\"%s\" version=\"%s\" description=\"%s\" id=\"%s\">%s</external_database>" %(current_eD.get_name(), current_eD.get_version(), current_eD.get_description(), current_eD.get_id(), "".join(attrStrList)))

                ontologies = dbaccess.get_available_ontology_names()
                ontList = ["<ontology_attribute attribute=\"%s\" ontology_name=\"%s\"/>" %(ontologies[x],x) for x in ontologies ]

                OutBianaInterface.send_data("<db_info dbname=\"%s\" dbhost=\"%s\" dbuser=\"%s\" dbpass=\"%s\">%s%s%s</db_info>" %(dbname,dbhost,dbuser,dbpassword,"\n".join(tstr),"\n".join(eDstr),"\n".join(ontList)))

            dbaccess.close()
        except:
            traceback.print_exc()
            OutBianaInterface.send_data("<not_available_database dbname=\"%s\" dbhost=\"%s\"/>" %(dbname,dbhost))
            
        OutBianaInterface.send_end_process_message()


    
    check_database = staticmethod(check_database)


    def create_unification_protocol(unification_protocol_name, list_unification_atom_elements, dbname,dbhost,dbuser,dbpassword,dbport=None):
        """
        Creates a new unification protocol
        
        "unification_protocol_name" is the name of the unification protocol (it can be a string). It cannot contain blank spaces. It must be unique (required)
        
        "listUnifiationAtomElements" is a list of tuples. Each tuple consists on two elements: the list of database to cross and the list of attributes as for example [([1,2],["uniprotAccession"]),([1],["geneSymbol"])]
        
        "dbname": Biana database name
        
        "dbhost" is the machine with the mysql server that holds the biana database (required)
        
        "dbuser" is the mysql user (not required in most systems)

        "dbpassword" is the mysql password (not required in most systems)

        "dbport" is the mysql port (not required in most systems)
        """

        import BianaDB, BianaObjects

        OutBianaInterface.send_process_message("Creating unification protocol. This process can take long time.")

        try:

            # create a temporal session to make all database checkings
            
            import BianaObjects.BianaSessionManager as BianaSessionManager

            import tempfile

            tfile = tempfile.TemporaryFile('w', bufsize=10000)
            temp_session = BianaSessionManager.BianaSessionManager( pSessionID = "temp_session",
                                                                    unification_protocol_name = "No unification",
                                                                    dbname = dbname,
                                                                    dbhost = dbhost,
                                                                    dbuser = dbuser,
                                                                    dbport = dbport,
                                                                    dbpassword = dbpassword,
                                                                    out_method = tfile.write )

            tfile.close()

            temp_session.close()
            del temp_session

            # Commented because it generated an error in the GUI, as we only allow a session. If a session was started when creating a unification protocol, it generated some conflicts
            #create_new_session("temp_session", dbname,dbhost,dbuser,dbpassword,unification_protocol="No unification",dbport=None)

            dbaccess = BianaDB.BianaDBaccess( dbname = dbname,
                                              dbhost = dbhost,
                                              dbuser = dbuser,
                                              dbpassword = dbpassword,
                                              dbport = dbport )

            uProtocol = BianaObjects.UnificationProtocol(unification_protocol_name, "x")  #in X should go the biana database version

            for db_ids, attr_list in list_unification_atom_elements:
                # Add all databases, for adding those databases that are not unified with anything else
                for current_db in db_ids:
                    uProtocol.add_database(current_db)

                if len(attr_list)>0:
                    for current_db1_pos in xrange(len(db_ids)):
                        for current_db2_pos in xrange(current_db1_pos+1): 
                        #for current_db2_pos in xrange(current_db1_pos): # before was not unifying the database with itself (emre)
			    uProtocol.add_unification_atom_elements( BianaObjects.UnificationAtomElement(externalDatabaseID_A=db_ids[current_db1_pos],
                                                                                                         externalDatabaseID_B=db_ids[current_db2_pos],
                                                                                                         externalAttribute=attr_list) )        
            dbaccess.create_new_user_entities(uProtocol)

            OutBianaInterface.send_info_message("New unification protocol successfully created. You can start a working session with it by using \"Create session\" option and selecting database and unification protocol.")

            #OutBianaInterface.send_data(uProtocol.get_xml()) # TO CHECK WHY IS IT COMMENTED

            dbaccess.close()
        except:
            OutBianaInterface.send_error_notification("Error while creating unification protocol",traceback.format_exc())

        OutBianaInterface.send_end_process_message()

        return

    create_unification_protocol = staticmethod(create_unification_protocol)





    def delete_unification_protocol(unification_protocol_name, dbname,dbhost,dbuser,dbpassword,dbport=None):
        """
        Creates a new unification protocol

        "unification_protocol_name" is the name of the unification protocol (it can be a string). It cannot contain blank spaces. It must be unique (required)

        "dbname": Biana database name

        "dbhost" is the machine with the mysql server that holds the biana database (required)

        "dbuser" is the mysql user (not required in most systems)

        "dbpassword" is the mysql password (not required in most systems)

        "dbport" is the mysql port (not required in most systems)
        """

        import BianaDB

        OutBianaInterface.send_process_message("Deleting unification protocol")

        try:

            dbaccess = BianaDB.BianaDBaccess( dbname = dbname,
                                              dbhost = dbhost,
                                              dbuser = dbuser,
                                              dbpassword = dbpassword,
                                              dbport = dbport )

            dbaccess.drop_unification_protocol( unification_protocol_name = unification_protocol_name )

            dbaccess.close()
            
        except:
            OutBianaInterface.send_error_notification("Error while deleting unification protocol",traceback.format_exc())

        OutBianaInterface.send_end_process_message()

        return

    delete_unification_protocol = staticmethod(delete_unification_protocol)



    def get_unification_protocol_atoms(unification_protocol_name, dbname,dbhost,dbuser,dbpassword,dbport=None):
        """
        Returns unification protocol atom information of given unification protocol
        """
        import BianaDB, BianaObjects.output_utilities

        try:
            dbaccess = BianaDB.BianaDBaccess( dbname = dbname,
                                              dbhost = dbhost,
                                              dbuser = dbuser,
                                              dbpassword = dbpassword,
                                              dbport = dbport )

            unification_protocol_atoms = dbaccess.get_unification_protocol_atoms( unification_protocol_name = unification_protocol_name )
            OutBianaInterface.send_data( BianaObjects.output_utilities.get_html_table(columns=["External DB 1", "External DB 2", "Crossed Attributes"], values = [(dbaccess.get_external_database(atom.get_externalDatabaseID_A()), dbaccess.get_external_database(atom.get_externalDatabaseID_B()), atom.get_external_attribute_list()) for atom in unification_protocol_atoms ], title="Unification Protocol Details" ) )

            dbaccess.close()

        except:
            OutBianaInterface.send_error_notification("Error while fetching unification protocol information",traceback.format_exc())

    get_unification_protocol_atoms = staticmethod(get_unification_protocol_atoms)




    def get_available_parsers():
        """
        Returns the description of available parsers
        """
    
        import BianaParser
        xml_list = ["<available_parsers>"]
        xml_list.extend(["<parser name=\"%s\" description=\"%s\" external_entity_definition=\"%s\"/>" %(x[0],x[1],x[2]) for x in BianaParser.get_available_parsers_info()])
        xml_list.append("</available_parsers>")
        OutBianaInterface.send_data("\n".join(xml_list))

    get_available_parsers = staticmethod(get_available_parsers)




    def get_available_external_entity_types():
        """
        Returns all the possible external entity types
        """

        #from biana.BianaObjects.ExternalEntity import ExternalEntity
        #valid_eE_types_list = ExternalEntity.get_valid_external_entity_types()

        import biana.biana_globals as BIANA_GLOBALS
        valid_eE_types_list = list(BIANA_GLOBALS.EXTERNAL_ENTITY_TYPES)
        valid_eE_types_list.sort()
        OutBianaInterface.send_data( "".join( ["<external_entity_type name=\"%s\" />" %x for x in valid_eE_types_list] ) )

    get_available_external_entity_types = staticmethod(get_available_external_entity_types)




    def get_available_external_entity_relation_types():
        """
        Returns all the possible external entity relation types
        """
        import biana.biana_globals as BIANA_GLOBALS
        valid_eEr_types_list = list(BIANA_GLOBALS.EXTERNAL_ENTITY_RELATION_TYPES)
        valid_eEr_types_list.sort()
        OutBianaInterface.send_data( "".join( ["<external_entity_relation_type name=\"%s\" />" %x for x in valid_eEr_types_list] ) )

    get_available_external_entity_relation_types = staticmethod(get_available_external_entity_relation_types)
    

    def get_available_external_entity_attributes():
        """
        Returns all the possible external entity attributes
        """

        import biana.biana_globals as BIANA_GLOBALS
        OutBianaInterface.send_data( "".join( ["<external_entity_attribute name=\"%s\" />" % current_attribute[0] for current_attribute in BIANA_GLOBALS.EXTERNAL_ENTITY_IDENTIFIER_ATTRIBUTES ] ) )
        OutBianaInterface.send_data( "".join( ["<external_entity_attribute name=\"%s\" />" % current_attribute[0] for current_attribute in BIANA_GLOBALS.EXTERNAL_ENTITY_VERSIONABLE_IDENTIFIER_ATTRIBUTE_TYPES]))


    get_available_external_entity_attributes = staticmethod(get_available_external_entity_attributes)







############################
## BIANA SESSION METHODS ###
############################

def create_new_session(sessionID, dbname,dbhost,dbuser,dbpassword,unification_protocol,dbport=None):
    """
    Creates a new biana session with the specified ID. A session must be started in a populated biana database using a specified unification protocol.

    "sessionID" is the identifier for the session. It must be unique! (required)

    "unification_protocol"
    """

    import BianaObjects.BianaSessionManager as BianaSessionManager
    
    if sessionID in available_sessions:
        OutBianaInterface.send_error_notification("Trying to create two sessions with the same ID","Error in session creation")
        return

    try:
        available_sessions[sessionID] = BianaSessionManager.BianaSessionManager(pSessionID = sessionID,
                                                                                unification_protocol_name = unification_protocol,
                                                                                dbname = dbname,
                                                                                dbhost = dbhost,
                                                                                dbuser = dbuser,
                                                                                dbport = dbport,
                                                                                dbpassword = dbpassword,
                                                                                out_method = OutBianaInterface.send_data)
        return available_sessions[sessionID]
    except:
        OutBianaInterface.send_error_notification("Error in session creation.", traceback.format_exc())


def save_session(sessionID, file_name):
    """
    Saves the Session in the specified file

    "sessionID" is the identifier of the session

    "fileName" is the absolute or relative path where the session must be saved
    """
    import cPickle
    outfile_fd = open(file_name, 'w')
    objBianaSessionManager = available_sessions[sessionID]
    cPickle.dump(objBianaSessionManager, outfile_fd, cPickle.HIGHEST_PROTOCOL)
    outfile_fd.close()
    return


def remove_session(sessionID):
    """
    Removes a session from memory

    "sessionID" is the identifier of the session
    """
    
    available_sessions[sessionID].close()
    del available_sessions[sessionID]
    if( sessionID is not None ):
        OutBianaInterface.send_data("<close_session sessionID=\"%s\"/>" %sessionID)
    return



def load_session(file_name):
    """
    Loads a saved biana session
    """
    import cPickle
    infile_fd = open(file_name)
    session = cPickle.load(infile_fd)
    if available_sessions.has_key(session.sessionID):
        OutBianaInterface.send_error_notification("Load session error","Trying to load an existing session (same ID)")
    else:
        available_sessions[session.sessionID] = session

    infile_fd.close()
    return

def ping():
    OutBianaInterface.send_data("<biana_ping_response />")

def close():
    OutBianaInterface.close()
    sys.exit()
    



