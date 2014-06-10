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

# Class to manage/handle (create, analyze, operate vs..) networks 
# from userEntities and available interaction data in BianaDB (using BianaDBaccess) 

import sys

from biana.BianaDB.BianaDBaccess import BianaDBaccess
import UserEntitySet
import UserEntity

from ExternalEntityAttribute import ExternalEntityAttribute

import output_utilities
from biana.utilities import graph_utilities
from BianaReport import *
import copy
import traceback
from biana.OutBianaInterface import OutBianaInterface


class Enum(object):

    def __init__(self):
        self._enum_type_2_letter = {}
        self._letter_2_enum_type = {}
        self._enum_index = -1
        self._enum_letter = "abcdefghijklmnopqrstuvwxyz1234567890"

    def get_letter(self, enum):
        try:
            return self._enum_type_2_letter[enum.lower()]
        except:
            next = self._get_next()
            self._letter_2_enum_type[next] = enum
            self._enum_type_2_letter[enum.lower()] = next
            return next
    
    def _get_next(self):
        self._enum_index += 1
        return self._enum_letter[self._enum_index]

    def get(self, letter):
        return self._letter_2_enum_type[letter]


class BianaSessionManager(object):
    """
    A class to represent a BianaNetworkManager
    
    A BianaNetworkManager is composed of one or several userEntity(s) and one BianaDBaccess object 
    """
    
    # attribute substitution list for outputting identifiers of nodes in network when none of the external entities in the node has a valid value for given attribute
    substitution_list = [ "uniprotaccession", "geneid" ]

    #uE_types_enum = Enum()
    #eEr_types_enum = Enum()


    def __init__(self, pSessionID, unification_protocol_name, dbname, dbhost, dbuser=None, dbpassword=None, dbport=None, out_method=None):
        """
        Starts a new BIANA Working session
        ------
        pSessionID: identifier for current session object
        unification_protocol_name: name of the unification protocol to be used while retrieving data from database
        dbname: name of the mysql database to be used
        dbhost: host where mysql database resides
        dbuser: user for the given mysql database
        dbpassword: password for the given mysql database and host
        dbport: port for mysql database connection, if None default value is used
        out_method: is where biana session manager has to inform about changes
        """

	self.uE_types_enum = Enum()
	self.eEr_types_enum = Enum()

        self.sessionID = pSessionID
        self.unification_protocol_name = unification_protocol_name
        
        self.dbname = dbname
        self.dbhost=dbhost
        self.dbuser=dbuser
        self.dbpassword=dbpassword
        self.dbport=dbport

	OutBianaInterface.send_process_message("Checking database integrity... BIANA will remove data from unfinalized parsing attempts (if there had been any).")
        self.dbAccess = BianaDBaccess(dbname=dbname,
                                      dbhost=dbhost,
                                      dbuser=dbuser,
                                      dbpassword=dbpassword,
                                      dbport=dbport,
				      check_integrity=True)
	OutBianaInterface.send_end_process_message()
        
        
        self.dictUserEntity = {}                     # Dictionary to store user entity objects
        self.dictUserEntitySet = {}                  # dictionary to store user entity sets objects

        self.dictUserEntityInfo = {}                 # dictionary to store user entity basic info (not objects, as they use a lot of memory)
        self.dictExternalEntityInfo = {}             # dictionary to store external entity basic info (not objects, as they use a lot of memory)    

        self.uE_types_dict = {}
        self.eEr_types_dict = {}

        self.idLastUserEntitySet = 0

        self.selectedUserEntitySetsIds = set()

        self.outmethod = out_method

        self.report = None

        self.outmethod("<new_session id=\"%s\" dbname=\"%s\" dbhost=\"%s\" unification_protocol=\"%s\" description=\"Session description\"/>" %(self.sessionID,dbname,dbhost,self.unification_protocol_name))

        # Checks if database is optimized for running. If not, optimize it
        if self.dbAccess.isOptimizedForRunning()==False:
            OutBianaInterface.send_process_message("Optimizing database... Depending on BIANA database size, this process can take from few seconds to some hours.")
            self.dbAccess.optimize_database_for( mode="running" )
            OutBianaInterface.send_end_process_message()
        return


    # ----
    # methods required for using pickle with piana objects
    # ----
    def __getstate__(self):
        """
        Get state of the session for [c]pickle.dump()
        """
        odict = self.__dict__.copy() # copy the dict since we change it
        del odict['dbAccess']              # remove database entry
        del odict['outmethod']              # remove static outmethod entry
        return odict

    def __setstate__(self, dict):
        """
        Get state of the session for [c]pickle.load()
        ------
        dict: the dictionary of a previously created and saved BianaSessionManager object
        """
        self.__dict__.update(dict) # update attributes
        self.dbAccess = BianaDBaccess(dict["dbname"], dict["dbhost"], dict["dbuser"], dict["dbpassword"], dict["dbport"])
        self.outmethod = OutBianaInterface.send_data

        self.outmethod("<new_session id=\"%s\" dbname=\"%s\" dbhost=\"%s\" unification_protocol=\"%s\" description=\"Session description\"/>" %(self.sessionID,self.dbname,self.dbhost,self.unification_protocol_name))
        
        OutBianaInterface.send_process_message("Sending data...")
        for idUESet, user_entity_set in self.dictUserEntitySet.iteritems():
            self._send_complete_user_entity_set_info(user_entity_set=user_entity_set)
        OutBianaInterface.send_end_process_message()

        # Checks if database is optimized for running. If not, optimize it
        if self.dbAccess.isOptimizedForRunning()==False:
            OutBianaInterface.send_process_message("Optimizing database...")
            self.dbAccess.optimize_database_for( mode="running" )
            OutBianaInterface.send_end_process_message()

        return

    def close(self):
        self.dbAccess.close()

    
    def reconnect(self):
        self.dbAccess.reconnect()


    def _send_complete_user_entity_set_info(self, user_entity_set):
        """
        Sends information of given userEntitySet in xml format through self.outmethod
        ------
        user_entity_set: an instance of user entity set
        """

        self.outmethod(self._get_xml(inner_content="<new_user_entity_set id=\"%s\"/>" %(user_entity_set.id)))
        self.outmethod(self._get_xml(inner_content=self._get_user_entity_set_xml(user_entity_set)))
        if( user_entity_set.isNetworkCreated() == True ):
            self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content="<update_network_depth levels=\"%s\"/>" %user_entity_set.get_level())))


    def get_ontology(self, ontology_name, root_attribute_values=[]):
        """
        Fetchs the ontology structure identifier by ontology_name
        ------
        ontology_name: name of the ontology to be retrieved from database
        root_attribute_values: list of values used in selecting roots of ontology by attribute
        """

        #OutBianaInterface.send_process_message("Getting information from database...")

        #print "loading ontology from database"

        ontology_obj =  self.dbAccess.get_ontology( ontology_name = ontology_name, root_attribute_values = root_attribute_values, load_external_entities = False )  # Changed True to False

        #print "loaded"

        #OutBianaInterface.send_end_process_message()

        return ontology_obj


    
    
    def get_user_entity_set(self, user_entity_set_id):
        """
        returns the user entity set object with the name "user_entity_set_id"
        ------
        user_entity_set_id: identifier of user entity set to be returned
        """

        try:
            return self.dictUserEntitySet[user_entity_set_id]
        except:
            OutBianaInterface.send_error_notification( message = "ERROR", error = "Trying to get an unexisting user entity set: %s" %user_entity_set_id)

        return


    def select_user_entity_set(self, user_entity_set_id):
        """
        Method for selecting a user entity set
        ------
        user_entity_set_id: identifier of user entity set to be selected 

        # Used by BIANAGUI
        """

        if self.dictUserEntitySet.has_key(user_entity_set_id):
            self.selectedUserEntitySetsIds.add(user_entity_set_id)
            self.outmethod("<select_user_entity_set id=\"%s\"/>" %user_entity_set_id)
        else:
            OutBianaInterface.send_error_notification( message = "ERROR", error = "Trying to select an unexisting user entity set: %s" %user_entity_set_id)

        return


    def describe_user_entity_set(self, user_entity_set_id):
        """
        Method to print into the biana interface the description of the User Entity ID identified by "user_entity_set_id"

        "format" can be txt or html
        """
        
        if self.dictUserEntitySet.has_key(user_entity_set_id):
            OutBianaInterface.send_info_message(self.dictUserEntitySet[user_entity_set_id].__str__())
        else:
            OutBianaInterface.send_error_notification( message = "ERROR", error = "Trying to get an unexisting user entity set: %s" %user_entity_set_id)
        
        return


    def get_selected_user_entity_set_ids(self):
        """
        Returns list of identifiers of selected user entity sets

        # BIANAGUI method
        """
        return self.selectedUserEntitySetsIds


    def clear_selected_user_entity_sets(self):
        """
        Clears selection of user entity sets

        # BIANA GUI method
        """
        self.selectedUserEntitySetsIds.clear()
        self.outmethod("<clear_user_entity_set_selection/>")
        return

    def _convert_attribute_list_to_attribute_dictionary(self, identifier_description_list, attribute_name):
        """
        Returns dictionary containing identifier types mapped to identifiers created from given list of identifiers
        ------
        identifier_description_list: list of identifiers (or (identifier, id_type) tuples in case id_type is "embedded")
        id_type: type of the identifiers in the file, if "embedded" then file contains (attribute_name, identifier) tuples instead of just identifiers
        """
        dictTypeToName = {}
        for identifierDescription in identifier_description_list:
            if attribute_name == "embedded":
                identifierType = identifierDescription[0].lower()
                identifierString = identifierDescription[1]
            else:
                identifierType = attribute_name.lower()
                identifierString = identifierDescription
	    # Strip versioning parts of the values if they are versionable identifiers
	    if identifierType in self.dbAccess.get_versionable_external_entity_identifier_attributes():
		index = max(identifierString.find("-"), identifierString.find("."))
		if index > 0:
		    identifierString = identifierString[:index]
            dictTypeToName.setdefault(identifierType, set()).add(identifierString)
        return dictTypeToName

    def _get_next_uEs_id(self):
        """
        Private method to obtain automatically the next user entity set default id
        """
        self.idLastUserEntitySet += 1
        return "uEs_%s" %self.idLastUserEntitySet


    def duplicate_user_entity_set(self, user_entity_set_id, new_user_entity_set_id):
        """
        Returns an exact copy of user entity set object with given id
        ------
        user_entity_set_id: id of user entity set to be duplicated
        user_entity_set_id_new: id of created copy of user entity set
        """
        original_uEs = self.get_user_entity_set(user_entity_set_id)

        if original_uEs is not None:
            newObj = copy.deepcopy(original_uEs)
            newObj.id = new_user_entity_set_id
            self.dictUserEntitySet[newObj.id] = newObj
            self._send_complete_user_entity_set_info(user_entity_set=newObj)
        
        return newObj


    def create_randomized_user_entity_set(self, user_entity_set_id, new_user_entity_set_id, type_randomization):
        """
        Creates a new user entity set with randomized network from given user entity set
        ------
        user_entity_set_id: id of the user entity set whose copy with random network is going to be created
        new_user_entity_set_id: id for the created copy of user entity set
        type_randomization: randomization type to be used in network randomization, can be one of the following: "random", "preserve_topology", "preserve_topology_and_node_degree", "preserve_degree_distribution", "preserve_degree_distribution_and_node_degree"
        """
        original_uEs = self.get_user_entity_set(user_entity_set_id)

        if original_uEs is not None:
            if original_uEs.isNetworkCreated():
                newObj = copy.deepcopy(original_uEs)
                newObj.id = new_user_entity_set_id
                original_network = original_uEs.getNetwork()
		use_self_relations = original_uEs.getRestrictions(restriction_type="use_self_relations")
                random_network = graph_utilities.randomize_graph(graph = original_network, randomization_type = type_randomization, allow_self_edges = use_self_relations)
                #print type_randomization, original_network.number_of_edges(), random_network.number_of_edges()
                #newObj.setNetwork(graph_utilities.randomize_graph(original_uEs.getNetwork(), type_randomization))
                newObj.setNetwork(random_network)
                self.dictUserEntitySet[newObj.id] = newObj
                self._send_complete_user_entity_set_info(user_entity_set=newObj)
                return newObj
            else:
                OutBianaInterface.send_error_notification( message = "Cannot randomize without a created network", error = "Set is not created as it does not have a network" )
        else:
            OutBianaInterface.send_error_notification( message = "Randomization not done!", error = "Cannot randomize with an unexisting set: %s" %user_entity_set_id)

        return


    def randomize_user_entity_set_network(self, user_entity_set_id, type_randomization):
        """
        Randomizes network of a given user entity set.
        ------
        user_entity_set_id: id of the user entity set whose network will be randomized
        type_randomization: randomization type to be used in network randomization, can be one of the following: "random", "preserve_topology", "preserve_topology_and_node_degree", "preserve_degree_distribution", "preserve_degree_distribution_and_node_degree"
        """

        OutBianaInterface.send_process_message("Creating new user entity set. \nGetting information from database...")

        original_uEs = self.get_user_entity_set(user_entity_set_id)
        if original_uEs is not None:
            if original_uEs.isNetworkCreated():
		use_self_relations = original_uEs.getRestrictions(restriction_type="use_self_relations")
                original_uEs.setNetwork(graph_utilities.randomize_graph(original_uEs.getNetwork(), type_randomization, allow_self_edges = use_self_relations))
                #self._send_complete_user_entity_set_info(user_entity_set=original_uEs)
                #self.outmethod(self._get_xml(inner_content="<remove_user_entity_set id=\"%s\"/>" % user_entity_set_id))
                #self._send_user_entity_set_network_info(original_uEs)
                self.outmethod(self._get_xml(inner_content="<remove_user_entity_set id=\"%s\"/>" % user_entity_set_id))
                self._send_complete_user_entity_set_info(user_entity_set=original_uEs)

        OutBianaInterface.send_end_process_message()

        return


    def create_new_user_entity_set(self, identifier_description_list, id_type="embedded", new_user_entity_set_id=None, attribute_restriction_list=[], negative_attribute_restriction_list=[], external_database_restrictions=[]): 
        """
        create userEntity objects from given externalEntity ids (or get if already existing) then add them to a userEntitySet
        ------
        identifier_description_list: list with external identifiers of the nodes in the set provided by user 
        id_type: either list of attrubutes defined in EXTERNAL_ENTITY_TABLES or "embedded" meaning identifier_description_list is a list of (id_type, identifier) tuple in the form [(type1, id1), (type2, id2), ..., (typeN, idN)] 
        new_user_entity_set_id: identifier of the set provided by the user
        attribute_restriction_list: list of tuples provided by the user containing external entity attribute restrictons to be applied. They will be used always in the set (network creation, etc)
        negative_attribute_restriction_list: list of tuples provided by the user containing attributes that may neve appear in the user entity set. They will be used alwyas in the set.
        """
        
        OutBianaInterface.send_process_message("Creating new user entity set. \nGetting information from database...")

        try:

            # Check for missing parameters and set it to default
            if new_user_entity_set_id is None:
                new_user_entity_set_id = self._get_next_uEs_id()

            user_entity_set = UserEntitySet.UserEntitySet(new_user_entity_set_id )


            # Javi added: transform restrictions (for transferred attributes)
            if id_type=="embedded":
                identifier_description_list = self.dbAccess.transform_expanded_attribute_restrictions(identifier_description_list)
            else:
                identifier_description_list = self.dbAccess.transform_expanded_attribute_restrictions([ (id_type,x) for x in identifier_description_list ])

            attribute_restriction_list = self.dbAccess.transform_expanded_attribute_restrictions(attribute_restriction_list)
            negative_attribute_restriction_list = self.dbAccess.transform_expanded_attribute_restrictions(negative_attribute_restriction_list)


            ## prepare id_type classified sets to be used in data fetching step
            # Gets the input to build the user entity set (root or seed user entities)
            #dictTypeToName = self._convert_attribute_list_to_attribute_dictionary(identifier_description_list, id_type)
            dictTypeToName = self._convert_attribute_list_to_attribute_dictionary(identifier_description_list, "embedded")

            user_entity_set.userProvidedExternalIdsDict = dictTypeToName

            # Fill userProvidedExtIdsLowerDict:
            for k in user_entity_set.userProvidedExternalIdsDict:
                user_entity_set.userProvidedExtIdsLowerDict[k] = {}

                for extId in user_entity_set.userProvidedExternalIdsDict[k]:
                    user_entity_set.userProvidedExtIdsLowerDict[k][str(extId).lower()] = extId


            # selectedExtIdsLowerDict is filled by method UserEntitySet._update_selected_external_ids_dict()

            ## for each different type of id_types fetch userEntity ids associated with given externalEntity ids
            for identifierType, setIdentifierName in dictTypeToName.iteritems():

                listIdUserEntity = []

                # javi added:
                if( identifierType.lower()=="userentityid" ):
		    uEId_to_type = self.dbAccess.get_user_entity_type(unification_protocol_name = self.unification_protocol_name, user_entity_ids = setIdentifierName)
		    listIdUserEntity = uEId_to_type.items() 
                else:
                    field_values = []
                    for identifierName in setIdentifierName:
                        field_values.append(("value",identifierName))
                    listIdUserEntity = self.dbAccess.get_list_user_entities_IDs_by_attribute( unification_protocol_name = self.unification_protocol_name, 
			    attribute_identifier = identifierType,
			    field_values = field_values, 
			    attribute_restrictions = attribute_restriction_list, 
			    negative_attribute_restrictions = negative_attribute_restriction_list, 
			    include_type = True )


                ## create new userEntity objects if not created (by another external entity that belongs to that userEntity) before 
                for (user_entity_id, type) in listIdUserEntity:
                    user_entity_set.addUserEntityId(idUserEntity = user_entity_id, level = 0)
                    self.uE_types_dict[user_entity_id] = self.uE_types_enum.get_letter(type)
                    self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(user_entity_id, type))

                user_entity_set.addRestriction("negative_attribute_restrictions", negative_attribute_restriction_list)
                user_entity_set.addRestriction("attribute_restrictions", attribute_restriction_list)
                user_entity_set.addRestriction("externalDatabase_restrictions", external_database_restrictions)
        
        except:
            OutBianaInterface.send_error_notification( message = "New set not created. BIANA ERROR:", error = traceback.format_exc() )
            OutBianaInterface.send_end_process_message()
            return

        OutBianaInterface.send_end_process_message()
        
        if user_entity_set.getSize()==0:
            OutBianaInterface.send_error_notification(message="New set not created", error = "Any element has been found with given identifiers")
            return

        self.dictUserEntitySet[user_entity_set.id] = user_entity_set

        OutBianaInterface.send_process_message("Sending data...")
        self._send_complete_user_entity_set_info( user_entity_set = user_entity_set )
        OutBianaInterface.send_end_process_message()

        return user_entity_set


    def _get_user_entity_set_xml(self, user_entity_set_obj, level=None):
        
        user_entity_ids_list = user_entity_set_obj.get_user_entity_ids(level=level) 
        relations_str_list = []
        for (id1, id2) in user_entity_set_obj.getRelations():
            eEr_ids_list = user_entity_set_obj.get_external_entity_relation_ids(id1,id2)
            types = set()
            for current_eEr_id in eEr_ids_list:
                types.add(self.eEr_types_dict[current_eEr_id])
            for current_type in types:
                relations_str_list.append("<user_entity_relation node1=\"%s\" node2=\"%s\" type=\"%s\" relation_id=\"%s\"/>" %(id1, id2, self.eEr_types_enum.get(current_type), "0" ))  #TO CHECK IF IT IS NECESSARY TO PRINT THE RELATION ID...

        default_attributes_dict = dict([ (x,x) for x in user_entity_set_obj.get_groups_ids() ])
        default_attributes_dict.update(self.dbAccess.get_default_external_entity_ids( externalEntityIDsList=user_entity_set_obj.get_groups_ids() ))

        uE_tags_xml = [ "<new_tag tag=\"%s\"/>" %tag for tag in user_entity_set_obj.get_all_user_entity_tags() ]
        uEr_tags_xml = [ "<new_relation_tag tag=\"%s\"/>" %tag for tag in user_entity_set_obj.get_all_user_entity_relation_tags() ]
        

        return "%s%s%s%s%s%s%s" %(user_entity_set_obj._get_xml_header(),
                                  "".join(["<user_entity id=\"%s\" type=\"%s\" />" %(id,self.uE_types_enum.get(self.uE_types_dict[id])) for id in user_entity_ids_list]),
                                  "".join(relations_str_list),
                                  user_entity_set_obj._get_groups_xml(only_not_printed=False, group_identifiers = default_attributes_dict),
                                  uE_tags_xml,
                                  uEr_tags_xml,
                                  user_entity_set_obj._get_xml_foot() )

    



    def remove_user_entity_set(self, user_entity_set_id):
        """
        Remove user entity set from the session
        ------
        user_entity_set_id: identifier of user entity set to be removed
        """
        if self.dictUserEntitySet.has_key(user_entity_set_id):
            del self.dictUserEntitySet[user_entity_set_id]
            self.outmethod(self._get_xml(inner_content="<remove_user_entity_set id=\"%s\"/>" % user_entity_set_id))
        return

    def get_user_entity(self, user_entity_id):
        """
        Fetch all externalEntity ids associated with this user entity or retrieve it from dictionary if exists
        ------
        user_entity_id: identifier of user entity 
        """
        if not self.dictUserEntity.has_key(user_entity_id):
            listIdExternalEntity = self.dbAccess._get_list_eE_for_uE(self.unification_protocol_name, user_entity_id)
            objUserEntity = UserEntity.UserEntity(user_entity_id, listIdExternalEntity)
            self.dictUserEntity[user_entity_id] = objUserEntity
        else:
            objUserEntity = self.dictUserEntity[user_entity_id]
        return objUserEntity

    def expand_user_entity_set(self, user_entity_set_id, is_last_level=False):
        """
        Fetchs interactions of userEntities in the last level of the userEntitySet
        ------
        user_entity_set_id: identifier of the user entity set
        """
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        
        #if( len(user_entity_set.getRestrictions("relation_type_restrictions"))==0 ):
        if user_entity_set.isNetworkCreated() == False:
            OutBianaInterface.send_error_notification( message = "Cannot expand network!", error = "Before expanding a network, it is mandatory to initialize a network with the create network command")
            return
            
        OutBianaInterface.send_process_message("Getting information from database...")


        self.outmethod(self._get_xml_header())
        self.outmethod(user_entity_set._get_xml_header())

        try:

            temp_inserted = 0

            userEntityIds = user_entity_set.get_user_entity_ids(level="last")

            # Only creates the network of relations if types are specified
            if( len(user_entity_set.getRestrictions("relation_type_restrictions"))>0 ):
                listTupleIdUserEntity = self.dbAccess.get_user_entity_relations(unification_protocol_name = self.unification_protocol_name,
                                                                                userEntityID_list = userEntityIds,
                                                                                attribute_restrictions = user_entity_set.getRestrictions("attribute_restrictions"),
                                                                                negative_attribute_restrictions = user_entity_set.getRestrictions("negative_attribute_restrictions"),
                                                                                listRelationType = user_entity_set.getRestrictions("relation_type_restrictions"),
                                                                                dictRelationAttributeRestriction = user_entity_set.getRestrictions("relation_attribute_restrictions"),
                                                                                use_self_relations = user_entity_set.getRestrictions("use_self_relations"),
                                                                                limit_to_userEntityID_list = is_last_level) # ramon removes self. that was placed # before each user_entity_set

                OutBianaInterface.send_process_message("Processing information...")
             
                for (idUserEntity1, idUserEntity2, externalEntityRelationID, relation_type, partner_type) in listTupleIdUserEntity:
                    if not user_entity_set.has_user_entity(idUserEntity1):
                        self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity1, partner_type))
                    elif not user_entity_set.has_user_entity(idUserEntity2):
                        self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity2, partner_type))
		    # addUserEntityRelation adds nodes so above should be before a call to it
                    inserted = user_entity_set.addUserEntityRelation(idUserEntity1 = idUserEntity1,
                                                                     idUserEntity2 = idUserEntity2, 
                                                                     externalEntityRelationID = externalEntityRelationID)
                    self.eEr_types_dict[externalEntityRelationID] = self.eEr_types_enum.get_letter(relation_type)
                    #self.uE_types_dict.setdefault(idUserEntity1,self.uE_types_enum.get_letter(partner_type))
                    self.uE_types_dict.setdefault(idUserEntity2,self.uE_types_enum.get_letter(partner_type))

                    if inserted:
                        self.outmethod("<user_entity_relation node1=\"%s\" node2=\"%s\" type=\"%s\" relation_id=\"%s\"/>" %(idUserEntity1, idUserEntity2, relation_type, externalEntityRelationID))
                        temp_inserted += 1

                OutBianaInterface.send_end_process_message()

            # Creates the relation groups if types for groups are specified
            relation_ids_set  = set()
            if( len(user_entity_set.getRestrictions("group_relation_type"))>0 ):
                listTupleIdUserEntity = self.dbAccess.get_user_entity_relations(unification_protocol_name = self.unification_protocol_name,
                                                                                userEntityID_list = userEntityIds,
                                                                                attribute_restrictions = user_entity_set.getRestrictions("attribute_restrictions"),
                                                                                negative_attribute_restrictions = user_entity_set.getRestrictions("negative_attribute_restrictions"),
                                                                                listRelationType = user_entity_set.getRestrictions("group_relation_type"),
                                                                                dictRelationAttributeRestriction = user_entity_set.getRestrictions("relation_attribute_restrictions"),
                                                                                use_self_relations = user_entity_set.getRestrictions("use_self_relations"),
                                                                                limit_to_userEntityID_list = is_last_level) # ramon removes self. that was placed # before each user_entity_set

                OutBianaInterface.send_process_message("Processing information...")

                for (idUserEntity1, idUserEntity2, externalEntityRelationID, relation_type, partner_type) in listTupleIdUserEntity:
                    self.uE_types_dict.setdefault(idUserEntity2, self.uE_types_enum.get_letter(partner_type))
                    relation_ids_set.add(externalEntityRelationID)
                    if not user_entity_set.has_user_entity(idUserEntity1):
                        self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity1, partner_type))
                    elif not user_entity_set.has_user_entity(idUserEntity2):
                        self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity2, partner_type))
                    user_entity_set.addUserEntitiesToGroup(group_id = externalEntityRelationID, 
                                                           userEntityID = idUserEntity2,
                                                           representative_userEntityID = idUserEntity1,
                                                           group_type = relation_type )


                    
                # Get groups hierarchy information
                if len(listTupleIdUserEntity)>0:
		    #eE_dict = self.dbAccess.get_external_entities_dict( attribute_list=["name"], externalEntityIdsList = user_entity_set.get_groups_ids() )
		    default_attributes_dict = dict([ (x,x) for x in user_entity_set.get_groups_ids() ])
		    default_attributes_dict.update(self.dbAccess.get_default_external_entity_ids( externalEntityIDsList=user_entity_set.get_groups_ids() ))
                    user_entity_set.setGroupsHierarchy( list_hierarchy = self.dbAccess.get_relations_hierarchy( externalEntityRelationIDs = relation_ids_set ) )
                    self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content=user_entity_set._get_groups_xml(only_not_printed=True,  group_identifiers = default_attributes_dict ))))     
                    
                OutBianaInterface.send_end_process_message()

            if( len(user_entity_set.getRestrictions( restriction_type = "expansionAttributesList"))>0 ):

                # Get inferred relations
                OutBianaInterface.send_process_message("Getting information from database...")

                listTupleIdUserEntity = self.dbAccess.get_expanded_entity_relations(unification_protocol_name = self.unification_protocol_name,
                                                                                    userEntityID_list = userEntityIds,
                                                                                    listRelationType = user_entity_set.getRestrictions("expansionListRelationType"),
                                                                                    use_self_relations = user_entity_set.getRestrictions("use_self_relations"),
                                                                                    expansionLevel = user_entity_set.getRestrictions("expansionLevel"),
                                                                                    expansionAttributesList = user_entity_set.getRestrictions("expansionAttributesList"),
                                                                                    attribute_restrictions = user_entity_set.getRestrictions("attribute_restrictions"),
                                                                                    negative_attribute_restrictions = user_entity_set.getRestrictions("negative_attribute_restrictions"),
                                                                                    dictRelationAttributeRestriction = user_entity_set.getRestrictions("relation_attribute_restrictions"), # Added to restrict also expansion with respect to original relation
                                                                                    limit_to_userEntityID_list = is_last_level)

                for (idUserEntity1, idUserEntity2, externalEntityRelationID, relation_type, partner_type) in listTupleIdUserEntity:
                    self.uE_types_dict.setdefault(idUserEntity2, self.uE_types_enum.get_letter(partner_type))
                    #self.uE_types_dict.setdefault(idUserEntity1, self.uE_types_enum.get_letter(partner_type))
                    self.eEr_types_dict.setdefault(externalEntityRelationID, self.eEr_types_enum.get_letter(relation_type))
		    if not user_entity_set.has_user_entity(idUserEntity1):
			self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity1, partner_type))
		    elif not user_entity_set.has_user_entity(idUserEntity2):
			self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity2, partner_type))
                    inserted = user_entity_set.addUserEntityRelation(idUserEntity1=idUserEntity1, idUserEntity2=idUserEntity2, externalEntityRelationID=externalEntityRelationID)
                    if inserted:
                        self.outmethod("<user_entity_relation node1=\"%s\" node2=\"%s\" type=\"%s\" relation_id=\"%s\"/>" %(idUserEntity1, idUserEntity2, relation_type, externalEntityRelationID))
                OutBianaInterface.send_end_process_message()

            if( len(user_entity_set.getRestrictions( restriction_type = "attributeNetworkList"))>0 ):

                # Get attribute relations
                OutBianaInterface.send_process_message("Getting attribute relations...")

                listTupleIdUserEntity = self.dbAccess.get_user_entity_relations_by_sharing_attributes( unification_protocol_name = self.unification_protocol_name,
                                                                                                       userEntityID_list = userEntityIds,
                                                                                                       listAttributes = user_entity_set.getRestrictions("attributeNetworkList"),
                                                                                                       limit_to_userEntityID_list = is_last_level,
                                                                                                       attribute_restrictions = user_entity_set.getRestrictions("attribute_restrictions"),
                                                                                                       negative_attribute_restrictions = user_entity_set.getRestrictions("negative_attribute_restrictions") )

                for (idUserEntity1, idUserEntity2, type) in listTupleIdUserEntity:
                    self.uE_types_dict.setdefault(idUserEntity2, self.uE_types_enum.get_letter(type))
                    #self.uE_types_dict.setdefault(idUserEntity1, self.uE_types_enum.get_letter(type))
                    self.eEr_types_dict[0] = self.eEr_types_enum.get_letter("common_attribute")
		    if not user_entity_set.has_user_entity(idUserEntity1):
			self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity1, type))
		    elif not user_entity_set.has_user_entity(idUserEntity2):
			self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity2, type))
                    inserted = user_entity_set.addUserEntityRelation(idUserEntity1=idUserEntity1, idUserEntity2=idUserEntity2, externalEntityRelationID=0)
                    if inserted:
                        self.outmethod("<user_entity_relation node1=\"%s\" node2=\"%s\" type=\"%s\" relation_id=\"%s\"/>" %(idUserEntity1, idUserEntity2, "common_attribute", 0))
                OutBianaInterface.send_end_process_message()

        except:
            OutBianaInterface.send_error_notification( message = "Network not expanded. BIANA ERROR:", error = traceback.format_exc() )

        self.outmethod("<update_network_depth levels=\"%s\"/>" %user_entity_set.get_level())
        self.outmethod(user_entity_set._get_xml_foot())
        self.outmethod(self._get_xml_foot())

        OutBianaInterface.send_end_process_message()


    def _get_xml(self, inner_content):
        """
        Returns given information for current session in XML format.
        ------
        inner_content: Information to be encapsulated with current session.
        """
        return "%s%s%s" %(self._get_xml_header(), inner_content, self._get_xml_foot())

    def _get_xml_header(self):
        return "<session id=\"%s\">" %(self.sessionID)

    def _get_xml_foot(self):
        return "</session>"


    def create_relation_type_set(self, new_user_entity_set_id, relation_type_list, attribute_restriction_list=[], negative_attribute_restriction_list=[], relation_attribute_restriction_list=[], use_self_relations=True):
        """
        Creates a user entity set with all the user entities and relations of relation_type_list
        """
        
        OutBianaInterface.send_process_message("Creating new user entity set. \nGetting information from database...")

        try:

            # Check for missing parameters and set it to default
            if new_user_entity_set_id is None:
                new_user_entity_set_id = self._get_next_uEs_id()

            user_entity_set = UserEntitySet.UserEntitySet(new_user_entity_set_id )

            dictAttributeToValues = self._convert_attribute_list_to_attribute_dictionary(relation_attribute_restriction_list, "embedded")

            user_entity_set.addRestriction("negative_attribute_restrictions", negative_attribute_restriction_list)
            user_entity_set.addRestriction("attribute_restrictions", attribute_restriction_list)
            user_entity_set.addRestriction("relation_type_restrictions", relation_type_list)
            user_entity_set.setRestriction("relation_attribute_restrictions", dictAttributeToValues)
            user_entity_set.setRestriction(restriction_type="use_self_relations",restriction=use_self_relations)
            
            attribute_restriction_list = self.dbAccess.transform_expanded_attribute_restrictions(attribute_restriction_list)
            negative_attribute_restriction_list = self.dbAccess.transform_expanded_attribute_restrictions(negative_attribute_restriction_list)

            listTupleIdUserEntity = self.dbAccess.get_relations( unification_protocol_name= self.unification_protocol_name, 
                                                                 attribute_restrictions = user_entity_set.getRestrictions("attribute_restrictions"), 
                                                                 negative_attribute_restrictions = user_entity_set.getRestrictions("negative_attribute_restrictions"),
                                                                 listRelationType = user_entity_set.getRestrictions("relation_type_restrictions"), 
                                                                 dictRelationAttributeRestriction = user_entity_set.getRestrictions("relation_attribute_restrictions"),
                                                                 use_self_relations = user_entity_set.getRestrictions("use_self_relations") )

            for (idUserEntity1, idUserEntity2, externalEntityRelationID, relation_type, partner_type) in listTupleIdUserEntity:
                if not user_entity_set.has_user_entity(idUserEntity1):
                    user_entity_set.addUserEntityId( idUserEntity = idUserEntity1 )
                    self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity1, partner_type))
                if not user_entity_set.has_user_entity(idUserEntity2):
                    user_entity_set.addUserEntityId( idUserEntity = idUserEntity2 )
                    self.outmethod("<user_entity id=\"%s\" type=\"%s\"/>" %(idUserEntity2, partner_type))
                    
                self.eEr_types_dict[externalEntityRelationID] = self.eEr_types_enum.get_letter(relation_type)
                self.uE_types_dict.setdefault(idUserEntity2,self.uE_types_enum.get_letter(partner_type))
                    
                inserted = user_entity_set.addUserEntityRelation(idUserEntity1 = idUserEntity1,
                                                                 idUserEntity2 = idUserEntity2, 
                                                                 externalEntityRelationID = externalEntityRelationID)

                if inserted:
                    self.outmethod("<user_entity_relation node1=\"%s\" node2=\"%s\" type=\"%s\" relation_id=\"%s\"/>" %(idUserEntity1, idUserEntity2, relation_type, externalEntityRelationID))

            OutBianaInterface.send_end_process_message()

        except:
            OutBianaInterface.send_error_notification( message = "New set not created. BIANA ERROR:", error = traceback.format_exc() )
            OutBianaInterface.send_end_process_message()
            return

        OutBianaInterface.send_end_process_message()

        if user_entity_set.getSize()==0:
            OutBianaInterface.send_error_notification(message="New set not created", error = "Any element has been found with given identifiers")
            return

        self.dictUserEntitySet[user_entity_set.id] = user_entity_set

        return user_entity_set



    def create_network(self, user_entity_set_id, level=0, include_relations_last_level = True, source_database_version_list=[], relation_type_list=[], relation_attribute_restriction_list=[], use_self_relations=True, expansion_attribute_list=[], expansion_relation_type_list=[], expansion_level=2, attribute_network_attribute_list=[], group_relation_type_list=[]):
        """
        Creates network of given user entity set adding relations of nodes as edges.
        ------
        user_entity_set_id: identifier of user entity set for which network will be created
        level: level of the network to be created, network will be expanded till that level
        include_relations_last_level: include relations between the nodes residing at the last level
        source_database_version_list: list of (biological database, version) tuples
        relation_type_list: type of the relations to be used in expansion 
        relation_attribute_restriction_list: tuples of (attribute, value) corresponding to restrictions to be applied on attributes of relations
        use_self_relations: include relations within the node itself
        expansion_attribute_list: tuples of (attribute, value_dictionary) corresponding to attributes to be used in relation inference between nodes based on shared attributes - value_dictionary is empty if attribute is not parameterizable
        expansion_relation_type_list: type of relations to be used in shared attribute based relation inference 
        expansion_level: number of relations (edges) to look further while inferring relations based on shared attributes
        attribute_network_attribute_list: tuples of (attribute, value) corresponding to attributes to be used while associating nodes with common attributes - value_dictionary is empty if attribute is not parameterizable
        group_relation_type_list: type of relations that are going to be treated as a group (like pathway, complex, cluster..)

        """
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]

        if user_entity_set.isNetworkCreated() == True:
            sys.stderr.write("Cannot create network on a set with created network. To expand use \"expand_user_entity_set\" command\n");
            OutBianaInterface.send_error_notification( message = "Network already created", error = "You cannot create a network on a set with created network, you can only expand it. In order to create a network with different parameters, create a new user entity set with selected nodes and create a network on it" )
            return

        OutBianaInterface.send_process_message("Creating network...")
        
        try:

            if ( len(expansion_attribute_list)==0 and len(expansion_relation_type_list)>0 ) or (len(expansion_attribute_list)>0 and len(expansion_relation_type_list)==0 ):
                raise ValueError("For expansions it is mandatory to specify attributes and relation types")

            user_entity_set.setIsNetworkCreated()

            dictAttributeToValues = self._convert_attribute_list_to_attribute_dictionary(relation_attribute_restriction_list, "embedded")

            # Fill the restrictions dictionary
            user_entity_set.addRestriction("relation_type_restrictions", relation_type_list)
            user_entity_set.addRestriction("externalDatabase_restrictions", source_database_version_list)
            user_entity_set.setRestriction("relation_attribute_restrictions", dictAttributeToValues)
            user_entity_set.setRestriction(restriction_type="use_self_relations",restriction=use_self_relations)
            user_entity_set.setRestriction(restriction_type="include_relations_last_level",restriction=include_relations_last_level)
            user_entity_set.setRestriction("group_relation_type", group_relation_type_list)

            #expansion related
            user_entity_set.addRestriction("expansionAttributesList", expansion_attribute_list)
            user_entity_set.addRestriction("expansionListRelationType", expansion_relation_type_list)
            user_entity_set.setRestriction("expansionLevel", expansion_level)

            #attribute network related
            user_entity_set.addRestriction("attributeNetworkList", attribute_network_attribute_list)

            levelStart = user_entity_set.get_level()

            if levelStart == level:
                if level == 0 and include_relations_last_level:
                    #self._getLastLevelRelations(user_entity_set)
                    self.expand_user_entity_set(user_entity_set_id = user_entity_set_id, is_last_level=True)
            else:
                for l in xrange(levelStart, level):
                    self.expand_user_entity_set(user_entity_set_id = user_entity_set_id)

                if include_relations_last_level:
                    self.expand_user_entity_set(user_entity_set_id = user_entity_set_id, is_last_level=True)
                    #self._getLastLevelRelations(user_entity_set)

        except:
            OutBianaInterface.send_error_notification( message = "Network already created", error = "You cannot create a network on a set with created network, you can only expand it. In order to create a network with different parameters, create a new user entity set with selected nodes and create a network on it" )
            
        OutBianaInterface.send_end_process_message()

        return


    def get_union_of_user_entity_set_list(self, user_entity_set_list, include_relations=False, new_user_entity_set_id=None): 
        """
        Does the union between of the user entities belonging to the user entities sets in the list and returns a new UserEntitySet
        ------
        user_entity_set_list: list of user entity set objects/ids to be combined
        include_relations: flag to whether or not include relations in union
        new_user_entity_set_id: identifier of new user entity set to be created as a result of the union
        """

        if new_user_entity_set_id is None:
            new_user_entity_set_id = self._get_next_uEs_id()

        # Check if the user_entity_set_list contains the objects or the names
        if isinstance(user_entity_set_list[0],str) or isinstance(user_entity_set_list[0],int):
            user_entity_set_list = [self.get_user_entity_set(x) for x in user_entity_set_list ]

        if len(user_entity_set_list) <2:   # TO CHECK IF THIS SHOULD BE DONE LIKE THIS... This line is here to prevent errors when only one userentityset is given
            new_user_entity_set = copy.deepcopy(user_entity_set_list[0])
            new_user_entity_set.id = new_user_entity_set_id
            self.dictUserEntitySet[new_user_entity_set.id] = new_user_entity_set
            return
                
        (listLevelSetId, listRelations) = user_entity_set_list[0].getUnionWithGivenUserEntitySet(user_entity_set_list[1], include_relations)

        new_user_entity_set = UserEntitySet.UserEntitySet(id = new_user_entity_set_id, setIdUserEntity=None, listRelations = listRelations, listLevelSetIdUserEntity = listLevelSetId)

        lenList = len(user_entity_set_list)
        
        for i in xrange(lenList):
            if i < 2:
                continue
            user_entity_set = user_entity_set_list[i]
            (listLevelSetId, listRelations) = new_user_entity_set.getUnionWithGivenUserEntitySet(user_entity_set, include_relations)
            new_user_entity_set = UserEntitySet.UserEntitySet(id = new_user_entity_set_id, setIdUserEntity=None, listRelations = listRelations, listLevelSetIdUserEntity = listLevelSetId) 

            
        self.dictUserEntitySet[new_user_entity_set.id] = new_user_entity_set

        self._send_complete_user_entity_set_info(user_entity_set=new_user_entity_set)

        return new_user_entity_set
    

    def get_intersection_of_user_entity_set_list(self, user_entity_set_list, include_relations=False, new_user_entity_set_id=None): 
        """
        Does the intersection between of the user Entities belonging to the user entities sets in the list and returns a new UserEntitySet
        ------
        user_entity_set_list: list of user entity set objects to be combined
        include_relations: flag to whether or not include relations in intersection
        new_user_entity_set_id: identifier of new user entity set to be created as a result of the intersection
        """

        if new_user_entity_set_id is None:
            new_user_entity_set_id  = self._get_next_uEs_id()

        # Check if the user_entity_set_list contains the objects or the names
        if isinstance(user_entity_set_list[0],str ) or isinstance(user_entity_set_list[0],int):
            user_entity_set_list = [self.get_user_entity_set(x) for x in user_entity_set_list ]

        if len(user_entity_set_list) <2:
            sys.stderr.write("It is necessary to have at least 2 sets in order to do the intersection")
            return


        (listLevelSetId, listRelations) = user_entity_set_list[0].getIntersectionWithGivenUserEntitySet(user_entity_set_list[1], include_relations)

        user_entity_set_new = UserEntitySet.UserEntitySet(id = new_user_entity_set_id, setIdUserEntity=None, listRelations = listRelations, listLevelSetIdUserEntity = listLevelSetId)

        lenList = len(user_entity_set_list)

        for i in xrange(lenList):
            if i < 2:
                continue
            user_entity_set = user_entity_set_list[i]
            (listLevelSetId, listRelations) = user_entity_set_new.getIntersectionWithGivenUserEntitySet(user_entity_set, include_relations)
            user_entity_set_new = UserEntitySet.UserEntitySet(id = new_user_entity_set_id, setIdUserEntity=None, listRelations = listRelations, listLevelSetIdUserEntity = listLevelSetId) 
            
        self.dictUserEntitySet[user_entity_set_new.id] = user_entity_set_new
        self._send_complete_user_entity_set_info(user_entity_set=user_entity_set_new)

        return user_entity_set_new

    def remove_selected_user_entities(self, user_entity_set_id):
        """
        Removes selected user entitites from given user entity set
        ------
        user_entity_set_id: identifier of user entity set for which selected user entities will be removed 
        """

        user_entity_set = self.get_user_entity_set(user_entity_set_id)

        selected = list(user_entity_set.getSelectedUserEntities()) # Needed to transform it to a list, because if not the selected set is modified during iteration

        [ user_entity_set.remove_node(current_node) for current_node in selected ]
        
        self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content="<remove_user_entities ids=\"%s\"/>" %",".join(map(str,selected)))))

        return

    def remove_selected_relations(self, user_entity_set_id):
        """
        Removes selected relations from given user entity set
        ------
        user_entity_set_id: identifier of user entity set for which selected user entities will be removed 
        """

        user_entity_set = self.get_user_entity_set(user_entity_set_id)
        selected = list(user_entity_set.getSelectedUserEntityRelations()) # Needed to transform it to a list, because if not the selected set is modified during iteration
        
        #[ user_entity_set.remove_external_entity_relation(current_edge) for current_edge in selected ]

        xml = []
        for current_eErID in selected:
            participants = user_entity_set.get_external_entity_relation_participants(current_eErID)
            for x in xrange(len(participants)):
                for y in xrange(x):
                    xml.append( "<remove_user_entity_relations ids=\"%s,%s,%s\"/>" %(participants[x],participants[y],self.eEr_types_enum.get(self.eEr_types_dict[current_eErID])) )
            user_entity_set.remove_external_entity_relation(current_eErID)
                    
        
        self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content="".join(xml))))
        return


    ######################
    ### OUTPUT METHODS ###
    ######################

    def output_ontology(self, ontology_object, additional_attributes = [], out_method=None):
        """
        Outputs the ontology in XML Format
        ------
        ontology_name: Ontology to be loaded
        additional_attributes: Not implemented yet. TO DO!!! (Attribute "id" is the visual ID)
        out_method: output method to be used if None overwritten by instance default output method
        """
        
        if out_method is None:
            out_method = self.outmethod

        eEID_list = ontology_object.get_all_external_entity_ids()

        out_method(ontology_object.get_xml())

        return

    def output_external_entity_relation_participant_details(self, external_entity_relation_id, out_method=None ):
        """
        Outputs details of participants involved in a given external entry relation 
        ------
        external_entity_relation_id: relation identifier for which participant details will be outputted
        out_method: output method to be used if None overwritten by instance default output method

        Called by BianaGUI/showTableDialog.java via output_external_entity_relation_details method
        """
        # TO FINISH

        if out_method is None:
            out_method = self.outmethod

        columns = ["External Entity ID", "Source Database"]

        values = []

        eEr_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = [external_entity_relation_id] )

        eEr_obj = eEr_dict[external_entity_relation_id]

        eE_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = eEr_obj.get_participant_external_entity_ids_list() )
        
        for current_eE in eE_dict.values():
            
            current_values = [current_eE.get_id()]
            current_values.append( "%s" %self.dbAccess.get_external_database( database_id = current_eE.get_source_database() ) )

            values.append(current_values)
        
        out_method(output_utilities.get_html_table(columns,values,attributes=[("title","Details for relation %s" %external_entity_relation_id)]))
        return
        

    def output_external_entity_relation_details(self, out_method=None, external_entity_relation_id_list=[], node_attributes=[], relation_attributes=[], participant_attributes=[]):
        """
        Outputs details of external entry relations with given identifiers
        ------
        external_entity_relation_id_list: list of relation identifiers for which details will be outputted
        node_attributes: attributes of user entities connected by these relations for which information will be fetched
        relation_attributes: attributes of external entity relations for which information will be fetched
        participant_attributes: attributes of external entity relation participants for which information will be fetched
        out_method: output method to be used if None overwritten by instance default output method
        """

        if out_method is None:
            out_method = self.outmethod

        #columns = ["External Entity Relation ID","External Database","Relation Type", "Participants"]
        #columns.extend([ "Participant 1 %s" %x for x in node_attributes])
	#columns.extend([ "Participant 1 %s" %x for x in participant_attributes])
        #columns.extend([ "Participant 2 %s" %x for x in node_attributes])
	#columns.extend([ "Participant 2 %s" %x for x in participant_attributes])
        #columns.extend([ "Relation %s" %x for x in relation_attributes])

        # New columns
        columns = ["External Entity Relation ID", "External Database", "Relation Type", "Number of participants", "Participant External Entity ID"]
        columns.extend(relation_attributes)
        columns.extend(node_attributes)
        columns.extend(participant_attributes)

        values = []

        command = "output_external_entity_relation_participant_details( external_entity_relation_id = ? )"

        rowIDs = []

        eEr_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = external_entity_relation_id_list,
                                                             attribute_list = node_attributes,
                                                             relation_attribute_list = relation_attributes,
                                                             participant_attribute_list = participant_attributes )

        
        value_separator = ","

        for current_eEr in eEr_dict.values():

            eE_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = current_eEr.get_participant_external_entity_ids_list(), attribute_list = node_attributes )

            for current_participant in current_eEr.get_participant_external_entity_ids_list():

                rowIDs.append(current_eEr.get_id())
                current_values = [ current_eEr.get_id() ]
                current_values.append( "%s" %self.dbAccess.get_external_database( database_id = current_eEr.get_source_database()) )
                current_values.append( current_eEr.get_relation_type() )
                current_values.append( "%s" %len(current_eEr.get_participant_external_entity_ids_list()) )
                current_values.append( "%s" %current_participant )

                for current_attribute in relation_attributes:
                    inner_values = {}
                    for y in current_eEr.get_attribute(attribute_identifier = current_attribute):
                        n = inner_values.setdefault(str(y.value), 0)
                        inner_values[str(y.value)] = n+1
                    if len(inner_values) > 0:
                        current_values.append(value_separator.join([ "%s(%s)" % (i,j) for i,j in inner_values.iteritems()]))
                    else:
                        current_values.append("-")

                for current_attribute in node_attributes:
                    #defined_attributes = self.get_defined_node_attributes(user_entity_set.id, current_participant, current_attribute, output_1_value_per_attribute, substitute_node_attribute_if_not_exists = substitute_node_attribute_if_not_exists, return_set = True)
                    attribute_values = [ str(y.value).replace("\n"," ") for y in eE_dict[current_participant].get_attribute(attribute_identifier=current_attribute) ]
                    
                    temp_str = value_separator.join(attribute_values)
                    if temp_str == "":
                        current_values.append("-")
                    else:
                        current_values.append( temp_str )

                for current_attribute in participant_attributes:
                    participant_attribute_values = current_eEr.get_participant_attribute( participantExternalEntityID =  current_participant,
                                                                                    attribute_identifier = current_attribute )

                    temp_str = value_separator.join(participant_attribute_values)
                    if temp_str == "":
                        current_values.append("-")
                    else:
                        current_values.append( temp_str )

                values.append(current_values)

        out_method(output_utilities.get_html_table(columns,values,rowIDs,attributes=[("title","External Entity Relation Details"),("command",command),("session",self.sessionID)]))
        return



    def get_defined_node_attributes(self, user_entity_set_id, user_entity_id, attribute, output_1_value_per_attribute, return_set=False, substitute_node_attribute_if_not_exists=False):
        """
        Gets user defined and BIANA defined node external names - adds substitution functionality to helper
        ------
        user_entity_set_id: identifier for the user entity set information for whose user entity will be fetched 
        user_entity_id: user_entity_id for the node for wich defined attributes are obtained
        attribute: attribute name corresponding to the attributes to be outputed
        output_1_value_per_attribute: Boolean. Defines wether 1 or multiple values are outputed per each attribute
        """

        defined_attributes = self.get_defined_node_attributes_helper(user_entity_set_id, user_entity_id, attribute, output_1_value_per_attribute)

        if substitute_node_attribute_if_not_exists:
            current_attribute = attribute
            copy_list = copy.deepcopy(self.substitution_list)
            while len(defined_attributes)==0:
                if current_attribute in copy_list:
                    copy_list.remove(current_attribute)
                if len(copy_list) == 0:
                    break
                current_attribute = copy_list[0]
                defined_attributes = self.get_defined_node_attributes_helper(user_entity_set_id, user_entity_id, current_attribute, output_1_value_per_attribute)
                defined_attributes = [ current_attribute+":"+i for i in defined_attributes ]

        if return_set: 
            defined_attributes = set(defined_attributes )

        return defined_attributes

 
    def get_defined_node_attributes_helper(self, user_entity_set_id, user_entity_id, attribute, output_1_value_per_attribute):
        """
        Gets user defined and BIANA defined user_entity_id external names
        ------
        user_entity_set_id: identifier for the user entity set information for whose user entity will be fetched 
        user_entity_id: user_entity_id for the user_entity_id for wich defined attributes are obtained
        attribute: attribute name corresponding to the attributes to be outputed
        output_1_value_per_attribute: Boolean. Defines wether 1 or multiple values are outputed per each attribute
        """

        if False:      # Change to true for debugging purposes
            sys.stderr.write("* Entering in get_defined_node_attributes_helper():\n")
            sys.stderr.write("*\tuser_entity_set_id           : %s\n" %user_entity_set_id)
            sys.stderr.write("*\tnode                        : %s\n" %node)
            sys.stderr.write("*\tattribute                   : %s\n" %attribute)
            sys.stderr.write("*\toutput_1_value_per_attribute: %s\n" %output_1_value_per_attribute)

	def get_more_represented(aList, preference_function=None): #lambda x: x):
            """
                Helper method for get_defined_node_attributes_helper
		aList: list of values for an attribute
		preference_function: if not None, among more represented values, the value giving higher return value with the preference_function will be selected
            """
            #aList.sort()
            aDict = {}
            for i in aList:
                if not aDict.has_key(i): aDict[i] = 0
                aDict[i] += 1

            more_represented_list = []
	    n_more_represented = 0
            for k in aDict.keys():
		n = aDict[k]
                if n > n_more_represented: 
		    more_represented_list = [k]
		    n_more_represented = n
		elif n == n_more_represented:
		    more_represented_list.append(k)

	    more_represented = None
	    if preference_function is None and len(more_represented_list)>0:
		more_represented = more_represented_list[0]
	    else:
		more_represented_score = -1
		for i in more_represented_list:
		    score = preference_function(i)
		    if score > more_represented_score: # Select better scoring value
			more_represented_score = score
			more_represented = i

            return more_represented

        attribute = attribute.lower() 
	user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        new_values = []

        attribute_values_dict = self.dbAccess.get_user_entity_attributes( unification_protocol_name = self.unification_protocol_name,
                                                                          listUserEntityID = [user_entity_id],
                                                                          attribute_identifier = attribute )

        attribute_values = attribute_values_dict.setdefault(user_entity_id, [])

        #if output_1_value_per_attribute is False or ExternalEntityAttribute.isIdentifierType(attribute, self.dbAccess.biana_database) is False:  # It should be done in other way... it should not access biana_database from dbAccess from here...
        if output_1_value_per_attribute is False: 
            new_values.extend(attribute_values)
        else:
            user_attribute_values = []
            # Check attribute node names provided by user
            if user_entity_set.userProvidedExternalIdsDict.has_key(attribute):
                temp_attribute_values = set([i.lower() for i in attribute_values]).intersection(set([i.lower() for i in user_entity_set.userProvidedExternalIdsDict[attribute]]))
                user_attribute_values = [user_entity_set.userProvidedExtIdsLowerDict[attribute][i] for i in temp_attribute_values]

            # Check attribute user_entity_id names previously chosen by BianaSessionManager
            elif user_entity_set.selectedExternalIdsDict.has_key(attribute):
                temp_attribute_values    = set([i.lower() for i in attribute_values]).intersection(set([i.lower() for i in user_entity_set.selectedExternalIdsDict[attribute]]))
                user_attribute_values = [user_entity_set.selectedExtIdsLowerDict[attribute][i] for i in temp_attribute_values]

            if len(user_attribute_values) > 0: 
                new_values = map(str, user_attribute_values)
            else:
		if attribute=="proteinsequence":
		    more_represented = get_more_represented(aList = attribute_values, preference_function = len)
		else:
		    more_represented = get_more_represented(attribute_values)
                user_entity_set._update_selected_external_ids_dict(attribute, more_represented)
                if more_represented is not None:
                    new_values = [ more_represented ]

        #print new_values
        return new_values

    def get_user_entity_set_attribute_network(self, user_entity_set_id, node_attribute):
	"""
	Gets a network of node attributes instead of user entity nodes
        ------
        Deprecated - use get_defined_node_attributes instead
	"""

        sys.stderr.write("DEPRECATED!\nMust use BianaSessionManager.get_defined_node_attributes()\n") # line added by Joan. 

	user_entity_set = self.dictUserEntitySet[user_entity_set_id]
	nodes = user_entity_set.get_user_entity_ids()
	edges = user_entity_set.getRelations()

	added_nodes = set()

	#print len(edges), " edges"

	new_network = graph_utilities.create_graph(allow_multi_edges=False)

	new_edges_set = set()

	for current_edge in edges:
			

		uEobj1 = self.get_user_entity(user_entity_id=current_edge[0])
		uEobj2 = self.get_user_entity(user_entity_id=current_edge[1])

		added_nodes.add(current_edge[0])
		added_nodes.add(current_edge[1])

		eE_dict1 = self.dbAccess.get_external_entities_dict( externalEntityIdsList = uEobj1.get_externalEntitiesIds_set(),
									attribute_list = [node_attribute] )
		
		eE_dict2 = self.dbAccess.get_external_entities_dict( externalEntityIdsList = uEobj2.get_externalEntitiesIds_set(),
                                                                        attribute_list = [node_attribute] )

		for current_eE1 in uEobj1.get_externalEntitiesIds_set():
			for current_attribute1 in eE_dict1[current_eE1].get_attribute(attribute_identifier=node_attribute):
				for current_eE2 in uEobj2.get_externalEntitiesIds_set():
		                        for current_attribute2 in eE_dict2[current_eE2].get_attribute(attribute_identifier=node_attribute):
						new_edges_set.add((current_attribute1.value.lower(),current_attribute2.value.lower()))
	
	new_network.add_edges_from(new_edges_set)
	for current_node in nodes:
		if current_node not in added_nodes:
			uEobj = self.get_user_entity(user_entity_id=current_node)
			eE_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = uEobj.get_externalEntitiesIds_set(),
                                                                        attribute_list = [node_attribute] )
			for current_eE in uEobj.get_externalEntitiesIds_set():
				for current_attribute in eE_dict[current_eE].get_attribute(attribute_identifier=node_attribute):
					new_network.add_node(current_attribute.value.lower())

	return new_network


    def output_user_entity_set_group(self, user_entity_set_id, group_ids):
        """
        "group_ids" must be the external entity relation that represents the group
        """

        OutBianaInterface.send_process_message("Getting information from database...")


        # First, check the parameters are correct
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]

        info = []

        for group_id in group_ids:
            if user_entity_set.has_group(group_id) is False:
                OutBianaInterface.send_error_notification("Group %s not in %s" %(group_id, user_entity_set_id))
            else:
                external_entity_relation_id = int(group_id)

                external_entity_relation_obj = self.dbAccess.get_external_entities_dict( externalEntityIdsList = [external_entity_relation_id],
                                                                                         attribute_list = [],
                                                                                         relation_attribute_list = [] )[external_entity_relation_id]

                external_database_obj = self.dbAccess.get_external_database( database_id = external_entity_relation_obj.get_source_database() )
                
            
                external_entity_relation_obj = self.dbAccess.get_external_entities_dict( externalEntityIdsList = [external_entity_relation_id],
                                                                                         attribute_list = external_database_obj.get_valid_external_entity_relation_attribute_type(),
                                                                                         relation_attribute_list = [] )[external_entity_relation_id]
                info.append(external_entity_relation_obj.__str__())

            OutBianaInterface.send_info_message("<html><br><br>".join(info)+"</html>")
            
        # Then, print for the relation that is representing the group all its attributes
        
	OutBianaInterface.send_end_process_message()

	return

        
    def output_user_entity_set_network_in_sif_format(self, user_entity_set_id, output_path = "./", output_prefix = "", node_attributes = [], participant_attributes = [], relation_attributes=[], include_tags=True, output_1_value_per_attribute=True, only_selected=False):
	"""
        Outputs information of the network of a given user entity set in sif format. User entity id will be taken as the node identifier.
	------
        output_path: directory in which Cytoscape attribute files will be generated
        output_prefix: prefix to be added to the begining of Cytoscape attribute files
	"""
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        unconnected_nodes = user_entity_set.get_unconnected_nodes()
        if only_selected:
            nodes = user_entity_set.getSelectedUserEntities()
            unconnected_nodes &= nodes 
            edges = user_entity_set.getRelationsOfSelectedUserEntities()
        else:
            edges = user_entity_set.getRelations()

	import os
	if not os.path.exists(os.path.abspath(output_path)):
	    sys.stderr.write("output_user_entity_set_network: given output path does not exist\n")
	    return
	else:
	    sif_file = open("%s/%s.sif" % (os.path.abspath(output_path), output_prefix), 'w')
	    source_file = open("%s/%s_%s.eda" % (os.path.abspath(output_path), output_prefix, "source"), 'w')
	    source_file.write("%s\n" % "SourceDB")
	    node_attribute_files = {}
	    for current_attribute in node_attributes:
		node_attribute_files[current_attribute] = open("%s/%s_%s.noa" % (os.path.abspath(output_path), output_prefix, current_attribute), 'w')
		node_attribute_files[current_attribute].write("%s\n" % current_attribute)
	    for current_tag in user_entity_set.get_all_user_entity_tags():
		node_attribute_files["tag_%s" % current_tag] = open("%s/%s_tag_%s.noa" % (os.path.abspath(output_path), output_prefix, current_tag), 'w')
		node_attribute_files["tag_%s" % current_tag].write("tag_%s\n" % current_tag)
	    edge_attribute_files = {}
	    for current_attribute in relation_attributes:
		edge_attribute_files[current_attribute] = open("%s/%s_%s.eda" % (os.path.abspath(output_path), output_prefix, current_attribute), 'w')
		edge_attribute_files[current_attribute].write("%s\n" % current_attribute)
	    for current_tag in  user_entity_set.get_all_user_entity_relation_tags():
		edge_attribute_files["tag_%s" % current_tag] = open("%s/%s_tag_%s.eda" % (os.path.abspath(output_path), output_prefix, current_tag), 'w')
		edge_attribute_files["tag_%s" % current_tag].write("tag_%s\n" % current_tag)
	    participant_attribute_files = {}
	    for current_attribute in participant_attributes:
		participant_attribute_files[current_attribute] = open("%s/%s_%s.eda" % (os.path.abspath(output_path), output_prefix, current_attribute), 'w')
		participant_attribute_files[current_attribute].write("%s\n" % current_attribute)

        ## Unconnected nodes
	for current_node in unconnected_nodes :
	    sif_file.write("%s\n" %current_node)
	    ## Node attributes
	    for current_attribute in node_attributes:
		for current_value in self.get_defined_node_attributes(user_entity_set.id, current_node, current_attribute, output_1_value_per_attribute, substitute_node_attribute_if_not_exists = False, return_set = True):
		    node_attribute_files[current_attribute].write("%s = %s\n" % (current_node, current_value))
	    ## Node tags
	    if include_tags:
		for current_tag in user_entity_set.get_all_user_entity_tags():
		    node_attribute_files["tag_%s" % current_tag].write("%s = %s\n" % (current_node, user_entity_set.has_tag(current_node, current_tag)))

	## Edges and connected nodes
	eEr_dict = {}
	included_nodes = set()
        for current_edge in edges:
	    eEr_dict = {} #! testing memory consumption (emre)
            eErIDs_list = user_entity_set.get_external_entity_relation_ids(current_edge[0], current_edge[1])
            # Update dictionary
            for current_eErID in eErIDs_list:
                if not eEr_dict.has_key(current_eErID):
                    eEr_dict.update(self.dbAccess.get_external_entities_dict( externalEntityIdsList = [current_eErID],
                                                                              attribute_list = node_attributes,
                                                                              relation_attribute_list = relation_attributes,
                                                                              participant_attribute_list = participant_attributes ))
            types = set()
            types2eErIDs_list = {}

            for current_eEr_id in eErIDs_list:
                types.add(self.eEr_types_dict[current_eEr_id])
                types2eErIDs_list.setdefault(self.eEr_types_dict[current_eEr_id], []).append(current_eEr_id)

	    ## Edges
	    # Distinguish between relation type
	    for current_relation_type in types:
		sif_file.write("%s (%s) %s\n" % (current_edge[0], self.eEr_types_enum.get(current_relation_type), current_edge[1]))
		inner_values = {}
		for current_uErID in types2eErIDs_list[current_relation_type]:
		    if current_uErID > 0:
			source = eEr_dict[current_uErID].get_source_database()
			n = inner_values.setdefault(source, 0)
			inner_values[source] = n+1
			## Relation attributes
			for current_attribute in relation_attributes:
			    [ edge_attribute_files[current_attribute].write("%s (%s) %s = %s\n" % (current_edge[0], self.eEr_types_enum.get(current_relation_type), current_edge[1], str(y.value))) for y in eEr_dict[current_uErID].get_attribute(attribute_identifier = current_attribute) ]
			## Participant attributes
			for current_participant in [current_edge[0],current_edge[1]]:
			    for current_attribute in participant_attributes:
				[ participant_attribute_files[current_attribute].write("%s (%s) %s = %s\n" % (current_edge[0], self.eEr_types_enum.get(current_relation_type), current_edge[1], str(current_attr.value))) for current_eE in self.get_user_entity(user_entity_id = current_participant).get_externalEntitiesIds_set() for current_attr in eEr_dict[current_uErID].get_participant_attribute( participantExternalEntityID = current_eE, attribute_identifier = current_attribute ) ]
			## Relation tags
			if include_tags:
			    for current_tag in user_entity_set.get_all_user_entity_relation_tags():
				edge_attribute_files["tag_%s" % current_tag].write("%s (%s) %s = %s\n" % (current_edge[0], self.eEr_types_enum.get(current_relation_type), current_edge[1], user_entity_set.relation_has_tag(current_uErID, current_tag) ))
		[ source_file.write("%s (%s) %s = %s(%s)\n" % (current_edge[0], self.eEr_types_enum.get(current_relation_type), current_edge[1], self.dbAccess.get_external_database(i).get_name(),j)) for i,j in inner_values.iteritems() ]

	    ## Connected node attributes
	    for current_participant in [current_edge[0],current_edge[1]]:
		## Node attributes
		if current_participant not in included_nodes:
		    for current_attribute in node_attributes:
			[ node_attribute_files[current_attribute].write("%s = %s\n" % (current_participant, current_value)) for current_value in self.get_defined_node_attributes(user_entity_set.id, current_participant, current_attribute, output_1_value_per_attribute, substitute_node_attribute_if_not_exists = False, return_set = True) ]
		    included_nodes.add(current_participant)
	
		## Node tags
		if include_tags:
		    for current_tag in user_entity_set.get_all_user_entity_tags():
			node_attribute_files["tag_%s" % current_tag].write("%s = %s\n" % (current_participant, user_entity_set.has_tag(current_participant, current_tag)))
	    
	sif_file.close()
	source_file.close()
	for current_attribute in node_attributes:
	    node_attribute_files[current_attribute].close()
	for current_tag in user_entity_set.get_all_user_entity_tags():
	    node_attribute_files["tag_%s" % current_tag].close()
	for current_attribute in relation_attributes:
	    edge_attribute_files[current_attribute].close()
	for current_tag in  user_entity_set.get_all_user_entity_relation_tags():
	    edge_attribute_files["tag_%s" % current_tag].close()
	for current_attribute in participant_attributes:
	    participant_attribute_files[current_attribute].close()

        OutBianaInterface.send_end_process_message()
	return
    

    def output_user_entity_set_network(self, user_entity_set_id, out_method=None, node_attributes = [], participant_attributes = [], relation_attributes=[], allowed_relation_types="all", include_participant_tags=True, include_relation_tags=True, include_relation_ids=True, include_participant_ids=True, include_relation_type=True, include_relation_sources=True, output_1_value_per_attribute=True, output_format="xml", value_seperator=", ", only_selected=False, include_command_in_rows=False, substitute_node_attribute_if_not_exists=False, include_unconnected_nodes=True):
        """
        Outputs information of the network of a given user entity set
        ------
        output_1_value_per_attribute: Boolean. Defines wether 1 or multiple values are outputed per each attribute
        output_format: format for the output used in case format is "table"; can be "tabulated" or "xml"
        include_relation_ids: Boolean to whether display or not relation identifiers
        include_participant_ids: Boolean to whether display or not relation participant identifiers
        include_relation_type: Boolean to whether display or not types of relations
        include_relation_sources: Boolean to whether display or not relation sources
        include_participant_tags: Boolean to whether display or not tags of participants
        include_relation_tags: Boolean to whether display or not tags of relations
        value_seperator: string to seperate consequitive values in the same column
        only_selected: Boolean to decide whether to output only selected nodes or all nodes (and their interactions)
        include_command_in_rows: Include the command to output individual relation information at each row
        substitute_node_attribute_if_not_exists: In case the node does not have a value for a given attribute (s.t. uniprotaccession) this flag make it possible to output another attribute (e.g. geneid) in the same column indicated as attribute:value (e.g. geneid:123123)
        include_unconnected_nodes: Boolean to whether display or not unconnected nodes
        """

        OutBianaInterface.send_process_message("Getting information from database...")
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        
        unconnected_nodes = user_entity_set.get_unconnected_nodes()
        if only_selected:
            nodes = user_entity_set.getSelectedUserEntities()
            unconnected_nodes &= nodes 
            edges = user_entity_set.getRelationsOfSelectedUserEntities()
        else:
            edges = user_entity_set.getRelations()

        if out_method is None:
            out_method = self.outmethod

        # Check excluded relation types
        excluded_relation_types = {}

        if not allowed_relation_types == "all":
            [ excluded_relation_types.update({relation_type:None}) for relation_type in self.dbAccess.get_valid_external_entity_relation_types() if not relation_type in allowed_relation_types ]

	columns = []
	if include_relation_ids:
	    columns.append("Relation IDs")
	if include_relation_type:
	    columns.append("Relation Types")
	if include_relation_tags:
	    columns.append("Relation Tags")
	if include_participant_ids:
	    columns.append("Participant 1 User Entity")
	columns.extend([ "Participant 1 %s" %x for x in node_attributes])
	columns.extend([ "Participant 1 %s" %x for x in participant_attributes])
	if include_participant_tags:
	    columns.append("Participant 1 Tags")
	if include_participant_ids:
	    columns.append("Participant 2 User Entity")
	columns.extend([ "Participant 2 %s" %x for x in node_attributes])
	columns.extend([ "Participant 2 %s" %x for x in participant_attributes])
	if include_participant_tags:
	    columns.append("Participant 2 Tags")
	columns.extend([ "Relation %s" %x for x in relation_attributes])
	if include_relation_sources:
	    columns.append("Relation Source Databases")

	if include_command_in_rows:
	    #command = "output_external_entity_relation_details(external_entity_relation_id_list=[?])"
            command = "output_external_entity_relation_details(external_entity_relation_id_list=[?],node_attributes=%s,relation_attributes=%s,participant_attributes=%s)" %(node_attributes, relation_attributes, participant_attributes)
	else:
	    command = ""

	if output_format == "xml":
	    out_method(output_utilities.get_html_table_header( columns = columns, attributes = [ ("id", "user_entity_set_edges"), ("title","User Entity Set Network Details"),("command",command),("session",self.sessionID) ] ) )
	elif output_format == "tabulated":                                            # line added by Joan        # REQUIRED TO TAG THE PRINTED COLUMNS 
	    out_method( output_utilities.get_tabulated_table(values = [columns]) )    # line added by Joan        # REQUIRED TO TAG THE PRINTED COLUMNS 

        ## Unconnected node attributes  
        if include_unconnected_nodes:
	    for current_node in unconnected_nodes:
		new_values = []
		if include_relation_ids:
		    new_values.append("-")
		if include_relation_type:
		    new_values.append("-")
		if include_relation_tags:
		    new_values.append("-")
		if include_participant_ids:
		    new_values.append(current_node)

		current_rowID = ""  # Added by Joan.  
		for current_attribute in node_attributes:
		    defined_attributes = self.get_defined_node_attributes(user_entity_set.id, current_node, current_attribute, output_1_value_per_attribute, substitute_node_attribute_if_not_exists = substitute_node_attribute_if_not_exists, return_set = True)
		    temp_str = value_seperator.join(defined_attributes)
		    if temp_str == "":
			new_values.append("-")
		    else:
			new_values.append( temp_str ) #defined_attributes) )                                                                                                                                 
		[ new_values.append("-") for x in participant_attributes ]
		if include_participant_tags:
		    tags = user_entity_set.get_user_entity_tags(current_node)
		    if len(tags)==0:
			new_values.append("-")
		    else:
			new_values.append(value_seperator.join(tags))

		if include_participant_ids:
		    new_values.append("-")

		[ new_values.append("-") for x in node_attributes ]
		[ new_values.append("-") for x in participant_attributes ]

		if include_participant_tags:
		    new_values.append("-")

		for current_attribute in relation_attributes:
		    new_values.append("-")

		if include_relation_sources:
		    new_values.append("-")

		if output_format == "xml"        : out_method(output_utilities.append_html_table_values( rowIDs = [current_rowID], values = [new_values] )) # Added by Joan. 
		elif output_format == "tabulated": out_method( output_utilities.get_tabulated_table(values = [new_values]) )                                # Added by Joan. 
		else                             : raise ValueError("output_format is not valid. Valid output formats are ['xml', 'tabulated']\n")         # Added by Joan.   

        # PRINT THE EDGES
        # Get the objects themselves   
        eEr_dict = {}

        for current_edge in edges:
	    eEr_dict = {} #! otherwise memory could be a problem because of this dictionary (emre)
            eErIDs_list = user_entity_set.get_external_entity_relation_ids(current_edge[0], current_edge[1])
            # Update dictionary
            for current_eErID in eErIDs_list:
                if not eEr_dict.has_key(current_eErID):
                    eEr_dict.update(self.dbAccess.get_external_entities_dict( externalEntityIdsList = [current_eErID],
                                                                              attribute_list = node_attributes,
                                                                              relation_attribute_list = relation_attributes,
                                                                              participant_attribute_list = participant_attributes ))
            types = set()
            types2eErIDs_list = {}

            for current_eEr_id in eErIDs_list:
                types.add(self.eEr_types_dict[current_eEr_id])
                types2eErIDs_list.setdefault(self.eEr_types_dict[current_eEr_id], []).append(current_eEr_id)

            for current_relation_type in types:
                if excluded_relation_types.has_key(current_relation_type): continue
                new_values = []
                current_rowID = ""

                # RELATION HEADER 
		relation_id_list = map(str,types2eErIDs_list[current_relation_type])
		current_rowID = ",".join(relation_id_list)
                if include_relation_ids:
                    new_values.append( value_seperator.join(relation_id_list) )   # Relation ID
                    
		if include_relation_type:
		    new_values.append(self.eEr_types_enum.get(current_relation_type))

		if include_relation_tags:
		    tags = []
		    relation_id_list = types2eErIDs_list[current_relation_type]
		    [ tags.extend(user_entity_set.get_external_entity_relation_tags(i)) for i in relation_id_list ]
		    if len(tags)==0:
			new_values.append("-")
		    else:
			new_values.append(value_seperator.join( tags ))
		
		for current_participant in [current_edge[0],current_edge[1]]:
		    if include_participant_ids:
			new_values.append(current_participant)
		    for current_attribute in node_attributes:
			defined_attributes = self.get_defined_node_attributes(user_entity_set.id, current_participant, current_attribute, output_1_value_per_attribute, substitute_node_attribute_if_not_exists = substitute_node_attribute_if_not_exists, return_set = True)
			temp_str = value_seperator.join(defined_attributes)
			if temp_str == "":
			    new_values.append("-")
			else:
			    new_values.append( temp_str )

		    for current_attribute in participant_attributes:
			defined_attributes = [ str(current_attr.value) for current_eE in self.get_user_entity(user_entity_id = current_participant).get_externalEntitiesIds_set()
					       for current_eEr_id in current_edge[2] for current_attr in eEr_dict[current_eEr_id].get_participant_attribute( participantExternalEntityID =  current_eE,
																			     attribute_identifier = current_attribute ) ]
			temp_str = value_seperator.join(defined_attributes)
			if temp_str == "":
			    new_values.append("-")
			else:
			    new_values.append( temp_str )

		    if include_participant_tags:
			tags = user_entity_set.get_user_entity_tags(current_participant)
			if len(tags)==0:
			    new_values.append("-")
			else:
			    new_values.append(value_seperator.join(tags))
            
		# RELATION SPECIFIC 
		for current_attribute in relation_attributes:
		    inner_values = {}
		    for current_uErID in types2eErIDs_list[current_relation_type]:   # JAVI TO CHECK: WHY LOOP types2eErIDs???
			if current_uErID > 0:
			    for y in eEr_dict[current_uErID].get_attribute(attribute_identifier = current_attribute):
				n = inner_values.setdefault(str(y.value), 0)
				inner_values[str(y.value)] = n+1
		    if len(inner_values) > 0:
			new_values.append(value_seperator.join([ "%s(%s)" % (i,j) for i,j in inner_values.iteritems()]))
		    else:
			new_values.append("-")

		if include_relation_sources:
		    inner_values = {}
		    for current_uErID in types2eErIDs_list[current_relation_type]:
			if current_uErID > 0:
			    source = eEr_dict[current_uErID].get_source_database()
			    n = inner_values.setdefault(source, 0)
			    inner_values[source] = n+1
		    new_values.append(value_seperator.join([ "%s(%s)" % (self.dbAccess.get_external_database(i).get_name(),j) for i,j in inner_values.iteritems()]))

		if output_format == "xml"        : out_method(output_utilities.append_html_table_values( rowIDs = [current_rowID], values = [new_values] ))
		elif output_format == "tabulated": out_method( output_utilities.get_tabulated_table(values = [new_values]) )
		else                             : raise ValueError("output_format is not valid. Valid output formats are ['xml', 'tabulated']\n")

        if output_format == "xml":
            out_method(output_utilities.get_html_table_foot())
	OutBianaInterface.send_end_process_message()
        return
    

    def output_user_entity_details(self, user_entity_set, user_entity_id_list, out_method = None, attributes = [], include_level_info = True, include_degree_info=True, include_tags_info=True, include_tags_linkage_degree_info=[], substitute_node_attribute_if_not_exists=False, output_1_value_per_attribute=True, output_format="tabulated", include_command_in_rows=False):
        """
        Outputs given user entity. Called by output_user_entity_set_details.
        ------
        user_entity_set: Instace of the user entity set, user entities of which will be outputted
        user_entity_id_list: list of user entities 
        out_method: method to which output information will be directed
        attributes: list of attribute names corresponding to the attributes to be outputed
        include_level_info: Boolean. Defines inclussion of network level information in the output. 
        include_degree_info: Boolean. Defines inclussion of node degree information in the output.
        inlcude_tags_info: Boolean. Defines inclussion of node tags information in the output.
        include_tags_linkage_degree_info: Calculate linkage degree considering only tagged nodes where tags given as a list
        output_1_value_per_attribute: Boolean. Defines wether 1 or multiple values are outputed per each attribute
        output_format: format for the output. Either "tabulated" or "xml"
        include_command_in_rows: Include the command to output individual user entity information at each row
        """

        OutBianaInterface.send_process_message("Getting information from database...")

        # Define columns
        columns = ["User Entity ID"]
        tld = {}
        if include_level_info:
            columns.append("Level")
        if include_degree_info:
            columns.append("Degree")
        if include_tags_info:
            columns.append("Tags")
        for current_tag_linkage_degree in include_tags_linkage_degree_info:
            columns.append("%s linkage degree" %current_tag_linkage_degree)
            tld[current_tag_linkage_degree] = user_entity_set.getTagLinkageDegree(tag = current_tag_linkage_degree)
        
        # Add default identifier for the user entity
        # columns.append("Default Attribute")
        
        

        columns.extend([ str(x) for x in attributes])

        if include_command_in_rows:
            #command = "output_external_entity_details(user_entity_id_list=[%s], attributes=[\'%s\'])" %( ",".join( [ str(uE_id) for uE_id in user_entity_id_list ] ), "\',\'".join(attributes))
            #command = "output_external_entity_details(user_entity_id_list=[?], attributes=[\'%s\'])" %("\',\'".join(attributes))
            command = "output_external_entity_details(user_entity_id_list=[?]"
            if len(attributes)>0:
                command += ", attributes=[\'%s\'])" %("\',\'".join(attributes))
            else:
                command += ")"
        else:
            command = "" #"output_external_entity_details(user_entity_id_list=[...], attributes=[\'%s\'])" %(",".join(attributes))

        if out_method is None:
            out_method = self.outmethod

        if output_format == "xml":
            out_method(output_utilities.get_html_table_header( columns = columns, attributes = [("id", "user_entity_set_nodes"), ("command",command),("title","User Entity Set Details"),("session",self.sessionID)] ) )

        elif output_format == "tabulated":                                            # line added by Joan        # REQUIRED TO TAG THE PRINTED COLUMNS
            out_method( output_utilities.get_tabulated_table(values = [columns]) )     # line added by Joan        # REQUIRED TO TAG THE PRINTED COLUMNS


##
##         # ----------------------------------------------
##         # Modified by Joan
##         # Dict not needed any more as generalized in self.get_defined_node_attributes()
##         # ----------------------------------------------
##                
##         # THIS METHOD COULD BE SPEED UP IF ALL FOR EACH ATTRIBUTE WE EXECUTE ALL NODES AT THE SAME TIME, BUT IT MAY USE MORE MEMORY...
##         attributes_dict = {}
##         for current_attribute in attributes:
##             attributes_dict[current_attribute] = self.dbAccess.get_user_entity_attributes( unification_protocol_name = self.unification_protocol_name,
##                                                                                            listUserEntityID = user_entity_id_list,
##                                                                                            attribute_identifier = current_attribute )

        for current_node in user_entity_id_list:

            uEobj = self.get_user_entity(user_entity_id=current_node)

            # TEMP COMMENTED JAVI
            #eE_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = uEobj.get_externalEntitiesIds_set(), 
            #                                                    attribute_list = attributes )

            new_values = [current_node]

            if include_level_info:
                new_values.append(user_entity_set.getNodeLevel(current_node))

            if include_degree_info:
                new_values.append(user_entity_set.get_node_degree(current_node))

            if include_tags_info:
                tags = user_entity_set.get_user_entity_tags( user_entity_id  = current_node )
                if len(tags)>0:
                    new_values.append( ", ".join(tags) )
                else:
                    new_values.append("-")

            for current_tag_linkage_degree in include_tags_linkage_degree_info:
                new_values.append(user_entity_set.getUserEntityTagLinkageDegree( userEntityID = current_node, tag = current_tag_linkage_degree))

                
            for current_attribute in attributes:
                #new_values.extend( self.get_defined_node_attributes(user_entity_set.id, current_node, current_attribute, output_1_value_per_attribute, substitute_node_attribute_if_not_exists = substitute_node_attribute_if_not_exists, return_set = True) )
                attribute_values = self.get_defined_node_attributes(user_entity_set.id, current_node, current_attribute, output_1_value_per_attribute, substitute_node_attribute_if_not_exists = substitute_node_attribute_if_not_exists, return_set = True)
                temp_str = ", ".join([ x.replace("\n"," ") for x in attribute_values])  #Remove new line in attributes
                if temp_str != "": #len(attribute_values)>0:
                    new_values.append( temp_str )
                else:
                    new_values.append("-")
                
                # TEMP COMMENTED JAVI
                #current_attribute_values = [ str(y.value) for current_eE in uEobj.get_externalEntitiesIds_set() 
                #                             for y in eE_dict[current_eE].get_attribute(attribute_identifier=current_attribute) ]
                
                #current_attribute_values_dict = self.dbAccess.get_user_entity_attributes( unification_protocol_name = self.unification_protocol_name,
                #                                                                          listUserEntityID = [current_node],
                #                                                                          attribute_identifier = current_attribute )

            if output_format == "tabulated":
                out_method( output_utilities.get_tabulated_table(values = [new_values]) )
            elif output_format == "xml":
                out_method(output_utilities.append_html_table_values( rowIDs = [current_node], values = [new_values] ))
            else:
                raise ValueError("output_format is not valid. Valid output formats are ['xml', 'tabulated']\n")
        
        if output_format == "xml":
            out_method(output_utilities.get_html_table_foot())

        OutBianaInterface.send_end_process_message()


    def output_user_entity_set_details(self, user_entity_set_id, out_method = None, attributes=[], include_level_info=True, include_degree_info=True, level=None, only_selected=False, output_format="tabulated", include_tags_info = True, include_tags_linkage_degree_info=[], substitute_node_attribute_if_not_exists=False, output_1_value_per_attribute=True, include_command_in_rows=False, output_only_native_values=False):
        """
        Outputs given user entity set. 
        ------
        user_entity_set_id: Identifier of the user entity set, user entities of which will be outputted
        out_method: method to which output information will be directed
        attributes: list of attribute names corresponding to the attributes to be outputed
        include_level_info: Boolean. Defines inclussion of network level information in the output. 
        include_degree_info: Boolean. Defines inclussion of node degree information in the output.
        inlcude_tags_info: Boolean. Defines inclussion of node tags information in the output.
        include_tags_linkage_degree_info: Calculate linkage degree considering only tagged nodes where tags given as a list.
        output_1_value_per_attribute: Boolean. Defines wether 1 or multiple values are outputed per each attribute.
        output_format: format for the output. Either "tabulated" or "xml".
        only_selected: Boolean to decide whether to output only selected nodes or all nodes. "level" parameter is used only it is False.
	level: if different than None, information for the nodes given in that level are outputed.
        include_command_in_rows: Include the command to output individual user entity information at each row.
	output_only_native_values: If true, for each attribute, shows only values coming from the database providing that attribute by default.
        """

        OutBianaInterface.send_process_message("Getting information from database...")

        user_entity_set = self.dictUserEntitySet[user_entity_set_id]

        if only_selected:
            nodes = user_entity_set.getSelectedUserEntities()
        else:
            nodes = user_entity_set.get_user_entity_ids(level=level)
        
	if output_only_native_values:
	    f_method = self.output_user_entity_details_with_only_defaults
	else:
	    f_method = self.output_user_entity_details

	f_method(user_entity_set = user_entity_set,
                 user_entity_id_list=nodes,
                 out_method = out_method, 
                 attributes = attributes, 
                 include_level_info = include_level_info, 
                 include_degree_info=include_degree_info,
                 output_format=output_format,
                 include_tags_info = include_tags_info,
                 include_tags_linkage_degree_info = include_tags_linkage_degree_info,
                 substitute_node_attribute_if_not_exists = substitute_node_attribute_if_not_exists, 
                 output_1_value_per_attribute = output_1_value_per_attribute,
                 include_command_in_rows = include_command_in_rows)

        OutBianaInterface.send_end_process_message()


    def output_external_entity_details(self, out_method=None, attributes=[], user_entity_id_list=[]):
        """
        Outputs a summary of the external entities given as a list
        ------
        user_entity_id_list: list of user entities 
        out_method: method to which output information will be directed
        attributes: list of attribute names corresponding to the attributes to be outputed

        Called by BianaGUI/showTableDialog.java via output_user_entity_details method
        """

        OutBianaInterface.send_process_message("Getting information from database...")

        if out_method is None:
            out_method = self.outmethod

        if not isinstance(attributes,list):
            attributes = attributes.split(",")
            
        columns = ["External Entity ID","User Entity ID","Type","Default attribute","External Database"]
        columns.extend(attributes)
        values = []

        list_eE_to_search_dict = dict([ (eE_id,x) for x in user_entity_id_list for eE_id in self.get_user_entity(x).get_externalEntitiesIds_set() ])

        if (len( list_eE_to_search_dict ) > 0 ):
            eE_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = list_eE_to_search_dict.keys(), attribute_list = attributes )
        else:
            raise ValueError("LIST CANNOT BE EMPTY FOR OUTPUT EXTERNAL ENTITY DETAILS...")

        default_ids_dict = dict([ (x,"-") for x in eE_dict.keys() ])
        default_ids_dict.update(self.dbAccess.get_default_external_entity_ids( externalEntityIDsList=eE_dict.keys() ))
            
        for current_eE in eE_dict.values():
            current_values = [ current_eE.get_id() ]
            current_values.append(list_eE_to_search_dict[current_eE.get_id()])
            current_values.append( current_eE.get_type() )
            current_values.append( default_ids_dict[current_eE.get_id()] )

            current_values.append("%s" %self.dbAccess.get_external_database( database_id = current_eE.get_source_database() ) )

            for current_attribute in attributes:
                attribute_values = [ str(y.value).replace("\n"," ") for y in current_eE.get_attribute(attribute_identifier=current_attribute) ]
                temp_str = ",".join(attribute_values)
                #if len(attribute_values)>0:
                if temp_str != "":
                    current_values.append( temp_str )
                else:
                    current_values.append("-")

            values.append(current_values)

        out_method(output_utilities.get_html_table(columns,values,attributes=[("title","External Entity Details")]))

        OutBianaInterface.send_end_process_message()

    def output_user_entity_details_with_only_defaults(self, user_entity_set, user_entity_id_list, out_method = None, attributes = [], include_level_info = True, include_degree_info=True, include_tags_info=True, include_tags_linkage_degree_info=[], substitute_node_attribute_if_not_exists=False, output_1_value_per_attribute=True, output_format="tabulated", include_command_in_rows=False):
        """
        Outputs given user entity (for each attribute, shows only values coming from the database providing that attribute by default). Called by output_user_entity_set_details.
        ------
        user_entity_set: Instace of the user entity set, user entities of which will be outputted
        user_entity_id_list: list of user entities 
        out_method: method to which output information will be directed
        attributes: list of attribute names corresponding to the attributes to be outputed
        include_level_info: Boolean. Defines inclussion of network level information in the output. 
        include_degree_info: Boolean. Defines inclussion of node degree information in the output.
        inlcude_tags_info: Boolean. Defines inclussion of node tags information in the output.
        include_tags_linkage_degree_info: Calculate linkage degree considering only tagged nodes where tags given as a list
        output_1_value_per_attribute: Boolean. Defines wether 1 or multiple values are outputed per each attribute
        output_format: format for the output. Either "tabulated" or "xml"
        include_command_in_rows: Include the command to output individual user entity information at each row
        """

        OutBianaInterface.send_process_message("Getting information from database...")

        # Define columns
        columns = ["User Entity ID"]
        tld = {}
        if include_level_info:
            columns.append("Level")
        if include_degree_info:
            columns.append("Degree")
        if include_tags_info:
            columns.append("Tags")
        for current_tag_linkage_degree in include_tags_linkage_degree_info:
            columns.append("%s linkage degree" %current_tag_linkage_degree)
            tld[current_tag_linkage_degree] = user_entity_set.getTagLinkageDegree(tag = current_tag_linkage_degree)
        
        columns.extend([ str(x) for x in attributes])

        if include_command_in_rows:
            command = "output_external_entity_details(user_entity_id_list=[?]"
            if len(attributes)>0:
                command += ", attributes=[\'%s\'])" %("\',\'".join(attributes))
            else:
                command += ")"
        else:
            command = "" 

        if out_method is None:
            out_method = self.outmethod

        if output_format == "xml":
            out_method(output_utilities.get_html_table_header( columns = columns, attributes = [("id", "user_entity_set_nodes"), ("command",command),("title","User Entity Set Details"),("session",self.sessionID)] ) )

        elif output_format == "tabulated":                                            
            out_method( output_utilities.get_tabulated_table(values = [columns]) )   

	#! commented (emre)
	#list_eE_to_search_dict = dict([ (eE_id,x) for x in user_entity_id_list for eE_id in self.get_user_entity(x).get_externalEntitiesIds_set() ])

	#if (len( list_eE_to_search_dict ) > 0 ):
	#    eE_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = list_eE_to_search_dict.keys(), attribute_list = attributes )
	#else:
	#    raise ValueError("LIST CANNOT BE EMPTY FOR OUTPUT EXTERNAL ENTITY DETAILS...")

	#default_ids_dict = dict([ (x,"-") for x in eE_dict.keys() ])
	#default_ids_dict.update(self.dbAccess.get_default_external_entity_ids( externalEntityIDsList=eE_dict.keys() ))

        for current_node in user_entity_id_list:

            uEobj = self.get_user_entity(user_entity_id=current_node)

            # TEMP COMMENTED JAVI #! uncommented and default dict moved here (emre)
            eE_dict = self.dbAccess.get_external_entities_dict( externalEntityIdsList = uEobj.get_externalEntitiesIds_set(), 
                                                                attribute_list = attributes )
	    default_ids_dict = dict([ (x,"-") for x in eE_dict.keys() ])
	    default_ids_dict.update(self.dbAccess.get_default_external_entity_ids( externalEntityIDsList=eE_dict.keys() ))

            new_values = [current_node]

            if include_level_info:
                new_values.append(user_entity_set.getNodeLevel(current_node))

            if include_degree_info:
                new_values.append(user_entity_set.get_node_degree(current_node))

            if include_tags_info:
                tags = user_entity_set.get_user_entity_tags( user_entity_id  = current_node )
                if len(tags)>0:
                    new_values.append( ", ".join(tags) )
                else:
                    new_values.append("-")

            for current_tag_linkage_degree in include_tags_linkage_degree_info:
                new_values.append(user_entity_set.getUserEntityTagLinkageDegree( userEntityID = current_node, tag = current_tag_linkage_degree))

                
            for current_attribute in attributes:
		attribute_values = []
		for current_eE in eE_dict.values():
		    #if current_node != list_eE_to_search_dict[current_eE.get_id()]: #! commented (emre)
		    #	continue
		    default_id = default_ids_dict[current_eE.get_id()]
		    if default_id.startswith(current_attribute):
			attribute_values.append(default_id[default_id.find(":")+1:])
			#attribute_values.extend([ str(y.value) for y in current_eE.get_attribute(attribute_identifier=current_attribute) ])
                temp_str = ", ".join(attribute_values)
                if temp_str != "": #len(attribute_values)>0:
                    new_values.append( temp_str )
                else:
                    new_values.append("-")
                
                # TEMP COMMENTED JAVI
                #current_attribute_values = [ str(y.value) for current_eE in uEobj.get_externalEntitiesIds_set() 
                #                             for y in eE_dict[current_eE].get_attribute(attribute_identifier=current_attribute) ]
                
                #current_attribute_values_dict = self.dbAccess.get_user_entity_attributes( unification_protocol_name = self.unification_protocol_name,
                #                                                                          listUserEntityID = [current_node],
                #                                                                          attribute_identifier = current_attribute )

            if output_format == "tabulated":
                out_method( output_utilities.get_tabulated_table(values = [new_values]) )
            elif output_format == "xml":
                out_method(output_utilities.append_html_table_values( rowIDs = [current_node], values = [new_values] ))
            else:
                raise ValueError("output_format is not valid. Valid output formats are ['xml', 'tabulated']\n")
        
        if output_format == "xml":
            out_method(output_utilities.get_html_table_foot())

        OutBianaInterface.send_end_process_message()



    def __str__(self):
        return "BIANA Session with %d User Entity Sets: \n%s\n" % ( len(self.dictUserEntitySet), "\n".join([ "%s" % user_entity_set for user_entity_set_id, user_entity_set in self.dictUserEntitySet.iteritems() ]) )


    def deliver_error_message(self, message):
        OutBianaInterface.send_error_notification(error="BIANA Error!", message=message)


    def output_user_entity_set_sequences_in_fasta(self, user_entity_set_id, only_selected = False, out_method = None, type="proteinsequence", attributes = [], include_tags_info=True, one_line_per_sequence=False, output_1_value_per_attribute=False, output_only_native_values=False):
        """
        Outputs sequences associated with given user entity set in fasta format. 
        ------
        user_entity_set_id: Id of user entity set containing user entities whose sequence will be outputted
        only_selected: Includes all user entites if False, otherwise only selected ones
        out_method: method to which output information will be directed
        type: either "proteinsequence" or "nucleotidesequence"
        attributes: list of attribute names corresponding to the attributes to be outputed as sequence header
        inlcude_tags_info: Boolean. Defines inclussion of node tags information in the output.

        Called by BianaGUI/ViewUserEntityDetails.java 
        """

	if output_only_native_values is True or output_1_value_per_attribute is True:
	    OutBianaInterface.send_error_notification(message="Output Error!" ,error="Displaying only single/native value feature is not available in FASTA formatted output. Use All")
	    return

        user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        
        if only_selected: # get selected user entities of this set
            user_entity_id_list = user_entity_set.getSelectedUserEntities()
        else: # get all user entities of this set
            user_entity_id_list = user_entity_set.get_user_entity_ids()

        if out_method is None:
            out_method = self.outmethod
       
        columns = ["externalentityid","userentityid","type","defaultattribute","externaldatabase"]
        if include_tags_info:
            columns.append("tags")

        new_attributes = []
        for x in attributes:
            if x.lower() == "proteinsequence":
                pass
            elif x.lower() == "nucleotidesequence":
                pass
            else:
                columns.append(x.lower())
                new_attributes.append(x)

        if type == "proteinsequence" or "nucleotidesequence":
            new_attributes.append(type)
        else:
            raise Exception("Unrecognized sequence type: %s !" % type)

        values = []
        sequences = []

        list_eE_to_search_dict = dict([ (eE_id,x) for x in user_entity_id_list for eE_id in self.get_user_entity(x).get_externalEntitiesIds_set() ])

	default_ids_dict = dict([ (x,"-") for x in list_eE_to_search_dict.keys() ])
	default_ids_dict.update(self.dbAccess.get_default_external_entity_ids( externalEntityIDsList=list_eE_to_search_dict.keys() ))

        for current_eE_id, uE_id in list_eE_to_search_dict.iteritems():

            current_eE = self.dbAccess.get_external_entities_dict( externalEntityIdsList = [current_eE_id], attribute_list = new_attributes )[current_eE_id]

            #default_ids_dict = self.dbAccess.get_default_external_entity_ids( externalEntityIDsList=[current_eE.get_id()] )

            current_values = [ current_eE.get_id() ]
            current_values.append(uE_id)
            current_values.append( current_eE.get_type() )
            default_id = default_ids_dict[current_eE.get_id()]
            current_values.append(default_id)
            current_values.append("%s" %self.dbAccess.get_external_database( database_id = current_eE.get_source_database() ) )

            skip = False

            if include_tags_info:
                tags = user_entity_set.get_user_entity_tags( user_entity_id  = uE_id )
                if len(tags)>0:
                    current_values.append( ",".join(tags) )
                else:
                    current_values.append("-")

            for current_attribute in new_attributes:
                attribute_values = [ str(y.value).replace("\n"," ") for y in current_eE.get_attribute(attribute_identifier=current_attribute) ]

                if current_attribute == "proteinsequence" or current_attribute == "nucleotidesequence":
                    if len(attribute_values)>0:
                        sequences.append(",".join(attribute_values))
                    else:
                        skip = True
                    #else:
                    #    sequences.append("-")
                else:
                    temp_str = ",".join(attribute_values)

                    if temp_str != "":
                        current_values.append( temp_str )
                    else:
                        current_values.append("-")

            if skip is False:
                values.append(current_values)
            
        import biana.utilities.FastaWriter as FastaWriter
        fasta_writer = FastaWriter.FastaWriter(out_method = out_method, one_line_per_sequence=one_line_per_sequence)

        for current_values, sequence in zip(values, sequences):

            sequence_header = "|".join( map((lambda x, y: x+'|'+str(y)), columns, current_values) )
            [ fasta_writer.output_sequence(sequence_header = sequence_header, sequence = seq) for seq in sequence.split(",")]
            
        return
 

    # BIANA REPORT METHODS #

    def initReport(self, reportOutputMethod, reportTitle, reportFormat, url):
        """
        Inits the Biana report where a description of the experiment and results will be written

        "reportTitle" is the title that will be written on top of the report (default is 'BIANA report')
        
        "reportOutputMethod" is the write method where the report will be written (e.g. file_object.write)
              --> if no streamOutputMethod is given, report is written to stdout

        "reportFormat" is the formatting that will be applied to the report
             - 'html': report in html
             - 'txt': report in raw text

        "url" is the URL where the interface to BIANA is placed
        """
        self.report = BianaReport(title=reportTitle, streamOutputMethod=reportOutputMethod, 
                                  format=reportFormat, url=url)

        return

    def closeReport(self, reportOutputMethod, reportFormat):
        """
        Closes the Biana report. This method is in charge of writing the closing HTML tags of the report

        "reportOutputMethod" is the write method where the report will be written (e.g. file_object.write)
              --> if no streamOutputMethod is given, report is written to stdout

        "reportFormat" is the formatting that will be applied to the report
             - 'html': report in html
             - 'txt': report in raw text

        Attention!!! This method does not close the file object associated to the report. The user is responsible
                     for closing the file (if needed).
        """

        if self.report is None:

            raise ValueError("This BIANA session does not have an associated report. You should init one using initReport()")
        self.report.closeReport()

        return

    def addTextToReport(self, theText=None, type="regular" ):
        """
        This method adds text to the report.

        Depending on the 'type' of text, it will be written one way or another

          - 'regular' prints standard text
          - 'title' prints a section title
          - 

        """

        if self.report is None:

            raise ValueError("This BIANA session does not have an associated report. You should init one using initReport()")
        
        self.report.addText(theText=theText, type=type )


    def addResultToReport(self, resultFileName=None, resultFilePath=None, associatedText=None, bulletedList=[] ):
        """
        This method adds information to the resport about a result file.

        'resultFileName' is the name you want to give to this result file

        'resultFilePath' is the path (or URL) of the result file (e.g. http://localhost/this_file.html

        'associatedText' is the text shown next to the results file name (i.e description of the file)

        'bulletedList' is a list of strings that will be shown under the results name as a bulleted list

        """
 
        if self.report is None:

            raise ValueError("This BIANA session does not have an associated report. You should init one using initReport()")
        self.report.addResult(resultFileName=resultFileName, 
                              resultFilePath=resultFilePath, 
                              associatedText=associatedText, 
                              bulletedList=bulletedList )

    # OTHERS #
    def get_new_user_entity_set_from_user_entities(self, user_entity_id_list, new_user_entity_set_id=None, external_entity_attribute_restriction_list=[]):
        """
        Creates a new user entity set from given user entities.
        ------
        user_entity_id_list: user entities to be included new user entity set
        user_entity_set_id: identifier of user entity set to be created
        external_entity_attribute_restriction_list: Restrictions to not to include user entities with external entities having those restriction values. List of tuples provided by the user containing external entity attribute restrictons to be applied. They will be used always in the set (network creation, etc)

        # Used by BIANAGUI
        """

        if new_user_entity_set_id is None:
            new_user_entity_set_id = self._get_next_uEs_id()
            
        user_entity_set_new = UserEntitySet.UserEntitySet( id = new_user_entity_set_id, 
                                                           setIdUserEntity = user_entity_id_list )
        
        user_entity_set_new.addRestriction("attribute_restrictions", external_entity_attribute_restriction_list)
        
        self.dictUserEntitySet[user_entity_set_new.id] = user_entity_set_new

        return user_entity_set_new

    def get_sub_user_entity_set(self, user_entity_set_id, include_relations=False, include_restrictions=True, new_user_entity_set_id=None):
        """
        Creates a new user entity set using selected user entities
        ------
        user_entity_set_id: identifier of user entity set from whose user entities a new one will be created
        new_user_entity_set_id: identifier of user entity set to be created
        include_relations: flag to whether or not include relations in the new user entity set
        """

        if new_user_entity_set_id is None:
            new_user_entity_set_id = self._get_next_uEs_id()

        user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        
        #self.idLastUserEntitySet += 1
        selected_user_entities = user_entity_set.getSelectedUserEntities()
        selected_relations = user_entity_set.getRelationsOfSelectedUserEntities()
        if include_relations:
            user_entity_set_new = UserEntitySet.UserEntitySet( id = new_user_entity_set_id, 
                                                               setIdUserEntity = selected_user_entities,
                                                               listRelations = selected_relations)
        else:
            user_entity_set_new = UserEntitySet.UserEntitySet( id = new_user_entity_set_id, 
                                                               setIdUserEntity = selected_user_entities)

	if include_restrictions:
	    user_entity_set_new.addRestriction("attribute_restrictions", user_entity_set.getRestrictions("attribute_restrictions"))

        self.dictUserEntitySet[user_entity_set_new.id] = user_entity_set_new

        #self.outmethod(self._get_xml(inner_content="<new_user_entity_set id=\"%s\"/>" %(new_user_entity_set_id)))
        #self.outmethod(self._get_xml(inner_content=user_entity_set_new._get_xml(inner_content="all")))
        self._send_complete_user_entity_set_info( user_entity_set_new )

        return user_entity_set_new


    def unselect_user_entities_from_user_entity_set(self, user_entity_set_id):
        """
        Clears selection of user entities in a given user entity set 
        ------
        user_entity_set_id: identifier of user entity set 
        """
        user_entity_set = self.get_user_entity_set(user_entity_set_id)
        user_entity_set.clear_user_entity_selection()
        self.outmethod(self._get_xml(user_entity_set._get_xml(inner_content="<clear_user_entity_selection/>")))
        return

    def unselect_user_entity_relations_from_user_entity_set(self, user_entity_set_id):
        """
        Clears selection of relations in a given user entity set 
        ------
        user_entity_set_id: identifier of user entity set 
        """
        user_entity_set = self.get_user_entity_set(user_entity_set_id)
        user_entity_set.clear_user_entity_relation_selection()
        self.outmethod(self._get_xml(user_entity_set._get_xml(inner_content="<clear_user_entity_relation_selection/>")))
        return

    def select_all_user_entities(self, user_entity_set_id):
        """
        Selects all user entities in a given user entity set 
        ------
        user_entity_set_id: identifier of user entity set 
        """
        user_entity_set = self.get_user_entity_set(user_entity_set_id)
        user_entity_id_list = user_entity_set.get_user_entity_ids()
        self.select_user_entities_from_user_entity_set( user_entity_set_id = user_entity_set_id,
                                                  user_entity_id_list = user_entity_id_list )
        return


    def select_all_user_entity_relations(self, user_entity_set_id):
        """
        Selects all relations in a given user entity set 
        ------
        user_entity_set_id: identifier of user entity set 
        """
        user_entity_set = self.get_user_entity_set(user_entity_set_id)
        #user_entity_relation_id_list = user_entity_set.getRelations()  
        self.select_user_entity_relations_from_user_entity_set( user_entity_set_id = user_entity_set_id,
                                                                #user_entity_relation_id_list = user_entity_relation_id_list )
                                                                external_entity_relation_ids_list = user_entity_set.get_external_entity_relation_ids() )

        #self.select_user_entity_relations_from_user_entity_set( user_entity_set_id = user_entity_set_id,
        #                                                        user_entity_relation_id_list = user_entity_relation_id_list )
        return


    def select_user_entities_from_user_entity_set(self, user_entity_set_id, user_entity_id_list, clear_previous_selection = False):
        """
        Select user entities with given id in a given user entity set 
        ------
        user_entity_set_id: identifier of user entity set 
        user_entity_id_list: list of user entity identifiers
        clear_previous_selection: Boolean to clear or not previous selection
        """
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]

        xml = []

        if clear_previous_selection:
            self.unselect_user_entities_from_user_entity_set(user_entity_set_id)


        user_entity_set.select_user_entities(user_entity_id_list)

        # To send it one by one
        #xml.extend([ "<select_user_entity id=\"%s\"/>" %x for x in user_entity_id_list ])
        xml.append( "<select_user_entity id=\"%s\"/>" %",".join(map(str,user_entity_id_list)) )

        self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content="".join(xml))))
        return

    def select_user_entity_relations_from_user_entity_set(self, user_entity_set_id, user_entity_relation_id_list=[], clear_previous_selection = False, external_entity_relation_ids_list = []):
        """
        Select user entities with given id in a given user entity set 
        ------
        user_entity_set_id: identifier of user entity set 
        user_entity_relation_id_list: list of relation identifiers: (userEntityID1, userEntityID2, relation_type)
        clear_previous_selection: Boolean to clear or not previous selection
        """

        user_entity_set = self.dictUserEntitySet[user_entity_set_id]


        if clear_previous_selection:
            self.unselect_user_entity_relations_from_user_entity_set(user_entity_set_id)

        eErID_list = []

        for node1, node2, relation_type in user_entity_relation_id_list:
            relation_type_enum = self.eEr_types_enum.get_letter(relation_type)
            for current_eErID in user_entity_set.get_external_entity_relation_ids(node1,node2):
                if self.eEr_types_dict[current_eErID]==relation_type_enum:
                    eErID_list.append(current_eErID)

        eErID_list.extend(external_entity_relation_ids_list)

        user_entity_set.selectUserEntityRelations(eErID_list)

        #user_entity_set.selectUserEntityRelations(user_entity_relation_id_list)

        # To send it one by one
        #xml = [ "<select_user_entity_relation id=\"%s\"/>" %",".join([",".join(map(str,x)) for x in user_entity_relation_id_list]) ]
        xml = []
        for current_eErID in eErID_list:
            participants = user_entity_set.get_external_entity_relation_participants(current_eErID)
            for x in xrange(len(participants)):
                for y in xrange(x):
                    xml.append( "<select_user_entity_relation id=\"%s,%s,%s\"/>" %(participants[x],participants[y],self.eEr_types_enum.get(self.eEr_types_dict[current_eErID])) )
        
        self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content="".join(xml))))

        return


    def tag_selected_user_entities( self, user_entity_set_id, tag):
        """
        Tags selected user entities in a given user entity set with the given tag
        ------
        user_entity_set_id: identifier of user entity set 
        tag: tag string
        """
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        
        #print user_entity_set.hasTag( tag = tag )
        if user_entity_set.hasTag( tag = tag ) is False:
            self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content="<new_tag tag=\"%s\"/>" %tag)))
            
        user_entity_set.addTagToSelectedUE(tag=tag)
        return

    def tag_selected_user_entity_relations( self, user_entity_set_id, tag):
        """
        Tags selected relations in a given user entity set with the given tag
        ------
        user_entity_set_id: identifier of user entity set 
        tag: tag string
        """
        
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]
        
        #print user_entity_set.hasTag( tag = tag )
        if user_entity_set.relationHasTag( tag = tag ) is False:
            self.outmethod(self._get_xml(inner_content=user_entity_set._get_xml(inner_content="<new_relation_tag tag=\"%s\"/>" %tag)))
            
        user_entity_set.addTagToSelectedUER(tag=tag)
        return

    def select_user_entities_using_attributes(self, user_entity_set_id, identifier_description_list, id_type="embedded", external_entity_attribute_restriction_list=[], negative_external_entity_attribute_restriction_list=[], clear_previous_selection = False ):
        """
        Selects user entities in a given user entity set by given attributes
        ------
        user_entity_set_id: identifier of user entity set 
        identifier_description_list: list of identifiers (or (identifier, id_type) tuples in case id_type is "embedded")
        id_type: type of the identifiers in the file, if "embedded" then file contains (attribute_name, identifier) tuples instead of just identifiers
        external_entity_attribute_restriction_list: Restrictions to not to include user entities with external entities having those restriction values. List of tuples provided by the user containing external entity attribute restrictons to be applied. 
        clear_previous_selection: Boolean to clear or not previous selection
        """
        user_entity_set = self.dictUserEntitySet[user_entity_set_id]

        if clear_previous_selection:
            self.unselect_user_entities_from_user_entity_set(user_entity_set_id)

        # Gets the input to build the user entity set (root or seed user entities)
        dictTypeToName = self._convert_attribute_list_to_attribute_dictionary(identifier_description_list, id_type)

        OutBianaInterface.send_process_message("Getting information from database...")

        ## for each different type of id_type's fetch userEntity ids associated with given externalEntity ids
        for identifierType, setIdentifierName in dictTypeToName.iteritems():
            
            # javi added:
            if( identifierType.lower()=="userentityid" ):
                user_entity_id_list = [int(identifierName)]

            field_values = []
            
            for identifierName in setIdentifierName:
                field_values.append(("value",identifierName))

            user_entity_id_list = self.dbAccess.get_list_user_entities_IDs_by_attribute(unification_protocol_name = self.unification_protocol_name,
                                                                                     attribute_identifier = identifierType,
                                                                                     field_values = field_values,
                                                                                     attribute_restrictions = external_entity_attribute_restriction_list,
                                                                                     negative_attribute_restrictions = negative_external_entity_attribute_restriction_list,
                                                                                     restrict_to_user_entity_ids_list = user_entity_set.get_user_entity_ids() )

            
            self.select_user_entities_from_user_entity_set( user_entity_set_id = user_entity_set_id, user_entity_id_list = user_entity_id_list, clear_previous_selection = False)


        OutBianaInterface.send_end_process_message()

        return


    def select_user_entity_relations_using_attributes(self, user_entity_set_id, identifier_description_list, id_type="embedded", relation_attribute_restriction_list=[], relation_type_list=[], clear_previous_selection = False ):
        """
        Selects relations in a given user entity set by given attributes
        ------
        user_entity_set_id: identifier of user entity set 
        identifier_description_list: list of identifiers (or (identifier, id_type) tuples in case id_type is "embedded")
        id_type: type of the identifiers in the file, if "embedded" then file contains (attribute_name, identifier) tuples instead of just identifiers
        realation_attribute_restrictions: Restrictions to not to include relations with external entity relations having those restriction values. List of tuples provided by the user containing external entity attribute restrictons to be applied. 
        relation_type_list: types of relations to be included in selection

        # Not called by any method - for user scripts
        """

        user_entity_set = self.dictUserEntitySet[user_entity_set_id]

        if clear_previous_selection:
            self.unselect_user_entity_relations_from_user_entity_set(user_entity_set_id)

        # Gets the input to build the user entity set (root or seed user entities)
        dictTypeToName = self._convert_attribute_list_to_attribute_dictionary(identifier_description_list, id_type )

        OutBianaInterface.send_process_message("Getting information from database...")

        ## for each different type of id_type's fetch userRelatinoEntity ids associated with given externalRelationEntity ids
        for identifierType, setIdentifierName in dictTypeToName.iteritems():
            
            # javi added:
            #if( identifierType.lower()=="userentityrelationid" ):
            #    user_entity_relation_id_list = [int(identifierName)]

            field_values = []
            
            for identifierName in setIdentifierName:
                field_values.append(("value",identifierName))

            relation_type_list.extend(user_entity_set.getRestrictions("relation_type_restrictions"))
            #realation_attribute_restrictions.extend(user_entity_set.getRestrictions("relation_attribute_restrictions"))
            dictAttributeToValues = self._convert_attribute_list_to_attribute_dictionary(realation_attribute_restrictions, "embedded")
            dictAttributeToValuesExisting = user_entity_set.getRestrictions("relation_attribute_restrictions")
            for attr, setValue in dictAttributeToValuesExisting:
                dictAttributeToValues.setdefault(attr, set()).union(setValue)

            listTupleIdUserEntity = self.dbAccess.get_user_entity_relations(unification_protocol_name = self.unification_protocol_name,
                                                                            userEntityID_list = user_entity_set.get_user_entity_ids(),
                                                                            attribute_restrictions = user_entity_set.getRestrictions("attribute_restrictions"),
                                                                            negative_attribute_restrictions = user_entity_set.getRestrictions("negative_attribute_restrictions"),
                                                                            listRelationType = relation_type_list,
                                                                            dictRelationAttributeRestriction = dictAttributeToValues,
                                                                            use_self_relations = user_entity_set.getRestrictions("use_self_relations"),
                                                                            limit_to_userEntityID_list = True)

            self.select_user_entity_relations_from_user_entity_set( user_entity_set_id = user_entity_set_id, user_entity_relation_id_list = listTupleIdUserEntity, clear_previous_selection = False)

        OutBianaInterface.send_end_process_message()

        return
    
    
    def cluster_user_entities(self, user_entity_set_id, attribute):
        """
        Clusters the user Entities according to an attribute

        # Not called by any other method
        """

        # TO CHECK

        eE_uE_dict = {}

        #import networkx
        import biana.ext.networkx as networkx
        g = networkx.Graph()
        
        #for currentUE in self.dictUserEntity.values():
        #    eE_set = currentUE.get_externalEntitiesIds_set()
        #    g.add_node(currentUE.id)
        #    for current_eE_id in eE_set:
        #        eE_uE_dict[current_eE_id] = currentUE.id

        uEs = self.get_user_entity_set(user_entity_set_id)

        user_entity_id_list = uEs.get_user_entity_ids()
        for currentUEid in user_entity_id_list: #getListUserEntityId():
            currentUE = self.get_user_entity(currentUEid)
            eE_set = currentUE.get_externalEntitiesIds_set()
            g.add_node(currentUE.id)
            for current_eE_id in eE_set:
                eE_uE_dict[current_eE_id] = currentUEid

        self.userEntityCluster = {}

        equivalent_eE = self.dbAccess.get_equivalent_external_entities_from_list( externalEntitiesList = eE_uE_dict.keys(), attribute=attribute )

        for current_equivalent_list in equivalent_eE:
            for x in xrange(len(current_equivalent_list)-1):
                g.add_edge(eE_uE_dict[current_equivalent_list[x]],eE_uE_dict[current_equivalent_list[x+1]])

        uEs.clusterUserEntities(clusters=networkx.connected_components(g))

        return
    
    def create_new_user_entity_set_and_network_from_sif_file(self, sif_file_name, new_user_entity_set_id=None): 
        """
        create userEntity objects and their network from given sif file
        ------
	sif_file_name: name of the sif file
        new_user_entity_set_id: identifier of the set provided by the user
        """
        
        OutBianaInterface.send_process_message("Creating new user entity set and its network.\nProcessing information...")

        try:
            # Check for missing parameters and set it to default
            if new_user_entity_set_id is None:
                new_user_entity_set_id = self._get_next_uEs_id()

            user_entity_set = UserEntitySet.UserEntitySet(new_user_entity_set_id )

            # Fill userProvidedExtIdsLowerDict:
            for k in user_entity_set.userProvidedExternalIdsDict:
                user_entity_set.userProvidedExtIdsLowerDict[k] = {}

                for extId in user_entity_set.userProvidedExternalIdsDict[k]:
                    user_entity_set.userProvidedExtIdsLowerDict[k][str(extId).lower()] = extId

	    f = open(sif_file_name)
	    user_entity_id = 0
	    externalEntityRelationID = 0
	    ## create new userEntity objects (if not created before) via reading from sif file
	    for line in f.readlines():
		words = line[:-1].split(" pp ")
		if len(words) == 1: 
		    user_entity_id = words[0] #-= 1
		    type = "protein"
                    user_entity_set.addUserEntityId(idUserEntity = user_entity_id, level = 0)
                    self.uE_types_dict[user_entity_id] = self.uE_types_enum.get_letter(type)
		elif len(words) == 2:
		    idUserEntity1, idUserEntity2 = words[0], words[1] 
		    externalEntityRelationID = "pp" #-= 1
		    relation_type = "interaction"
		    partner_type = "protein" 
		    if not user_entity_set.has_user_entity(idUserEntity1):
			user_entity_set.addUserEntityId(idUserEntity = idUserEntity1, level = 0)
			self.uE_types_dict.setdefault(idUserEntity1,self.uE_types_enum.get_letter(partner_type))
		    if not user_entity_set.has_user_entity(idUserEntity2):
			user_entity_set.addUserEntityId(idUserEntity = idUserEntity2, level = 0)
			self.uE_types_dict.setdefault(idUserEntity2,self.uE_types_enum.get_letter(partner_type))
		    # addUserEntityRelation adds nodes so above should be before a call to it
		    user_entity_set.addUserEntityRelation(idUserEntity1 = idUserEntity1,
							     idUserEntity2 = idUserEntity2, 
							     externalEntityRelationID = externalEntityRelationID)
		    self.eEr_types_dict[externalEntityRelationID] = self.eEr_types_enum.get_letter(relation_type)
		else:
		    raise Exception("SIF file format error: %s", line)
        except:
            OutBianaInterface.send_error_notification( message = "New set not created. BIANA ERROR:", error = traceback.format_exc() )
            OutBianaInterface.send_end_process_message()
            return

	user_entity_set.setIsNetworkCreated()

        OutBianaInterface.send_end_process_message()
        
        if user_entity_set.getSize()==0:
            OutBianaInterface.send_error_notification(message="New set not created", error = "Any element has been found with given identifiers")
            return

        self.dictUserEntitySet[user_entity_set.id] = user_entity_set

        OutBianaInterface.send_process_message("Sending data...")
	self.outmethod(self._get_xml_header())
        self._send_complete_user_entity_set_info( user_entity_set = user_entity_set )
	self.outmethod(self._get_xml_foot())
        OutBianaInterface.send_end_process_message()

        return user_entity_set



