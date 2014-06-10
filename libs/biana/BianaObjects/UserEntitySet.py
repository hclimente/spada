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

import sys
from biana.utilities import graph_utilities

import time
#import networkx
import biana.ext.biana_networkx as networkx
import copy


class UserEntitySet(object):
        """
        Class that represents a set of user entities.
        """

        def __init__(self, id, setIdUserEntity=None, listRelations=None, listLevelSetIdUserEntity=None):
		"""
		"""
                
		self.id = id

		self.current_level = 0

                self.is_network_created = False # can only be set create_network command of BianaSessionManager showing whether this command has been called or not

		# Tags related
		self.uE_tags = {}
		self.tag_uE = {}
		
		# Relation Tags related
		self.uER_tags = {}
		self.tag_uER = {}

                # For sending only interactions fetched with the last expansion to Cytoscape
		self._printed_groups = set()

		# Dictionary to store relation groups (relations taken as groups, not as edges)
		# As key stores the external entity relation ID and as a value a set of user entities ids
		self.relation_groups = {}
		self.uE_groups = {} # Stores the groups to which a user entity belongs
		self.groups_hierarchy = {}
		self.group_types = {}

		# List that contains the user entities ids at different levels
		if listLevelSetIdUserEntity is None:
			self.current_level = 0
			if setIdUserEntity is not None:
				self.size = len(setIdUserEntity)
				self.listLevelSetIdUserEntity = [ set(map(None,setIdUserEntity)) ]
			else:
				self.size = 0
				self.listLevelSetIdUserEntity = [ set() ]
		else:
			self.current_level = len(listLevelSetIdUserEntity)-1
			self.size = 0
			for setId in listLevelSetIdUserEntity:
				self.size += len(setId)
			self.listLevelSetIdUserEntity = listLevelSetIdUserEntity
			
		## Network of relations between the elements of this UserEntitySet
		self.network = graph_utilities.create_graph()

		# Add initial nodes and edges
                for current_level in self.listLevelSetIdUserEntity:
			self.network.add_nodes_from(current_level)
		if listRelations is not None:
			self.network.add_edges_from(listRelations)
                        #self.is_network_created = True

		# Selected user entity ids
		self.setIdUserEntitySelected = set()
		
		# Selected user entity relation ids
		self.setIdUserEntityRelationSelected = set()
		
		self.userEntityClusters = None        # Used to cluster user entities according to some criteria. Key: UserEntityID. Cluster: ClusterID
		self.clusteredNetwork = None

		self.nodeLevelsDict = {}              # Dictionary to store the level to which belongs node

		# Output related
		self.userProvidedExternalIdsDict = {} # Dictionary to store external Ids provided by the user.
		                                      # Corresponds to dictTypeToName in BianaSessionManager.create_new_user_entity_set()

		
		self.selectedExternalIdsDict = {}     # Dictionary to store external Ids selected previously by BianaSessionManager.
		                                      # It is updated in BianaSessionManager.outputUserEntityObjectsDetails

		self.userProvidedExtIdsLowerDict = {} # Stores correspondence between lower case and user provided codes in self.userProvidedExternalIdsDict
		self.selectedExtIdsLowerDict = {}     # Stores correspondence between lower case and selected codes in self.selectedExtneralIdsDict


		# Start the nodeLevelsDict in the case it is not empty
		# Start the tags dict in the case it is not empty
		for current_level_x in xrange(len(self.listLevelSetIdUserEntity)):
			for current_node in self.listLevelSetIdUserEntity[current_level_x]:
				self.nodeLevelsDict[current_node] = current_level_x

		self._initialize_restrictions_dict()


		self.eErIds2participants = {}        # Stores the participants in each eErID

		
		# TEMP
		#self.all_shortest_paths = None
		#self.calculated_shortest_paths = {}

		return

	def _initialize_restrictions_dict(self):
		"""
		"""

		self.restrictions_dict = { "negative_attribute_restrictions": [],
					   "attribute_restrictions": [],
					   "externaldatabase_restrictions": [],
					   "relation_type_restrictions": [],
					   "relation_attribute_restrictions": {},
					   "use_self_relations":  False,
					   "expansionlistrelationtype": [],
					   "expansionlevel":  [],
					   "expansionattributeslist": [],
					   "include_relations_last_level": False,
					   "attributenetworklist" : [],
					   "group_relation_type": []}

	def addRestriction(self, restriction_type, restrictions):
            if isinstance(self.restrictions_dict[restriction_type.lower()], list):
		self.restrictions_dict[restriction_type.lower()].extend(restrictions)
            elif isinstance(self.restrictions_dict[restriction_type.lower()], dict):
                # assuming that dictionary contains set as values
                for attr, setValue in restrictions:
                    self.restrictions_dict[restriction_type.lower()].setdefault(attr, set()).union(setValue)
            else:
                print "Warning: trying to add a restriction in non-supported format"


	def setRestriction(self, restriction_type, restriction):
		"""
		"""
		self.restrictions_dict[restriction_type.lower()] = restriction

	
	def getRestrictions(self, restriction_type):
		return self.restrictions_dict[restriction_type.lower()]


        def addTagToSelectedUE(self, tag):
             """
	     Adds a tag to selected user entities in this set
             """

             if not self.tag_uE.has_key(tag):
                  self.tag_uE[tag] = set()

	     if( len(self.setIdUserEntitySelected)==0 ):
		     print "Adding a tag to a set without any user entity selected"

             self.tag_uE[tag].update(self.setIdUserEntitySelected)

             [ self.uE_tags.setdefault(x, set()).add(tag) for x in self.setIdUserEntitySelected ]
	     

        def addTagToSelectedUER(self, tag):
             """
	     Adds a tag to selected user entities in this set
             """

             if not self.tag_uER.has_key(tag):
                  self.tag_uER[tag] = set()

	     if( len(self.setIdUserEntityRelationSelected)==0 ):
		     print "Adding a tag to a set without any user entity relation selected"

             self.tag_uER[tag].update(self.setIdUserEntityRelationSelected)

             [ self.uER_tags.setdefault(x, set()).add(tag) for x in self.setIdUserEntityRelationSelected ]


	def getTagLinkageDegree(self, tag):
		"""
		Calculates the linkage degree of all nodes to "tag"

		returns a dictionary node ids as keys and linkage degree as value
		"""

		if not self.tag_uE.has_key(tag):
		    sys.stderr.write("Tag %s is not defined" %tag)
		    return 

		tld = {}
		for current_node in self.network.nodes():
			value = 0
			for current_neighbor in self.network.neighbors(current_node):
				if current_neighbor!=current_node:
					if self.has_tag( userEntityID = current_neighbor, tag = tag ):
						value += 1
			tld[current_node] = value
					
		return tld

        def get_user_entity_ids_by_tag_linker_degree_cutoff(self, tag, linker_degree_cutoff):
            """
            Get user entities associated with the tag and having a linker degree equal or higher to a given cutoff.
            """
            result_set = set()
	    seed_set = self.get_user_entities_for_tag(tag)
            #print seed_set
            for seed in seed_set:
                for neighbor in self.network.neighbors(seed):
                    #print seed, neighbor, set(self.network.neighbors(neighbor))
                    if neighbor in result_set:
                        continue
                    elif len(set(self.network.neighbors(neighbor)) & seed_set) >= linker_degree_cutoff:
                        result_set.add(neighbor)
            return result_set

        def get_user_entity_ids_by_linker_degree_cutoff(self, linker_degree_cutoff):
            """
            Get user entities having a linker degree equal or higher to a given cutoff.
            """
            result_set = set()
            seed_set = self.listLevelSetIdUserEntity[0]
            #print seed_set
            for seed in seed_set:
                for neighbor in self.network.neighbors(seed):
                    #print seed, neighbor, set(self.network.neighbors(neighbor))
                    if neighbor in result_set:
                        continue
                    elif len(set(self.network.neighbors(neighbor)) & seed_set) >= linker_degree_cutoff:
                        result_set.add(neighbor)
            return result_set

	def getUserEntitySecondaryTagLinkageDegree(self, userEntityID, tag):
		"""
		Calculates a secondary linkage degree. It is defined as the mean of the linkage degree of the neighbors plus the node LD
		"""

                if not self.tag_uE.has_key(tag):
                        sys.stder.write("Tag %s is not defined" %tag)
                        return 0
                value = 0
                for current_neighbor in self.network.neighbors(userEntityID):
                        if current_neighbor!=userEntityID:
				value += self.getUserEntityTagLinkageDegree( userEntityID= current_neighbor, tag= tag)
                return ((float(value)/self.get_node_degree(userEntityID))+self.getUserEntityTagLinkageDegree(userEntityID=userEntityID, tag=tag))


	def getUserEntityTagLinkageDegree(self, userEntityID, tag):
		if not self.tag_uE.has_key(tag):
			sys.stder.write("Tag %s is not defined" %tag)
			return 0
		value = 0
		for current_neighbor in self.network.neighbors(userEntityID):
			if current_neighbor!=userEntityID:
				if self.has_tag( userEntityID = current_neighbor, tag = tag ):
					value += 1
		return value


	def getPathLengthNetwork(self, listUserEntityID, path_length_cutoff=10000):
		"""
		Returns an undirected weighted networkx graph containing only nodes in listUserEntityID

		Nodes are connected if exist a path between them, and weight of the edge consists in the length of shortest path
		
		"""

		return graph_utilities.get_path_network(self.network, listNodes = listUserEntityID, path_length_cutoff=10000 )


	def testfast_getUserEntityTagConnectedMetrics(self, userEntityIDList, tag):
		if not self.tag_uE.has_key(tag):
                        sys.stder.write("Tag %s is not defined" %tag)
                        raise ValueError("Tag not defined")

		return [ self.getUserEntityTagConnectedMetrics(x,tag) for x in userEntityIDList ]
		
		
	def _get_single_source_shortest_path(self, userEntityIDList):
		
		itime = time.time()
		for current_userEntityID in userEntityIDList:
			p = networkx.single_source_shortest_path(self.network, current_userEntityID)
			for current_connected in p:
				self.calculated_shortest_paths[(max(current_userEntityID,current_connected),min(current_userEntityID,current_connected))] = p[current_connected]


		print "All shortest distances calculated in %s seconds" %(time.time()-itime)
		print "Len: ",len(self.calculated_shortest_paths)
		

	def getUserEntityTagConnectedMetrics(self, userEntityID, tag):
		"""
		Test method
		"""

		#if not self.tag_uE.has_key(tag):
		#	sys.stder.write("Tag %s is not defined" %tag)
		#	raise ValueError("Tag not defined")

		#itime = time.time()

		if self.get_node_degree(userEntityID)==0:
			return (0,0,0,10000,0,0,0,0,0,'nan', 'nan', 'nan', 'nan', 'nan')

		num_connected = 0
		sum_path_connected = 0
		hub_path_weight = 0

		value = 0
		min_path = 10000
		max_path = 0

		direct_path_connected = 0 # Number of connected nodes with "tag", but without going though another node with the "tag"
		sum_path_direct_path_connected = 0
		max_path_direct_path_connected = 0
		hub_path_weight_direct_path_connected = 0

		for current_tagged in self.get_user_entities_for_tag( tag = tag ):

			if current_tagged == userEntityID:
				continue
			
			#ptime = time.time()
			path = self.getShortestPathBetween( idUserEntity1 = userEntityID,
							    idUserEntity2 = current_tagged )
			#print "shortest_path_time: ",time.time()-ptime

			#ltime = time.time()

			if path:
				
				num_connected += 1

				sum_path_connected += len(path)-1

				current_hub_value = 1
				
				connected_directly = True

				for current_path_component in path:
					current_hub_value *= 1.0/len(self.network.neighbors(current_path_component))
					if connected_directly == True:
						if current_path_component != userEntityID and current_path_component != current_tagged:
							if self.has_tag( userEntityID = current_path_component, tag = tag ):
								connected_directly = False
				
				if connected_directly:
					direct_path_connected += 1
					sum_path_direct_path_connected += len(path)-1
					hub_path_weight_direct_path_connected += current_hub_value
					max_path_direct_path_connected = max(len(path)-1,max_path_direct_path_connected)

				hub_path_weight += current_hub_value

				min_path = min(min_path,len(path)-1)
				max_path = max(max_path,len(path)-1)

			#print "Loop time: ",time.time()-ltime

		if num_connected>0:
			if direct_path_connected>0:
				return (num_connected, sum_path_connected, hub_path_weight, min_path, max_path, direct_path_connected, sum_path_direct_path_connected, max_path_direct_path_connected, hub_path_weight_direct_path_connected, float(sum_path_connected)/num_connected, hub_path_weight/num_connected, sum_path_direct_path_connected/direct_path_connected, hub_path_weight_direct_path_connected/direct_path_connected, direct_path_connected/num_connected )
			else:
				return (num_connected, sum_path_connected, hub_path_weight, min_path, max_path, direct_path_connected, sum_path_direct_path_connected, max_path_direct_path_connected, hub_path_weight_direct_path_connected, float(sum_path_connected)/num_connected, hub_path_weight/num_connected, 'nan', 'nan', direct_path_connected/num_connected)
		else:
			return (num_connected, sum_path_connected, hub_path_weight, min_path, max_path, direct_path_connected, sum_path_direct_path_connected, max_path_direct_path_connected, hub_path_weight_direct_path_connected, 'nan', 'nan', 'nan', 'nan', 'nan')

	def has_user_entity(self, user_entity_id):
		return self.network.has_node(user_entity_id)
	

        def has_tag(self, userEntityID, tag):
            if self.uE_tags.has_key(userEntityID):
                if tag in self.uE_tags[userEntityID]:
                    return True
            return False


        def relation_has_tag(self, userEntityRelationID, tag):
	    """
	    userEntityRelationID is a tuple (userEntityRelationA, userEntityRelationB) where userEntityRelationA <= userEntityRelationB
	    """
	    if self.uER_tags.has_key(userEntityRelationID):
	        if tag in self.uER_tags[userEntityRelationID]:
		    return True
	    return False

	def hasTag(self, tag):
		return self.tag_uE.has_key(tag)

	def relationHasTag(self, tag):
	    return self.tag_uER.has_key(tag)

        def get_user_entity_tags(self, user_entity_id):
                return self.uE_tags.setdefault(user_entity_id, set())

	#def get_user_entity_relation_tags(self, user_entity_relation_id):
	def get_external_entity_relation_tags(self, external_entity_relation_id):
		#return self.uER_tags.setdefault(user_entity_relation_id, set())
		return self.uER_tags.setdefault(external_entity_relation_id, set() )    ### ALERT: Possible performance.... it creates unnecessary sets
	
        def get_user_entities_for_tag(self, tag):
	    if not self.tag_uE.has_key(tag):
		sys.stderr.write("Tag %s is not defined" %tag)
		return set()
	    return self.tag_uE[tag]

	#def get_user_entity_relations_for_tag(self, tag):
	def get_external_entity_relations_for_tag(self, tag):
	    return self.tag_uER[tag]

        def get_all_user_entity_tags(self):
            return self.tag_uE.keys()

        def get_all_user_entity_relation_tags(self):
            return self.tag_uER.keys()

	def getUserEntityGroups(self, userEntityID):
		"""
		Returns the groups in which this user entity belongs
		"""
		if self.uE_groups.has_key(userEntityID):
			return self.uE_groups[userEntityID]
		return set()
		
	def get_group_user_entities(self, group_id):
		"""
		Returns the user entities belonging to this group
		"""
		if self.relation_groups.has_key(group_id):
			return self.relation_groups[group_id]
		return set()

	def get_groups_ids(self):
		"""
		Returns a collection of the ids of the groups (they should be external entity ids)
		"""
		return self.relation_groups.keys()


	def has_group(self, group_id):
		return self.relation_groups.has_key(group_id)


        def get_level(self):
             return self.current_level

	def getSize(self):
		return self.network.number_of_nodes()

	def getNumberEdges(self):
		return self.network.number_of_edges()

        def getNetwork(self):
             return self.network

        def setNetwork(self, netw):
             self.network = netw

        def setIsNetworkCreated(self):
             """
             """
             self.is_network_created = True

        def isNetworkCreated(self):
             """
             Checks if is_network_created is True or False #Checks if this user entity set contains any edge (so, if it contains any edge)
             """
             #if self.network.number_of_edges()==0:
             #     return False
             #else:
             #     return True
             return self.is_network_created


	def selectUserEntitiesByTag(self, tagList):
		"""
		"""
		for current_tag in tagList:
			if self.hasTag(current_tag):
				self.select_user_entities( idUserEntityList = self.get_user_entities_for_tag( tag = current_tag ) )
			else:
				print "%s tag does not exist" %current_tag

	def selectUserEntityRelationsByTag(self, tagList):
		"""
		"""
		for current_tag in tagList:
			if self.relationHasTag(current_tag):
				self.selectUserEntityRelations( idUserEntityRelationList = self.get_user_entity_relations_for_tag( tag = current_tag) )

	def inverse_user_entity_relations_selection(self):
		"""
		Reverses the current selection of user entity relations
		"""

		#temp_all_relations = self.network.edges()  # returns a list of [(interactor1, interactor2)]

		#all_relations = set()

		#for  in temp_all_relations:
		#	for current_type in dictio.keys():
		#		all_relations.add((int1,int2,current_type))

		all_relations = set()
		
		for edge in self.network.edges_iter():
			all_relations.update(self.network.get_edge(edge[0],edge[1]))

		self.setIdUserEntityRelationSelected = all_relations.difference(self.setIdUserEntityRelationSelected)


	def select_user_entity_relations_between_two_subsets(self, user_entity_ids_list1, user_entity_ids_list2, types_list):
		"""
		Selects the direct relations between the user entities belonging to the two lists
		
		For example, an interaction between a user entity in list 1 with an user entity in list 2 will be selected
		
		An interaction between two user entities in list 1 won't be selected

		"types" is a list of the relation types to be selected
		"""
		
		selected_relations = []

		for current_1 in user_entity_ids_list1:
			for current_2 in user_entity_ids_list2:
				t = current_1
				if t > current_2:
					t = current_2
					current_2 = current_1
				
				for current_type in types_list:
					rel_id = (current_1,current_2)
					if self.network.has_edge(rel_id):
						selected_relations.append((current_1,current_2,current_type))

		self.selectUserEntityRelations(selected_relations)
		
	def select_user_entities(self, user_entity_id_list):
             """
             Selects a user entity in this user entity set 

	     User entities in the list must exist in the set!!!
             """
	     for idUserEntity in user_entity_id_list:
		     if self.network.has_node(idUserEntity):
			     self.setIdUserEntitySelected.add(idUserEntity)
		     #else:
			#     sys.stderr.write("Trying to select an unexisting user entity (%s) in user entity set %s" %(idUserEntity,self.id))


	def selectUserEntityRelations(self, external_entity_relation_id_list):
             """
             Selects a user entity relation in this user entity set

	     "external_entity_relation_id_list": is a list of external entity relation ids to select
             """

	     for current_eErID in external_entity_relation_id_list:
		     if self.eErIds2participants.has_key(current_eErID):
			     self.setIdUserEntityRelationSelected.add(current_eErID)
		     #else:
			#     sys.stderr.write("Trying to select an unexisting edge: %s\n" %current_eErID)
	     

	def _get_user_entity_relation_id(self, node1, node2):
		"""
		"""
		return ( min(node1,node2), max(node1,node2) )
		

	def getSelectionSize(self):
		return len(self.setIdUserEntitySelected)
	
	
	def getRelationSelectionSize(self):
		return len(self.setIdUserEntityRelationSelected)


	def clear_user_entity_selection(self):
	    self.setIdUserEntitySelected.clear()


	def clear_user_entity_relation_selection(self):
	    self.setIdUserEntityRelationSelected.clear()


	def addUserEntityId(self, idUserEntity, level=0):
             """
             Adds a user entity to this user entity set to the level "level"
             """

             if not self.network.has_node(idUserEntity):
                  self.network.add_node(idUserEntity)

                  if level<=self.current_level:
                       self.listLevelSetIdUserEntity[level].add(idUserEntity)
                  elif level==self.current_level+1:
                       self.current_level += 1
                       self.listLevelSetIdUserEntity.append(set([idUserEntity]))
                  else:
                       sys.stderr.write("Cannot add a node in this level")
                       return

                  self.size += 1

		  self.nodeLevelsDict[idUserEntity] = level

	     #else:
	#	     print "Warning. Trying to add an existing node: %s" %idUserEntity
             return


	def setGroupsHierarchy(self, list_hierarchy):
		"""
		"list_hierarchy" is a list [(group_id1, group_id2) in which group_id1 is child of group_id2
		"""

		[ self.groups_hierarchy.setdefault(current_child,current_parent) for (current_child, current_parent) in list_hierarchy ]


	def addUserEntitiesToGroup(self, group_id, userEntityID, representative_userEntityID, group_type, parentGroupID=None):
		"""
		One of the two user entities must belong to the network
		"parentGroupID" is used to store groups hierarchy
		"""

		if not self.network.has_node(representative_userEntityID):
			raise ValueError("Trying to add a user entity to a group with an unexisting representative")

		if not self.network.has_node(userEntityID):
			self.addUserEntityId(idUserEntity=userEntityID,level = self.getNodeLevel(representative_userEntityID)+1)

		self.uE_groups.setdefault(userEntityID,set()).add(group_id)
		self.relation_groups.setdefault(group_id,set()).add(userEntityID)
		self.relation_groups.setdefault(group_id,set()).add(representative_userEntityID)
		self.group_types.setdefault(group_id,group_type)
		if parentGroupID is not None:
			if self.groups_hierarchy.has_key(group_id):
				print "Group %s has more than one parent?" %(group_id)
			self.groups_hierarchy.setdefault(group_id,parentGroupID)
			
		#self._not_printed_groups.add(group_id)


	def get_external_entity_relation_participants(self, external_entity_relation_id):
		return self.eErIds2participants[external_entity_relation_id]

		
	def addUserEntityRelation(self, idUserEntity1, idUserEntity2, externalEntityRelationID):
             """
	     returns True if it has been added, or False if it existed previously a relation between idUserEntity1 and idUserEntity2
             """

	     # First, it is necessary to add the user entity ids to the correct levels
	     if not self.network.has_node(idUserEntity1) and not self.network.has_node(idUserEntity2):
		     raise ValueError("Trying to add a relation between two user entities (%s,%s) and any of them belong to this set" %(idUserEntity1,idUserEntity2))
	     
	     if not self.network.has_node(idUserEntity1):
		     self.addUserEntityId(idUserEntity=idUserEntity1, level = self.getNodeLevel(idUserEntity2)+1)
	     elif not self.network.has_node(idUserEntity2):
		     self.addUserEntityId(idUserEntity=idUserEntity2, level = self.getNodeLevel(idUserEntity1)+1)

	     if not self.eErIds2participants.has_key(externalEntityRelationID):
		     self.eErIds2participants[externalEntityRelationID] = [idUserEntity1, idUserEntity2]
	     else:
		     lista = self.eErIds2participants[externalEntityRelationID]
		     if idUserEntity1 not in lista:
			     lista.append(idUserEntity1)
		     if idUserEntity2 not in lista:
			     lista.append(idUserEntity2)


	     if self.network.has_edge(idUserEntity1, idUserEntity2):
		     if externalEntityRelationID not in set(self.network.get_edge(idUserEntity1, idUserEntity2)):
			     self.network.get_edge(idUserEntity1, idUserEntity2).append(externalEntityRelationID)
			     return True
	     else:
		     if self.getRestrictions("use_self_relations") is False:
			 if idUserEntity1 != idUserEntity2:
			     self.network.add_edge( idUserEntity1, idUserEntity2, [externalEntityRelationID] )
			     return True
			 else:
			     return False
		     else:
			 self.network.add_edge( idUserEntity1, idUserEntity2, [externalEntityRelationID] )
			 return True

	     return False

	     
	def _getRelationHash(self, userEntityID1, userEntityID2):
		"""
		"""
		if userEntityID1 <= userEntityID2:
			userEntityID1, userEntityID2 = userEntityID2, userEntityID1

		return (userEntityID1, userEntityID2)


	def remove_node(self, nodeID):
		
		for current_neighbor in self.network.neighbors(nodeID):
			self.remove_edge(nodeID,current_neighbor)
		self.network.delete_node(nodeID)
		level = self.nodeLevelsDict[nodeID]
		del self.nodeLevelsDict[nodeID]
		self.listLevelSetIdUserEntity[level].remove(nodeID)
		self.setIdUserEntitySelected.discard(nodeID)
		if self.uE_tags.has_key(nodeID):
			tags = self.uE_tags[nodeID]
			for current_tag in tags:
                            #self.tag_uE.remove(current_tag)
                            self.tag_uE[current_tag].remove(nodeID)
			del self.uE_tags[nodeID]
		self.size -= 1
		

	def remove_selected_relations(self):

		list_selected = list(self.setIdUserEntityRelationSelected) # makes a copy of the content because it is modified in the following loop

		[ self.remove_external_entity_relation(x) for x in list_selected ]

		self.setIdUserEntityRelationSelected.clear()


	def get_external_entity_relation_ids(self, user_entity_id1, user_entity_id2):
		return self.network.get_edge(user_entity_id1, user_entity_id2)


	def remove_external_entity_relation(self, externalEntityRelationID):

		participants = self.eErIds2participants[externalEntityRelationID]
		
		for current_participant1 in participants:
			for current_participant2 in participants:
				if self.network.has_edge(current_participant1, current_participant2):
					self.remove_edge(current_participant1, current_participant2, externalEntityRelationID)

	
	def remove_edge(self, nodeID1, nodeID2, externalEntityRelationID=None):
		"""
		If externalEntityRelationID is None, remove all edges between these two nodes
		"""

		nodeID1, nodeID2 = self._get_user_entity_relation_id(nodeID1, nodeID2)
		
		eErIds_list = self.network.get_edge(nodeID1, nodeID2)

		if externalEntityRelationID is None:                                          # Possible performance...
			for current_eErID in eErIds_list:
				self.remove_edge( nodeID1, nodeID2, current_eErID )
			return


		# Possible performance problem... we are running a list
		number = eErIds_list.count(externalEntityRelationID)
		if number == 0:
			sys.stderr.write("Trying to remove an unexisting edge: %s" %externalEntityRelationID)
		else:
			[ eErIds_list.remove(externalEntityRelationID) for x in xrange(number) ]
			if len(eErIds_list)==0:
				self.network.delete_edge(nodeID1, nodeID2)

			if self.uER_tags.has_key(externalEntityRelationID):
				for current_tag in self.uER_tags[externalEntityRelationID]:
					self.tag_uER[current_tag].remove(externalEntityRelationID)
				del self.uER_tags[externalEntityRelationID]
				
				#self.tag_uER[current_tag].remove(externalEntityRelationID)
				#if self.uER_tags.has_key((nodeID1,nodeID2,relation_type)):
			#	for current_tag in self.uER_tags[(nodeID1,nodeID2,relation_type)]:
			#		self.tag_uER[current_tag].remove((nodeID1,nodeID2,relation_type))
			#		del self.uER_tags[(nodeID1,nodeID2,relation_type)]
			

			#self.setIdUserEntityRelationSelected.discard((nodeID1,nodeID2))
			self.setIdUserEntityRelationSelected.discard(externalEntityRelationID)



	def remove_unconnected_nodes(self):

		for current_node in self.network.nodes():
			if self.get_node_degree(current_node) == 0:
				self.remove_node(current_node)
				

        def getNodeLevelDict(self):
		return self.nodeLevelsDict
             
	def getNodeLevel(self, nodeID):
		
		return self.nodeLevelsDict[nodeID]
	
	def get_user_entity_ids(self, level=None):
             """
             "level" can be a number with the level or "last"
             """
             
             if level is not None:
                  if str(level).lower()=="last":
                       level = self.current_level
                  if level<len(self.listLevelSetIdUserEntity):
                       return list(self.listLevelSetIdUserEntity[level])
                  return None
             else:
                  return self.network.nodes()


	def getRelationIds(self):
             """
             """
	     set_relations = set()
	     
	     for edge in self.network.edges():
		     set_relations.update(self.network.get_edge(edge[0],edge[1]))

             return set_relations


        def getRelations(self, cardinality_check=False):
		"""
		o cardinality_check: Boolean. Lowest userEntityId is returned as first element of the tuple; highest as second.
		"""
		
		#t = [ (x[0], x[1]) for x in self.network.edges() ]
		t = self.network.edges()

		# cardinality_check added by Joan
		if cardinality_check is False: return t
		arranged_t = []

		for i in t:
			if cardinality_check is True and not i[0] <= i[1]: arranged_t.append( (i[1], i[0]) )

		return arranged_t

	def getSelectedUserEntities(self):
		return self.setIdUserEntitySelected		

	def getSelectedUserEntityRelations(self):
		return self.setIdUserEntityRelationSelected

	def getRelationsOfSelectedUserEntities(self):
		return self.getRelationsBetweenGivenUserEntities(self.setIdUserEntitySelected)

	def getRelationsBetweenGivenUserEntities(self, setIdUserEntity):
		"""
		Fetchs interactions between all the user entities provided in the set
		Returns a set of tuples containing (id1, id2) 
		"""
		nUserEntity = len(setIdUserEntity)
		listRelations = []

		for edge in self.network.edges_iter(setIdUserEntity):
			if edge[0] in setIdUserEntity and edge[1] in setIdUserEntity:
				listRelations.append((edge[0],edge[1],self.network.get_edge(edge[0],edge[1])))
		return listRelations
			
	def getShortestPathBetween(self, idUserEntity1, idUserEntity2):
		if idUserEntity1<idUserEntity2:
			idUserEntity2,idUserEntity1 = idUserEntity1, idUserEntity2
		return self.calculated_shortest_paths.setdefault((idUserEntity1, idUserEntity2),graph_utilities.get_shortest_path_between(self.network, idUserEntity1, idUserEntity2))
		
	def getAllPathsFrom(self, idUserEntity):
		return graph_utilities.get_all_paths_from(self.network, idUserEntity)

 
	def getUnionWithGivenUserEntitySet(self, objUserEntitySet, flagIncludeInteractions):

		union_user_entity_ids_set = set( self.get_user_entity_ids() ).union( set( objUserEntitySet.get_user_entity_ids() ) )

		listLevelSetId = [ union_user_entity_ids_set ]

		listInteraction = []

		if flagIncludeInteractions:

			# Merge relation ids
			dict_relations = {}
			all_edges = set()
			#for node1, node2, eErIds_list in self.network.edges():
			for node1, node2 in self.network.edges_iter():
				eErIds_list = self.network.get_edge(node1,node2)
				dict_relations[self._getRelationHash(node1,node2)] = copy.deepcopy(eErIds_list)

			#for node1, node2, eErIds_list in objUserEntitySet.network.edges():
			for node1, node2 in objUserEntitySet.network.edges():
				eErIds_list =  objUserEntitySet.network.get_edge(node1,node2)
				dict_relations.setdefault(self._getRelationHash(node1,node2), []).extend(eErIds_list)

			for key, value in dict_relations.iteritems():
				listInteraction.append((key[0],key[1],list(set(value))))
				
		return (listLevelSetId, listInteraction)

			
	def getIntersectionWithGivenUserEntitySet(self, objUserEntitySet, flagIncludeInteractions):
		"""
		Interactions for intersection nodes are included if flagIncludeInteractions is selected
		"""

		intersection_user_entity_ids_set = set( self.get_user_entity_ids() ).intersection( set( objUserEntitySet.get_user_entity_ids() ) )

		listLevelSetId = [ intersection_user_entity_ids_set ]

		listInteraction = []

		if flagIncludeInteractions:
			dict_relations = {}

			#for current_node1, current_node2, eErIdsList in self.network.edges_iter(intersection_user_entity_ids_set):
			for current_node1, current_node2 in self.network.edges_iter(intersection_user_entity_ids_set):
				if( objUserEntitySet.network.has_node(current_node1) and objUserEntitySet.network.has_node(current_node2) ):
					eErIdsList = self.network.get_edge(current_node1, current_node2)
					dict_relations[self._getRelationHash(current_node1,current_node2)] = copy.deepcopy(eErIdsList)
					#listInteraction.append(current_edge)

			#for current_node1, current_node2, eErIdsList in objUserEntitySet.network.edges_iter(intersection_user_entity_ids_set):
			for current_node1, current_node2 in objUserEntitySet.network.edges_iter(intersection_user_entity_ids_set):
				if( self.network.has_node(current_node1) and self.network.has_node(current_node2) ):
					eErIdsList = objUserEntitySet.network.get_edge(current_node1, current_node2)
					dict_relations.setdefault(self._getRelationHash(current_node1,current_node2), []).extend(eErIdsList)
					#listInteraction.append(current_edge)
			
			for key, value in dict_relations.iteritems():
				listInteraction.append((key[0],key[1],list(set(value))))
				

		return (listLevelSetId, listInteraction)
	

#	def outputUserEntitySet(self, fileName, format):
#		if format=="sif":
#			fileOut = open(fileName, 'w')
#			for actual_node in self.network.nodes():
#				if( len(self.network.neighbors(actual_node))==0 ):
#					fileOut.write("%s\n" %actual_node)
#			for id1, id2, x in self.network.edges():
#				fileOut.write("%d %s %d\n" % (id1, x, id2))
#			
#		elif format=="dot":
#			networkx.write_dot(self.network, fileName) 
#		return


        def get_node_degree(self, userEntityID):

             return self.network.degree(userEntityID)

        def get_unconnected_nodes(self):
             return_set = set()
             for current_uE in self.get_user_entity_ids(level=0):
                  if self.network.degree(current_uE)==0:
                       return_set.add(current_uE)

             return return_set


#        def outputNetwork(self, outmethod, format):
#             if format == "sif":
#                  for actual_node in self.network.nodes():
#                       if( len(self.network.neighbors(actual_node))==0 ):
#                            outmethod("%s\n" %actual_node)
#                  for id1, id2, x in self.network.edges():
#                       outmethod("%d %s %d\n" % (id1, x, id2))

#             elif format == "table":
#                  columns = ["interactor1","interactor2","interaction_type"]
#                  values = []
#                  for actual_node in self.network.nodes():
#                       if( len(self.network.neighbors(actual_node))==0 ):
#                            values.append((str(actual_node),"-","-"))
#                  for id1, id2, x in self.network.edges():
#                       values.append((str(id1), str(x), str(id2)))

#                  th_str = "<tr>%s</tr>" %"".join([ "<th>%s</th>" %x for x in columns ])
#                  data_str = "".join([ "<tr>%s</tr>" %"".join([ "<td>%s</td>" %current_column for current_column in current_row ]) for current_row in values ])
#                  
#                  outmethod("<table>%s%s</table>" %(th_str,data_str))

	def __repr__(self):
		return "%s / %s" % (self.listLevelSetIdUserEntity, self.network.edges())

        def __str__(self):
		str_list = ["User Entity Set %s" %self.id]
		str_list.append("\tNumber of user entities: %s" %self.network.number_of_nodes())
		str_list.append("\tNumber of relations: %s" %(self.network.number_of_edges()))
		str_list.append("\tLevels: %s" %(len(self.listLevelSetIdUserEntity)-1))
		str_list.append("\tNumber of root entities: %s" %(len(self.listLevelSetIdUserEntity[0])))

		str_list.append("\tRestrictions:")
		for current_restriction, restrictions in self.restrictions_dict.iteritems():
			if isinstance(restrictions,list):
				if len(restrictions)>0:
					str_list.append("\t\t%s: %s" %(current_restriction, restrictions))
				else:
					str_list.append("\t\t%s: Not applied" %current_restriction)
			else:
				str_list.append("\t\t%s: %s" %(current_restriction, restrictions))

		return "\n".join(str_list)

	def _get_xml_header(self):
		return "<user_entity_set id=\"%s\">" %self.id

	def _get_xml_foot(self):
		return "</user_entity_set>"

	def _get_xml(self, inner_content=""):
		return "%s%s%s" %(self._get_xml_header(),
				  inner_content,
				  self._get_xml_foot())

	def _get_groups_xml(self, only_not_printed=False, group_identifiers={}):
		"""
		group_identfiers is a dictionary with the id of the node as a key and a name for the node
		"""

		if only_not_printed is True:
			#toprint = list(self._not_printed_groups)
			toprint = set(self.relation_groups).difference(self._printed_groups)
		else:
			toprint = list(self.relation_groups)
			
		xml_list = []

		for current_group in toprint:
			#self._not_printed_groups.discard(current_group)
			#print self._not_printed_groups
			self._printed_groups.add(current_group)
			if self.groups_hierarchy.has_key(current_group):
				parent_str = "parentGroupID=\"%s\" parentGroupName=\"%s\"" %(self.groups_hierarchy[current_group],
											     str(group_identifiers.setdefault(self.groups_hierarchy[current_group],self.groups_hierarchy[current_group])).replace("<","\<").replace(">","\>"))
			else:
				parent_str = ""

			xml_list.append("<add_new_group groupID=\"%s\" groupName=\"%s\" groupType=\"%s\" node_ids=\"%s\" %s/>" %(current_group,str(group_identifiers.setdefault(current_group,current_group)).replace("<","\<").replace(">","\>"),self.group_types[current_group],",".join(map(str,self.relation_groups[current_group])),parent_str))

		return "".join(xml_list)


	def _update_selected_external_ids_dict(self, attribute, selected_ext_id_name):
		"""
		Updates self.selectedExternalIdsDict and self.selectedExtIdsLowerDict
		"""

		if selected_ext_id_name is None: return

		if not self.selectedExternalIdsDict.has_key(attribute): self.selectedExternalIdsDict[attribute] = set([])
		if not self.selectedExtIdsLowerDict.has_key(attribute): self.selectedExtIdsLowerDict[attribute] = {}
		self.selectedExternalIdsDict[attribute].add(selected_ext_id_name)
		self.selectedExtIdsLowerDict[attribute][selected_ext_id_name.lower()] = selected_ext_id_name
	
	
