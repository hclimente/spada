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
import ExternalEntity
import ExternalEntityAttribute


class Ontology(ExternalEntity.ExternalEntity):
    """
    Class to represent a general ontology

    The Ontology is an external entity itself. Each of its elemens is also an ExternalEntityObject

    Each ontology as a "linkedAttribute", which is the primary attribute to represent the ontology (for example, taxID for taxonomy), and a "descriptionAttribute", which is an attribute that describes the external entity element
    """
    
    def __init__(self, source_database, name, linkedAttribute, descriptionAttribute, id=None, levelAttribute=None ):
        """
        "source_database" is the source database id where this entity is described

        "name" is the name for the ontology. It must be UNIQUE! There cannot be different ontologies with the same name

        "linkedAttribute" is the attribute_identifier for the primary attribute of the ontology (for example, taxID for taxonomy ontology)
        
        "descriptionAttribute" is the attribute_identifier for representing a human readable description of the element. This attribute is used when showing the ontolgy to users

        "id" is the UNIQUE identifier in the database for this external entity (as the ontology is an External Entity)
        """

        self.name = name
        self.trees = {} # is_a in the first position, is_part_of in the second one
        self.hierarchy =  {}  # vertical hierarchy
        self.root_ids = set()
        self.linked_attribute = linkedAttribute
        self.level_attribute = levelAttribute
        self.linked_attribute_values = {}
        self._attrID2id = {}
        self.description_attribute = descriptionAttribute
        self.externalEntityObjs = {}
        self.all_ids = set()

        self.precalculated_descendants = {}

        ExternalEntity.ExternalEntity.__init__(self, source_database = source_database, type="ontology", id=id)
    

    def add_element(self, ontologyElementID, isA=[], isPartOf=[], linkedAttributeValue=None):
        """
        Adds an element to the ontology.

        "ontologyElementID": externalEntityID that identifies the externalEntityObject belonging to the ontology
        
        "isA": list with the externalEntityIDs of the parents of this element

        "isPartOf": list with the externalEntityIDs of the elements to which this element is part of

        "linkedAttributeValue" is the value of the main attribute of the added external entity. Not mandatory.
        """

        self.all_ids.add(ontologyElementID)
        self.linked_attribute_values[ontologyElementID] = linkedAttributeValue
        self._attrID2id[linkedAttributeValue]=ontologyElementID
        self.hierarchy.setdefault(ontologyElementID,[])

        if( len(isA)==0 ):
            self.root_ids.add(ontologyElementID)

        self.trees[ontologyElementID] = (isA,isPartOf)

        for current_parent in isA:
            self.hierarchy.setdefault(current_parent,[]).append(ontologyElementID)


    def _set_external_entities_dict(self, externalEntitiesDict):
        """
        Sets the external entity objects corresponding to the elements of the ontology

        "externalEntitiesDict": Dictionary with all the external entities. Key: externalEntityID. Value: externalEntity Object

        Objects are only required for printing the ontology
        """
        self.externalEntityObjs = externalEntitiesDict

    def get_all_external_entity_ids(self):
        return self.all_ids


    def linkedAttrID2ID(self, attr_id):
        return self._attrID2id[attr_id]
        
    
    def get_descendants(self, ontologyElementID):
        """
        Gets all the descendants, using the "is_a" relation
        """

        if self.precalculated_descendants.has_key(ontologyElementID):
            return self.precalculated_descendants[ontologyElementID]

        result = set()
        #result = []

        for current_descendant_id in self.hierarchy[ontologyElementID]:
            if current_descendant_id == ontologyElementID:
                sys.stderr.write("Ontology has a loop. An element %s [%s] is a child of itself?\n" %(current_descendant_id,self.linked_attribute_values[current_descendant_id]))
                return result
            else:
                if current_descendant_id not in result:
                    result.add(current_descendant_id)
                    result.update(self.get_descendants(current_descendant_id))
                    # result.update(self.get_descendants(current_descendant_id))
                else:
                    sys.stderr.write("Ontology has a loop, between %s [%s] and %s [%s]\n" %(current_descendant_id, 
                                                                                            self.linked_attribute_values[current_descendant_id],
                                                                                            ontologyElementID,
                                                                                            self.linked_attribute_values[ontologyElementID]))

        self.precalculated_descendants[ontologyElementID] = result

        return result


    def get_linked_attr_and_description_tuples(self):
        """
        Returns a list of tuples with the format: (linked_attr, descriptive_attr)
        """

        return [ ( self.linked_attribute_values[x],"|".join([ y.value for y in self.externalEntityObjs[x].get_attribute(self.description_attribute)]) ) for x in self.linked_attribute_values ]

    def get_all_linked_attributes(self):
        """
        Returns a list with the main attribute for all the elements in the ontology
        """
        return self.linked_attribute_values.values()

    def get_all_external_entity_ids(self):
        """
        Returns a list with the external entity ids of all the elements in the ontology
        """
        return self.linked_attribute_values.keys()

    def has_element(self, linkedAttributeID):
        """
        Returns a boolean indicating if an external entity with this attribute is found in the ontology
        """
        return linkedAttributeID in self.linked_attribute_values

    def get_parents_ids(self, elementID):
        """
        Returns a list with the parents of the element with this externalEntityID (using the relation is_a)
        """
        return self.trees[elementID][0]

    def get_part_parents_ids(self, elementID):
        """
        Returns a list with the parents of the element with this externalEntityID (using the relation is_part_of)
        """
        return self.trees[elementID][1]


    def _recursive_tree_print(self, id, outmethod, depth=0):
        """
        Prints recursively in stdout a tree representing the ontology, using the external entity id.

        Only for testing purposes
        """
        
        for x in xrange(depth):
            outmethod("\t")
        #outmethod(str(id))
        outmethod("%s [%s]" %('|'.join([ x.value for x in self.externalEntityObjs[id].get_attribute(self.description_attribute) ]),self.linked_attribute_values[id]))
        outmethod("\n")
        depth = depth+1
        for current_descendant_id in self.hierarchy[id]:
            self._recursive_tree_print(current_descendant_id, outmethod, depth)


    def print_tree(self, outmethod=sys.stdout.write):
        """
        Prints recursively in stdout a tree representing the ontology, using the external entity id.

        Only for testing purposes
        """

        #print self.root_ids
        
        for x in self.root_ids:
            self._recursive_tree_print( id = x,
                                        depth = 0,
                                        outmethod = outmethod )


    def _recursive_tree_xml(self, id):
            
        nodes_xml = [ "<node ontologyNodeID=\"%s\" id=\"%s\">" %(self.linked_attribute_values[id],
                                                                 '|'.join([ x.value for x in self.externalEntityObjs[id].get_attribute(self.description_attribute) ])) ]
        
        for current_descendant_id in self.hierarchy[id]:
            nodes_xml.extend(self._recursive_tree_xml(current_descendant_id))
            
        nodes_xml.append("</node>")
        return nodes_xml
        

    def get_xml(self):
        """
        Returns a String with an XML representation of the ontology

        To execute this method, is necessary to load in the ontology all the external Entity Objects
        """

        nodes_xml = ["<ontology name=\"%s\">" %self.name]
        
        for x in self.root_ids:
            nodes_xml.extend(self._recursive_tree_xml( id = x) )

        nodes_xml.append("</ontology>")

        return "\n".join(nodes_xml)


    def _traverse_tree_for_leafs(self, id):
        """
            Helper function to reach leafs traversing the tree 
        """
        if len(self.hierarchy[id]) == 0:
            #print self.linked_attribute_values[id], self.hierarchy[id], self.trees[id]
            nodes = [ ", ".join([ x.value for x in self.externalEntityObjs[id].get_attribute(self.description_attribute) ]) ]
        else:
            nodes = []
        
        for current_descendant_id in self.hierarchy[id]:
            nodes.extend(self._traverse_tree_for_leafs(current_descendant_id))
            
        return nodes
        

    def get_leafs(self):
        """
        Returns a list of leafs in the ontology tree
        To execute this method, is necessary to load in the ontology all the external Entity Objects
        """
        leaf_nodes = []
        
        for x in self.root_ids:
            leaf_nodes.extend(self._traverse_tree_for_leafs( id = x) )

        return leaf_nodes

