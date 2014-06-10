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
import biana.biana_globals as BIANA_GLOBALS


class ExternalEntity(object):
    """
    This general object represents an entry in an external database

    It can contain several attributes.
    """

    # to remove
    PROMISCUOUS_EXTERNAL_ENTITY_TYPES_DICT = dict(BIANA_GLOBALS.PROMISCUOUS_EXTERNAL_ENTITY_TYPES_DICT)


    def __init__(self, source_database, type, id=None):
        """
        "source_database" is the source database id wbere this entity is described

        "type" indicates the type of externalEntity.It can be:
                              'protein','gene','protein alignment','cluster','structure, SCOPElement'

        "id" is the UNIQUE identifier in the database for this external entity
        """

        type = type.lower()

        self.id = id           # This value will be "None" if still not inserted in the database. It has been inserted in the database, it MUST be a numeric value
        self.attributes = {}
        self.sourceDatabase = source_database
        self.type = type

        return


    def __str__(self):
        
        string = "External Entity ID: %s\n" %self.id
        
        for current_attribute in self.attributes:
            string += "%s:\n\t%s\n" %(current_attribute, "\t".join(map(str,self.attributes[current_attribute])))

        #print string
        return string



    def __eq__(self, other):

        return self.id==other.id



    def get_id(self):
        """
        Gets the UNIQUE identifier
        """
        return self.id



    def get_source_database(self):
        """
        Gets the ExternalDatabase Object where this externalEntity is contained
        """
        return self.sourceDatabase



    def get_type(self):
        """
        Gets a string representing the type of this externalEntity
        """
        return self.type



    def set_id(self, id_value):
        if self.id is None:
            self.id = id_value
        else:
            raise ValueError("Cannot change identifier from externalEntity")



    def add_attribute(self, externalEntityAttribute):
        """
        Adds an attribute to this ExternalEntity Object
        
        "externalEntityAttribute": A valid externalEntityAttribute Object. It cannot contain value None as its value
        """
        
        if externalEntityAttribute.value is None:
            sys.stderr.write("Attribute %s not added because it has None value\n" %(externalEntityAttribute.attribute_identifier))
            return

        self.attributes.setdefault(externalEntityAttribute.attribute_identifier.lower(),set([])).add(externalEntityAttribute)



    def get_attribute(self, attribute_identifier):
        """
        Returns a SET with the externalEntityAttribute objects corresponding to the attribute_identifier associated to this externalEntity.

        If the externalEntity object does not contain any attribute of "attribute_identifier" type, it returns an empty list
        """

        try:
            return self.attributes[attribute_identifier.lower()]
        except KeyError:
            return set()

        
    def get_attributes_dict(self):
        """
        Returns a dictionary with all the ExternalEntityAttributes associated to this externalEntity.

        Keys correspond to "attribute_identifier" and values to a Set with all external entity objects
        """
        return self.attributes



    def isValidType(type, biana_database):

        return type.lower() in biana_database.VALID_EXTERNAL_ENTITY_TYPES_DICT

    isValidType = staticmethod(isValidType)



    def get_valid_external_entity_types(biana_database):
        
        return biana_database.VALID_EXTERNAL_ENTITY_TYPES_DICT.values()

    get_valid_external_entity_types = staticmethod(get_valid_external_entity_types)




    
