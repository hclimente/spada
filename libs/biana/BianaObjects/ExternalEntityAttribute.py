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
import re # for checking given value's integrity


ATTRIBUTE_IDENTIFIER_TO_OBJECT = {}


class ExternalEntityAttribute(object):
    """
    Class that represents an external entity attribute
    """

    
    #attr_to_len = {}

    def __init__(self, attribute_identifier, value, type=None, version=None, additional_fields={}):

        # Check for a correct value
        self.value = value
        self.type = type
        self.version = version
        self.additional_fields = additional_fields
        self.attribute_identifier = attribute_identifier
	# Track lengths of attribute values
	#if self.attr_to_len.setdefault(attribute_identifier, 0) < len(str(value)):
	#    self.attr_to_len[attribute_identifier] = len(str(value))
	#    print self.attribute_identifier, len(str(value))
	#print self


    def get_field(self, field_name):

        if( field_name == "value" ):
            return self.value
        elif( field_name == "type" ):
            return self.type
        elif( field_name == "version" ):
            return self.version
        else:
            return self.additional_fields[field_name]


    def __str__(self):

        temp1 = ["value: %s" %self.value]
        if( self.type is not None ):
            temp1.append("type: %s" %self.type)
        if( self.version is not None ):
            temp1.append("version: %s" %self.version)
        temp2 = [ ("%s: %s") %(x,self.additional_fields[x]) for x in self.additional_fields ]

        return "Attribute %s. %s %s" %(self.attribute_identifier,"\t".join(temp1),"\t".join(temp2))

    def __repr__(self):
        return self.__str__()



    def __eq__(self, other):
        """Other can be either string of dictionary however comparison below is generic enough"""
        return str(self.value).lower() == str(other.value).lower()




    def isFullTextSearchable(attribute_identifier, biana_database):

        # Possible bug: it is only applyed for general attributes. If a specific attribute has to have a keyword index, this won't work
        return attribute_identifier.lower() in biana_database.VALID_EXTERNAL_ENTITY_DESCRIPTIVE_SEARCHABLE_ATTRIBUTE_TYPES_SET
 
    isFullTextSearchable = staticmethod(isFullTextSearchable)




    def isNumericAttribute(attribute_identifier, biana_database):
        return attribute_identifier.lower() in biana_database.VALID_EXTERNAL_ENTITY_NUMERIC_ATTRIBUTE_TYPES_SET

    isNumericAttribute = staticmethod(isNumericAttribute)


    def isSpecialAttribute(attribute_identifier, biana_database):
        return attribute_identifier.lower() in biana_database.VALID_EXTERNAL_ENTITY_SPECIAL_ATTRIBUTE_TYPES_SET

    isSpecialAttribute = staticmethod(isSpecialAttribute)


    def isValidAttribute(attribute_identifier, biana_database):

        return attribute_identifier.lower() in biana_database.VALID_EXTERNAL_ENTITY_ATTRIBUTE_TYPES_DICT

    isValidAttribute = staticmethod(isValidAttribute)



    def isIdentifierType(attribute_identifier, biana_database):


        return attribute_identifier.lower() in biana_database.VALID_EXTERNAL_ENTITY_IDENTIFIER_ATTRIBUTE_TYPES_SET or attribute_identifier.lower() in biana_database.VALID_EXTERNAL_ENTITY_VERSIONABLE_IDENTIFIER_ATTRIBUTE_TYPES_SET

    isIdentifierType = staticmethod(isIdentifierType)



    def isVersionableIdentifier(attribute_identifier, biana_database):

        return attribute_identifier.lower() in biana_database.VALID_EXTERNAL_ENTITY_VERSIONABLE_IDENTIFIER_ATTRIBUTE_TYPES_SET

    isVersionableIdentifier = staticmethod(isVersionableIdentifier)

    


    
    




