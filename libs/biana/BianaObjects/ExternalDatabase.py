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

import time
import sets


class ExternalDatabase(object):
    """
    Class to represent an external database (any database giving biologic relevant information)
    """


    def __init__(self, databaseName, databaseVersion, databaseFile, databaseDescription, defaultExternalEntityAttribute, databaseDate=None, externalDatabaseID=None, isPromiscuous=False):
        """
        Initializes a ExternalDatabase Object

        "databaseName" is the name for the external database (i.e. swissprot, IntAct, Reactome...)
        
        "datavaseVersion" is the version of the database. The combination of databaseName and databaseVersion must be unique!

        "databaseFile" is the parsed database file. If there are multiple files, it is an empty string

        "databaseDescription" is a longer description of the database (for example, "public repository of interactions of ncbi", etc.)

        "externalDatabaseID" is a unique identifier for the database. When parsing it should not be added, only when this information is persistent in the database

        "defaultExternalEntityAttribute" is the default attribute type of the data an external database is providing

        "databaseDate" Parsing date. If it is None, automatically assigns to current date to it
	"isPromiscuous" Flag deciding whether database gives information that is going to be added to more than one user entiries
        """

        self.databaseName = databaseName.lower()
        self.databaseVersion = databaseVersion.lower()
        self.databaseFile = databaseFile
        self.databaseDescription = databaseDescription
        self.defaultExternalEntityAttribute = defaultExternalEntityAttribute 
	self.isPromiscuous = isPromiscuous

        self.externalDatabaseID = externalDatabaseID         # Only it is setted to a value when it is inserted or loaded from the database. If it is not loaded, it is None

        if( databaseDate is None ):
            date = time.localtime()
            actual_date = "%s-%s-%s" %(date[0],date[1],date[2])
            self.databaseDate = actual_date
        else:
            self.databaseDate = databaseDate

        self.valid_eE_attributes = sets.Set()
        self.valid_eEr_attributes = sets.Set()
        self.valid_eE_types = sets.Set()
        self.valid_eEr_types = sets.Set()

        self.parsing_time = None   # Parsing time is set to None by default

    def add_valid_external_entity_attribute_type(self, attribute_identifier):
        self.valid_eE_attributes.add(attribute_identifier.lower())

    def add_valid_external_entity_relation_attribute_type(self, attribute_identifier):
        self.valid_eEr_attributes.add(attribute_identifier.lower())

    def add_valid_external_entity_type(self, eE_type):
        self.valid_eE_types.add(eE_type.lower())

    def add_valid_external_entity_relation_type(self, eEr_type):
        self.valid_eEr_types.add(eEr_type.lower())

    def get_valid_external_entity_attribute_type(self):
        return self.valid_eE_attributes

    def get_valid_external_entity_relation_attribute_type(self):
        return self.valid_eEr_attributes

    def get_valid_external_entity_type(self):
        return self.valid_eE_types

    def has_external_entity_type(self, eE_type):
        """
        Checks if this external database contains external entities of the type "eE_type"
        """
        return self.valid_eE_types.has_key(eE_type)

    def get_valid_external_entity_relation_type(self):
        return self.valid_eEr_types

    def get_name(self):
        return self.databaseName

    def get_version(self):
        return self.databaseVersion

    def get_parsed_file(self):
        return self.databaseFile

    def get_parsing_date(self):
        return self.databaseDate

    def get_description(self):
        return self.databaseDescription

    def get_default_eE_attribute(self):
        return self.defaultExternalEntityAttribute

    def get_id(self):
        return self.externalDatabaseID

    def set_id(self, externalDatabaseID):
        self.externalDatabaseID = externalDatabaseID

    def set_parsing_time(self, time):
        self.parsing_time = time

    def get_parsing_time(self):
        return self.parsing_time

    def get_promiscuity(self):
	return self.isPromiscuous

    def __str__(self):
        return "%s [%s]" %(self.databaseName,self.databaseVersion)
