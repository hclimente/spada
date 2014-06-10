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


# Each participant in the relation can have different attributes: role, detection method,...
import ExternalEntity

class ExternalEntityRelation(ExternalEntity.ExternalEntity):
    """
    This general object represents an entry represeting a relation between external entities in an external database

    An external Entity Relation is a relation of any type between two or more external entities

    It can contain several attributes.
    """


    def __init__(self, source_database, relation_type, id=None):
        """
        "source_databas_id" is the source database id where this entity is described

        "relation_type" indicates the type of externalEntityRelation

        "id" is the UNIQUE identifier in the database for this external entity
        """

        self.participants = {}   # This is a dictionary with the following format:
                                 # key: externalEntityID for the participant
                                 # values: dictionary with its attributes (role, detection method,...)

        self.relation_type = relation_type

        ExternalEntity.ExternalEntity.__init__(self, source_database = source_database, type="relation", id=id)


    def add_participant(self, externalEntityID):
        """
        Adds a participant to this relation

        "externalEntityID" externalEntityID corresponding to the participant
        """

        if self.participants.has_key(externalEntityID):
            return
        else:
            self.participants[externalEntityID] = {}




    def get_participant_external_entity_ids_list(self):
        """
        Returns a list with all the externalEntityIds of the participants
        """
        return self.participants.keys()




    def get_participant_attributes(self, participantExternalEntityID):
        """
        Returns a list of [(attribute_name, fieldValues dictionary)]
        """
        return self.participants[participantExternalEntityID]



    def get_participant_attribute(self, participantExternalEntityID, attribute_identifier):
        """
        """
        try:
            return self.participants[participantExternalEntityID][attribute_identifier.lower()]
        except:
            return []



    def add_participant_attribute(self, externalEntityID, participantAttribute ):
        """
        Adds an attribute to a participant in the relation

        If the participant didn't exist previously, it gives error. It is necessary to insert previously the participant!
        """
    
        if self.participants.has_key(externalEntityID):
            try:
                self.participants[externalEntityID][participantAttribute.attribute_identifier.lower()].append(participantAttribute)
            except KeyError:
                self.participants[externalEntityID][participantAttribute.attribute_identifier.lower()] = [participantAttribute]
        else:
            raise ValueError("Trying to add attributes to an unexisting participant. Participant: %s" %externalEntityID)




    def get_relation_type(self):
        """
        Returns the type of this relation (interaction, reaction...)
        """
        return self.relation_type




    def get_valid_relation_types(biana_database):
        """
        Static method

        Returns a list with all available relation types
        """
        return biana_database.VALID_EXTERNAL_ENTITY_RELATION_TYPES_DICT.values()

    get_valid_relation_types = staticmethod(get_valid_relation_types)




    
    def isValidType(type, biana_database):

        return type.lower() in ExternalEntityRelation.VALID_EXTERNAL_ENTITY_RELATION_TYPES_DICT

    isValidType = staticmethod(isValidType)
