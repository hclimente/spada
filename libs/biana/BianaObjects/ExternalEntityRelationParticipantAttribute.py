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

import sets
import ExternalEntityAttribute

class ExternalEntityRelationParticipantAttribute(ExternalEntityAttribute.ExternalEntityAttribute):

    
    def __init__(self, attribute_identifier, value):
        ExternalEntityAttribute.ExternalEntityAttribute.__init__(self, attribute_identifier = attribute_identifier, value = value )


    def isValidType(type, biana_database):

        return type.lower() in biana_database.VALID_EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TYPES_SET

    isValidType = staticmethod(isValidType)


