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


class UserEntity(object):
    """
    A container class to represent a User Entity

    A user entity is composed by one or several externalEntities

    It only contains the collection of external entity ids
    """

    def __init__(self, id, listExternalEntityId=None):
        """
        externalEntities
        """
        # Set that stores the identifiers of externalEntities

        self.id = id

        if listExternalEntityId is None:
            self.setExternalEntityId = set()
        else:
            self.setExternalEntityId = set(listExternalEntityId)

        return
    
    def addExternalEntity(self, externalEntityId):
        self.setExternalEntityId.add(externalEntityId)
        return
    
    def containsExternalEntity(self, externalEntityId):
        if externalEntityId in self.setExternalEntityId:
            return True
        return False
    
    def getSize(self):
        return len(self.setExternalEntityId)

    def __str__(self):
        return "User Entity with id %s,composed by %s external Entities" %(self.id, self.getSize())
        #return "%s" % self.id

    def get_externalEntitiesIds_set(self):
        return self.setExternalEntityId
    
