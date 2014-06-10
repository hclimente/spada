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

class UnificationAtomElement(object):
    """
    Class to represent a condition of unification
    """

    def __init__(self, externalDatabaseID_A, externalDatabaseID_B, externalAttribute, field_conditions_A=[], field_conditions_B=[],
                 field_cross_references = ["value"]):
        """
        "externalDatabaseID_A" and "externalDatabaseID_B" are the databaseIntegers that identify the databases we want to crossreference
                                          externalDatabaseID_A must be lower or equal to externalDatabaseID_B. If not, they will be flipped
        
        "externalAttribute" is the attribute that we want to cross-reference
                            It can be a list of attributes if they have to be taken into account together (for example, ["sequence","taxID"]
        
        "fields_A" and "fields_B" is a list with the fields we want to cross-reference (for example, ["value","type"]
        
        "field_conditions_A" and "field_conditions_B" is optional, and they must be a dictionary with restricted fields { "field": "restricted_value" }
        
        "field_cross_references" is a list of the fields that have to be cross-referenced. By default, it takes the value [ "value" ]
        """

        # I DON'T KNOW IF DOING THIS IN THIS WAY...
        if externalDatabaseID_A > externalDatabaseID_B:
            self.externalDatabaseID_A = externalDatabaseID_B
            self.externalDatabaseID_B = externalDatabaseID_A
            self.field_conditions_A = field_conditions_B
            self.field_conditions_B = field_conditions_A
        else:
            self.externalDatabaseID_A = externalDatabaseID_A
            self.externalDatabaseID_B = externalDatabaseID_B
            self.field_conditions_A = field_conditions_A
            self.field_conditions_B = field_conditions_B
            
        self.field_cross_references = field_cross_references

        if (isinstance(externalAttribute,list) ):
            self.externalAttribute = [ x.lower() for x in externalAttribute ]
        else:
            self.externalAttribute = [externalAttribute.lower()]

        #for actual_external_attribute in self.externalAttribute:
        #    if NOT_ALLOWED_CROSS_REFERENCES.has_key(actual_external_attribute):
        #        raise ValueError("%s cannot be used for integration purposes" %externalAttribute)

        # We sort them to be able to compare later distinct unification atom elements
        self.field_cross_references.sort()
        self.field_conditions_A.sort()
        self.field_conditions_B.sort()


    def __str__(self):
        """
        String representation of the class
        """

        return_list = ["Unification atom"]
        return_list.append("Cross databases: %s - %s" %(self.get_externalDatabaseID_A(),self.get_externalDatabaseID_B()))
        return_list.append("Using external attribute: %s" %self.get_external_attribute_list())
        return_list.append("Field cross_references: %s" %self.field_cross_references)
        
        return "%s" %("\n\t".join(return_list))


    def get_externalDatabaseIDs(self):
        """
        returns a tuple with both database identifiers
        """
        return (self.externalDatabaseID_A, self.externalDatabaseID_B)

    def get_externalDatabaseID_A(self):
        return self.externalDatabaseID_A

    def get_externalDatabaseID_B(self):
        return self.externalDatabaseID_B

    def get_external_attribute_list(self):
        return self.externalAttribute

    def get_field_conditions_A(self):
        """
        gets a list of tuples with the following format: (field,value)
        """
        return [ (x,self.field_conditions_A[x]) for x in self.field_conditions_A ]


    def get_field_conditions_B(self):
        """
        gets a list of tuples with the following format: (field,value)
        """
        return [ (x,self.field_conditions_B[x]) for x in self.field_conditions_B ]


    def get_field_cross_references(self):

        return self.field_cross_references


    def is_comparable(self, unificationAtomElement2 ):
        """
        Checks if two atom elements are comparable, by checking the equivalence of all its attributes except the source databases
        """

        return ( self.field_cross_references == unificationAtomElement2.get_field_cross_references() and
                 ( (self.field_conditions_A == unificationAtomElement2.get_field_conditions_A() and
                    self.field_conditions_B == unificationAtomElement2.get_field_conditions_B() ) or
                   (self.field_conditions_A == unificationAtomElement2.get_field_conditions_B() and
                    self.field_conditions_B == unificationAtomElement2.get_field_conditions_A() ) ) and
                 self.externalAttribute == unificationAtomElement2.get_external_attribute_list() )



class UnificationProtocol(object):
    """
    Class to represent
    """

    
    def __init__(self, description, BianaDatabaseVersion, id=None):
        """
        """

        self.unification_atom_elements = []

        self.description = description
        self.databaseVersion = BianaDatabaseVersion

        self.id = id  # Only has an ID if it has been inserted to the database

        self.use_databases = set()


    def add_database(self, externalDatabaseID):

        self.use_databases.add(externalDatabaseID)

    def get_database_ids(self):
        """
        Returns a set with the database ids it contains
        """

        return self.use_databases


    def add_unification_atom_elements(self, unificationAtomElement):
        """
        Adds a crossreference condition
        """
        self.unification_atom_elements.append(unificationAtomElement)
        self.use_databases.add(unificationAtomElement.externalDatabaseID_A)
        self.use_databases.add(unificationAtomElement.externalDatabaseID_B)

    def get_unification_atom_elements(self):
        return self.unification_atom_elements


    def get_description(self):
        return self.description
            

    def get_id(self):
        return self.id


    def get_summary_str_list(self):
        """
        """
        






