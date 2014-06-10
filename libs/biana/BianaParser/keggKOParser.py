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

"""
File        : keggkoParser.py
Author      : Javier Garcia Garcia
Creation    : January 2008
Contents    : fills up tables in database biana with information from kegg ko database
Called from : 

=======================================================================================================

This file implements a program that fills up tables in database biana with information of kegg ko databases

"""

from bianaParser import *


class KeggKOParser(BianaParser):
    """

    """

    name = "kegg_kO"
    description = "This file implements a program that fills up tables in database biana with information of kegg KO Database"
    external_entity_definition = "A external entity represents a KEGG KO"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "KEGG KO database",
                             default_script_name = "keggKOParser.py",
                             default_script_description = KeggKOParser.description )
        self.default_eE_attribute = "keggCode"
        self.initialize_input_file_descriptor()

    def parse_database(self):
        """
        """

        self.initialize_input_file_descriptor()

        # General regex
        continue_field_regex = re.compile("^\s{3,}([^;]+);*$")
        field_regex = re.compile("^(\w+)\s+([^;]+);*$")
        pathway_regex = re.compile("PATH\:\s+(map|rn)(\d+)\s+(.+)$")

        space_regex = re.compile("\s+")
        parenthesis_regex = re.compile("\(.+\)")  # used to eliminate extra information in sequence
        

        entry_regex = re.compile("^ENTRY\s+(\w+)\s+KO")

        dblink_split_regex = re.compile("(\w+)\:")

        kegg_ko_object = None

        temp_value = []           # List used to store the information of those fields that can have more than a single line
        current_field = None
        
        for line in self.input_file_fd:

            m = entry_regex.match(line)

            if m:
                if kegg_ko_object is not None:
                    self.biana_access.insert_new_external_entity( externalEntity = kegg_ko_object )
                
                kegg_ko_object = ExternalEntityRelation( source_database = self.database, relation_type="cluster" )
                kegg_ko_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value=m.group(1), type = "unique" ) )
                continue


            new_field = field_regex.match(line)
            if new_field:
                if current_field == "DBLINK":
                    all_db_links = " ".join(temp_value)
                    list_db_links = [ x.strip() for x in dblink_split_regex.split(all_db_links) ]

                    for actual_position in xrange(len(list_db_links)):
                        if list_db_links[actual_position] == "EC":
                            [ kegg_ko_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "EC", value=x, type="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]
                            
                        elif list_db_links[actual_position] == "COG":
                            [ kegg_ko_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "COG", value=x,type="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]
                            
                        elif list_db_links[actual_position] == "GO":
                            [ kegg_ko_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "GO", value=x,type="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]
                            
                current_field = new_field.group(1)
                temp_value = [new_field.group(2)]
            else:
                cont_value = continue_field_regex.match(line)
                if cont_value:
                    temp_value.append(cont_value.group(1))

                
        # Insert the last one
        if kegg_ko_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = kegg_ko_object )
            


        
