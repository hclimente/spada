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

from bianaParser import *


class TaxonomyParser(BianaParser):
    """
    Taxonomy Parser Class
    """

    name = "taxonomy"
    description = "This program fills up tables in database biana related with taxonomy ontology"
    external_entity_definition = "A external entity represents a taxonomy of any type"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

         BianaParser.__init__(self, default_db_description = "NCBI Species information",
                              default_script_name = "taxonomyParser",
                              default_script_description = "This program fills up tables in database piana related to external DB 'taxonomy'")
         self.default_eE_attribute = "taxid"

    def parse_database(self):
        """
        Method that implements the specific operations of taxonomy parser
        """

        self.biana_access.add_valid_external_entity_attribute_type( name = "TaxID_category",
                                                                    data_type = "ENUM(\"class\",\"family\",\"forma\",\"genus\",\"infraclass\",\"infraorder\",\"kingdom\",\"no rank\",\"order\",\"parvorder\",\"phylum\",\"species\",\"species group\",\"species subgroup\",\"subclass\",\"subfamily\",\"subgenus\",\"subkingdom\",\"suborder\",\"subphylum\",\"subspecies\",\"subtribe\",\"superclass\",\"superfamily\",\"superkingdom\",\"superorder\",\"superphylum\",\"tribe\",\"varietas\")",
                                                                    category = "eE descriptive attribute")

        
        self.biana_access.add_valid_external_entity_attribute_type( name = "TaxID_name",
                                                                    data_type = { "fields": [("value","varchar(255)"),
                                                                                             ("taxid_name_type", "ENUM(\"acronym\",\"anamorph\",\"blast name\",\"common name\",\"equivalent name\",\"genbank acronym\",\"genbank anamorph\",\"genbank common name\",\"genbank synonym\",\"includes\",\"in-part\",\"misnomer\",\"misspelling\",\"scientific name\",\"synonym\",\"teleomorph\",\"authority\",\"unpublished name\")", True)] },
                                                                    category = "eE special attribute")


        # IMPORTANT: As we have added new types and attributes that are not in the default BIANA distribution, we must execute the follwing command:
        self.biana_access.refresh_database_information()                                                                        
                                                                    
        # Add the possibility to transfer taxonomy name and taxonomy category using taxID as a key
        self.biana_access._add_transfer_attribute( externalDatabaseID = self.database.get_id(), # A single taxonomy element can have multiple names
                                                   key_attribute = "taxID",
                                                   transfer_attribute="TaxID_name" )

        self.biana_access._add_transfer_attribute( externalDatabaseID = self.database.get_id(),   # Category is stored in a different attribute, as a single taxonomy can have multiple names
                                                   key_attribute = "taxID",
                                                   transfer_attribute = "TaxID_category" )

        nodes_dmp_file = None
        names_dmp_file = None

        if os.path.isdir(self.input_file):
            if( not self.input_file.endswith(os.sep) ):
                self.input_file += os.sep
            nodes_dmp_file = os.path.dirname(self.input_file) + os.sep + "nodes.dmp"
            names_dmp_file = os.path.dirname(self.input_file) + os.sep + "names.dmp"
        
        ontology = Ontology( source_database = self.database, linkedAttribute="taxid", name="taxonomy", descriptionAttribute="TaxID_name", levelAttribute="TaxID_category" )
        
        specific_identifiers_and_parent = {}

        names_fd = file(names_dmp_file,"r")

        names_comment_dict = {}

        print "Reading names"
        for line in names_fd:

            line_fields = line.split("|")

            if line_fields:

                tax_id = line_fields[0].strip()
                tax_name = line_fields[1].strip().replace('"','\"').replace("\\","\\\\")[0:254]
                tax_comment = line_fields[3].strip().replace('"','\"').replace("\\","\\\\")

                if tax_id != "" and tax_name != "":
                    	names_comment_dict.setdefault(tax_id, set()).add((tax_name, tax_comment))
                   
        names_fd.close()

        print "Reading nodes"
        nodes_fd = file(nodes_dmp_file,"r")
        for line in nodes_fd:

            line_fields = line.split("|")
            tax_id = line_fields[0].strip()
            parent_tax_id = line_fields[1].strip()
            type = line_fields[2].strip()
            

            externalEntity = ExternalEntity( source_database = self.database, type = "taxonomyElement" )

            externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID", value = tax_id, type = "unique" ) )
            externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "TaxID_category", value = type, type="unique" ) )
            
            for (tax_name, tax_comment) in names_comment_dict[tax_id]:
                externalEntity.add_attribute( ExternalEntityAttribute( attribute_identifier = "TaxID_name",
                                                                       value = tax_name,
                                                                       additional_fields = {"taxid_name_type": tax_comment},
								       type="unique") )

            self.biana_access.insert_new_external_entity( externalEntity )

            specific_identifiers_and_parent[tax_id] = (externalEntity.get_id(),parent_tax_id)

        nodes_fd.close()

        print "Done!"

        # Set the ontology hierarch and insert elements to ontology
        for current_tax_ID in specific_identifiers_and_parent:
            ontology.add_element( ontologyElementID = specific_identifiers_and_parent[current_tax_ID][0],
                                  isA = [specific_identifiers_and_parent[specific_identifiers_and_parent[current_tax_ID][1]][0]] )
            
        print "Inserting ontology to database"
        self.biana_access.insert_new_external_entity( ontology )



                                                           
