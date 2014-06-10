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

class CogParser(BianaParser):
    """
    COG Parser Class
    """

    name = "cog"
    description = "Clusters of Orthologous Groups of proteins (COGs)"
    external_entity_definition = "An element in a COG"
    external_entity_relations = "A COG"

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "COG database",
                             default_script_name = "cogParser.py",
                             default_script_description = CogParser.description,
                             additional_optional_arguments = [])
        self.default_eE_attribute = "cog"
	#self.is_promiscuous = True


    def parse_database(self):

        # FIRST: Check that all the files exist
        if os.path.isdir(self.input_file):
            self.input_path = self.input_file
        else:
            raise ValueError("You must specify a path instead of a file")

        files = ["myva","myva=gb","org.txt","fun.txt","whog"] #"pa" for the moment it is not necessary

        for current_file in files:
            if os.path.exists(self.input_path+os.sep+current_file) is False:
                raise ValueError("File %s is missing in %s" %(current_file, self.input_path))

        
        # Read correspondence letters to TaxID for the external entities

        species_file_fd = open(self.input_path+os.sep+"org.txt",'r')
        specie_taxid_dict = {}
        sp_taxid_regex = re.compile("\s*(\S+)\s+(\S+)\s+")

        for line in species_file_fd:
            m = sp_taxid_regex.match(line)
            if m:
                specie_taxid_dict[m.group(1).lower()] = m.group(2)

        species_file_fd.close()


        # Read the functional information

        function_dict = {}
        function_file_fd = open(self.input_path+os.sep+"fun.txt",'r')
        funct_regex = re.compile("\s*\[(\w+)\]\s+(.+)$")
        
        for line in function_file_fd:
            m =  funct_regex.match(line)
            if m:
                function_dict[m.group(1)] = m.group(2)

        function_file_fd.close()

        
        # Read the name and gi correspondence
        name_to_gi_dict = {}
        name2gi_file_fd = open(self.input_path+os.sep+"myva=gb",'r')
        name2gi_regex = re.compile("\s*(\S+)\s+(\S+)\s+$")

        for line in name2gi_file_fd:
            m = name2gi_regex.match(line)
            if m:
                if m.group(2) != "gi?":
                    name_to_gi_dict[m.group(1).lower()] = m.group(2)

        name2gi_file_fd.close()


        # Obtain, from the COGs file, to which specie belongs each protein
        # Obtain also the information for the COGs, description, functional_classification...
        whog_file_fd = open(self.input_path+os.sep+"whog",'r')
        name2species_dict = {}
        cogs_components_dict = {}
        cogs_funct_dict = {}
        cogs_description_dict = {}
        name2cogs_dict = {}
        current_cog = None

        
        new_cog_regex = re.compile("\s*\[(\w+)\]\s+(\w+)\s+(.+)$")
        assignment_regex = re.compile("\s*(\w{3})\:\s+(.+)$")

        for line in whog_file_fd:

            m = new_cog_regex.match(line)
            if m:
                cogs_description_dict[m.group(2)] = m.group(3)
                cogs_funct_dict[m.group(2)] = m.group(1)
                cogs_components_dict.setdefault(m.group(2),[])
                current_cog = m.group(2)
                continue

            m = assignment_regex.match(line)
            
            if m:
                components = m.group(2).split(" ")
                for current_component in components:
                    cogs_components_dict[current_cog].append(current_component)
                    name2cogs_dict.setdefault(current_component.lower(),[]).append(current_cog)
                    name2species_dict.setdefault(current_component.lower(),set()).add(m.group(1).lower())


        whog_file_fd.close()


        def create_and_insert_eE():
            eE_object = ExternalEntity( source_database = self.database, type="protein" )
            eE_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "proteinsequence", 
                                                                     value = ProteinSequence("".join(sequence)),
								     type = "cross-reference"))
            if name_to_gi_dict.has_key(protein_name.lower()):
                eE_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "GI", 
                                                                 value = name_to_gi_dict[protein_name.lower()],
								 type="cross-reference" ))

                
            if name2species_dict.has_key(protein_name.lower()):
                species = name2species_dict[protein_name.lower()]

                if len(species)>1:
                    print "Protein %s has more than a single specie assigned!" %protein_name

                for current_specie in species:
                    eE_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "taxID",
                                                                     value = specie_taxid_dict[current_specie.lower()],
								     type = "cross-reference") )

                for current_cog in name2cogs_dict[protein_name.lower()]:
                    eE_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "COG",
                                                                      value = current_cog,
								      type="cross-reference" ) )
                    for current_function in cogs_funct_dict[current_cog]:
                        eE_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "function",
                                                                          value = function_dict[current_function],
									  type="cross-reference" ) )

                    # HOW SHOULD THE NAME BE INSERTED??? In NCBI, they appear as Locus Name... Insert it as Ordered Locus Name?
                    eE_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "OrderedLocusName",
                                                                      value = protein_name,
								      type="cross-reference" ) )
                

                self.biana_access.insert_new_external_entity( externalEntity = eE_object )
            

        # Read the sequences and insert the external entities

        fasta_file_fd = open(self.input_path+os.sep+"myva",'r')
        sequence = []
        protein_name_regex = re.compile(">(.+)$")
        protein_name = None

        for line in fasta_file_fd:

            m = protein_name_regex.match(line)
            if m:
                if len(sequence)>0:
                    create_and_insert_eE()

                sequence = []
                protein_name = m.group(1)
            else:
                sequence.append(line.strip())

        fasta_file_fd.close()

        if len(sequence)>0:
            create_and_insert_eE()


        # FINALLY, IT WOULD BE POSSIBLE TO INSERT THE COGS AS A RELATION... IS IT NECESSARY? In principle, it is not necessary, as attribute networks can be generated...
        # For the moment, not done
