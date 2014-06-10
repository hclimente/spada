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
File        : ipi2piana.py
Author      : Javier Garcia
Creation    : November 2007
Contents    : fills up tables in database piana with information from IPI
Called from : 
=======================================================================================================
"""

import re
from bianaParser import *


class IPIParser(BianaParser):
    """
    IPI Parser Class
    """

    name = "ipi"
    description = "Inserts information of IPI database into BIANA"
    external_entity_definition = "External entities are proteins"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "IPI. International Protein Index",
                             default_script_name = "ipi2piana.py",
                             default_script_description = IPIParser.description,
                             additional_compulsory_arguments = [])
        self.default_eE_attribute = "ipi"


    def parse_database(self):

        self.initialize_input_file_descriptor()

        if self.input_file_fd is not None:
            self.parse_file()
        else:  # is a directory
            dirname = os.path.dirname(self.input_file+os.sep)+os.sep
            files = os.listdir(self.input_file)
            for current_file in files:
                if re.search("ipi.\w+\.fasta",current_file):
                    if ( current_file.endswith(".gz") ):
                        self.input_file_fd = gzip.open(dirname+current_file,'r')
                    else:
                        self.input_file_fd = file(dirname+current_file, 'r')
                    self.parse_file()


    def parse_file(self):
        """
        Method that implements the specific operations of HGNC parser
        """

        # Example:

        #>IPI:IPI00000001.2|SWISS-PROT:O95793-1|TREMBL:Q59F99|ENSEMBL:ENSP00000360922;ENSP00000379466|REFSEQ:NP_059347|H-INV:HIT000329496|VEGA:OTTHUMP00000031233 Tax
        #_Id=9606 Gene_Symbol=STAU1 Isoform Long of Double-stranded RNA-binding protein Staufen homolog 1
        #MSQVQVQVQNPSAALSGSQILNKNQSLLSQPLMSIPSTTSSLPSENAGRPIQNSALPSAS
        #ITSTSAAAESITPTVELNALCMKLGKKPMYKPVDPYSRMQSTYNYNMRGGAYPPRYFYPF


        line_number = 0
        ipi_object = None
        ipi_object_number = 0
        actual_sequence = []

        for line in self.input_file_fd:

            line_number += 1

            line.strip()

            field_search_re = re.compile("([\w\-]+)\:(\S+)")
            tax_id_regex = re.compile("Tax_Id=(\d+)")
            gene_symbol_regex = re.compile("Gene_Symbol=(\S+)\s+(.*)")

            if line[0]=='>':


                # Insert the last entry to the database
                if ipi_object is not None:
                    ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="proteinSequence", 
                                                                     value = ProteinSequence("".join(actual_sequence)),
								     type = "cross-reference"))

                    self.biana_access.insert_new_external_entity( externalEntity = ipi_object )

                # Start new entry
                ipi_object = ExternalEntity( source_database = self.database, type="protein" )
                ipi_object_number += 1

                if self.time_control:
                    if ipi_object_number%20000==0:
                        sys.stderr.write("%s entries done in %s seconds\n" %(ipi_object_number,time.time()-self.initial_time))

                actual_sequence = []
                line_fields = line.lstrip(">").split("|")

                for actual_field in line_fields:
                    
                    search = field_search_re.search(actual_field)
                    
                    if search:
                        identifier_type = search.group(1)
                        values = search.group(2).split(";")

                        if( identifier_type == "IPI" ):
                            for actual_value in values:
                                #if actual_value.startswith("IPI"):
                                #    actual_value = actual_value[3:]
                                ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "ipi", 
                                                                                 value = actual_value,
                                                                                 type = "unique" ))
                    
                        elif( identifier_type == "ENSEMBL" ):
                            for actual_value in values:
                                ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="ensembl", 
                                                                                 value = actual_value,
                                                                                 type = "cross-reference" ))

                        elif( identifier_type == "REFSEQ" ):
                            for actual_value in values:
                                ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="refseq", 
                                                                                 value = actual_value,
                                                                                 type = "cross-reference" ))

                        elif( identifier_type == "TREMBL" ):
                            for actual_value in values:
                                ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="uniprotaccession", 
                                                                                 value = actual_value,
                                                                                 type = "cross-reference" ))

                        elif( identifier_type == "SWISS-PROT" ):
                            for actual_value in values:
                                ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="UniprotAccession", 
                                                                                 value = actual_value[0:6],
                                                                                 type = "cross-reference" ))

                        elif( identifier_type == "TAIR" ):
                            for actual_value in values:
                                ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="tair", 
                                                                                 value = actual_value,
                                                                                 type = "cross-reference" ))

                search = tax_id_regex.search(line)
                if search:
                    ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="taxID",
                                                                     value = search.group(1),
								     type = "cross-reference"))

                search = gene_symbol_regex.search(line)
                if search:
                    ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="geneSymbol",
                                                                     value = actual_value,
                                                                     type = "cross-reference" ))

                    search2 = re.search("[Emb|Gb]\|(\S+)",search.group(2))
                    if search2:
                        ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "accessionNumber",
                                                                         value = search2.group(1),
                                                                         type = "cross-reference" ))
                    else:
                        ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="description",
                                                                         value = search.group(2) ))

            else:
                # Sequence line
                actual_sequence.append(line.strip())

        # Insert the last entry
        if ipi_object is not None:
            ipi_object.add_attribute(ExternalEntityAttribute(attribute_identifier="proteinSequence", 
                                                             value = ProteinSequence("".join(actual_sequence)),
							     type = "cross-reference"))
            self.biana_access.insert_new_external_entity( externalEntity = ipi_object )
