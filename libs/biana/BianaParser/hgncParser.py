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
File        : hgnc2piana.py
Author      : Javier Garcia
Creation    : 14 November 2007
Contents    : fills up tables in database piana with information from HGNC
Called from : 
=======================================================================================================

"""


## STEP 1: IMPORT NECESSARY MODULES

from bianaParser import *


class HGNCParser(BianaParser):
    """
    HGNC Parser Class
    """

    name = "hgnc"
    description = "This file implements a program that fills up tables in BIANA database with information from HGNC"
    external_entity_definition = "A external entity represents a protein"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "HUGO Gene Nomenclature Committee",
                             default_script_name = "hgncParser.py",
                             default_script_description = HGNCParser.description,
                             additional_compulsory_arguments = [])
        self.default_eE_attribute = "hgnc"
        
    def parse_database(self):
        """
        Method that implements the specific operations of HGNC parser

        # Python generated dict
        0 :  HGNC ID
        1 :  Approved Symbol
        2 :  Approved Name
        3 :  Status
        4 :  Locus Type
        5 :  Previous Symbols
        6 :  Previous Names
        7 :  Aliases
        8 :  Name Aliases
        9 :  Chromosome
        10 :  Date Approved
        11 :  Date Modified
        12 :  Date Symbol Changed
        13 :  Date Name Changed
        14 :  Accession Numbers
        15 :  Enzyme IDs
        16 :  Entrez Gene ID
        17 :  Ensembl Gene ID
        18 :  Mouse Genome Database ID
        19 :  Specialist Database Links
        20 :  Specialist Database IDs
        21 :  Pubmed IDs
        22 :  RefSeq IDs
        23 :  Gene Family Name
        24 :  Record Type
        25 :  Primary IDs
        26 :  Secondary IDs
        27 :  CCDS IDs
        28 :  VEGA IDs
        29 :  Locus Specific Databases
        30 :  GDB ID (mapped data)
        31 :  Entrez Gene ID (mapped data supplied by NCBI)
        32 :  OMIM ID (mapped data supplied by NCBI)
        33 :  RefSeq (mapped data supplied by NCBI)
        34 :  UniProt ID (mapped data supplied by UniProt)
        35 :  Ensembl ID (mapped data supplied by Ensembl)
        36 :  UCSC ID (mapped data supplied by UCSC)
        37 :  Rat Genome Database ID (mapped data supplied by RGD)

        """
        

		
        # List of tables to lock. It is used to improve speed inserts, as the indices are not updated in each insert
        # Commented. Locking all tables for the moment
        # tables_to_lock = [PianaGlobals.crossReferenceSource_table,
        #                   PianaGlobals.crossReferences_table]


        # HGNC Fields are the following: (CD: Multiple values, comma delimited    QCD: Multiple Quoited values in a comma delimited list
        #  0: HGNC ID
        #  1: Approved Symbol (Oficial Gene Symbol)
        #  2: Approved Name (Oficial Gene Name)
        #  3: Status
        #  4: Locus Type
        #  5: Previous Symbols CD
        #  6: Previous Names QCD
        #  7: Aliases CD
        #  8: Name Aliases QCD         # EMPTY!!!!
        #  9: Chromosome
        # 10: Date Approved
        # 11: Date Modified
        # 12: Date Symbol changed      NOT EXISTS!!!!
        # 13: Date Name Changed
        # 14: Accession Numbers CD
        # 15: Enzyme ID CD
        # 16: Entrez Gene ID (Replaeced Locus Link)
        # 17: Ensembl Gene ID
        # 18: MGD ID
        # 19: Specialist Database Links (CD)
        # 20: Specialist Database IDs (CD)         NOT EXISTS!!!!
        # 21: Pubmed IDs (CD)
        # 22: RefSeq IDs (CD) Only One is selected!
        # 23: Gene Family Name (CD)
        # 24: Record Type
        # 25: Primary IDs
        # 26: Secondary IDs
        # 27: CCDS IDs
        # 28: VEGA IDs
        # 29: Locus Specific Databases
        # 30: GBD ID
        # 31: Entrez Gene ID
        # 32: OMIM ID
        # 33: RefSeq
        # 34: Uniprot ID
        # 35: EnsembL
        # 36: UCSC ID
        # 37: RGD ID
        

        self.initialize_input_file_descriptor()

        line_number=0
        header_columns = {}

        columns = 0

        for line in self.input_file_fd:

            line_number += 1
            

            # Read columns of header line into dictionary
            if line_number == 1:
                value_list = line.strip().split("\t")
                header_columns = dict([ (value_list[i], i) for i in xrange(len(value_list))])
                #sys.stderr.write("%s columns in header\n" %len(value_list))
                columns = len(value_list)
                                

            if line_number>1:
	
                try:
	            if line_number>1:

                        line.strip()
	            
	            	# Create a new external entity object
	            	hgnc_object = ExternalEntity( source_database = self.database, type="protein" )

                        # ADDING TAXID AS IT ONLY CONTAINS HUMAN GENES
                        hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier="taxid",
                                                                          value=9606,
                                                                          type="unique"))

	            	line_fields = line.split("\t")

                        
                        if len(line_fields) != columns:
                            sys.stderr.write("Incorrect fields number\n%s\n" %(line))

                        

	                column_index = header_columns["HGNC ID"]
                        column_value =  line_fields[column_index].strip()
                        if column_value.startswith("HGNC:"):
                            hgnc_id = column_value[5:]
                        else:
                            hgnc_id = column_value
	                # Take the values. Those that can be multiple values are stored in a list
                        hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "hgnc", 
                                                                          value = hgnc_id,
                                                                          type = "unique" ))
	                
	                column_index = header_columns["Approved Symbol"]
	                official_gene_symbol = line_fields[column_index].strip()
                        hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "geneSymbol", 
                                                                          value = official_gene_symbol,
                                                                          type = "unique" ))
	                
	                column_index = header_columns["Approved Name"]
	                official_gene_name = line_fields[column_index].strip()
	                # Oficial gene Name is entered as a description
                        hgnc_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "description",
                                                                           value = official_gene_name ))
	                
	                column_index = header_columns["Previous Symbols"]
	                previous_symbols = line_fields[column_index].strip()
	                if len(previous_symbols)>0:
                            previous_symbols = [ x.strip() for x in previous_symbols.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "geneSymbol", 
                                                                                value = x,
                                                                                type = "previous")) for x in previous_symbols ]

	                column_index = header_columns["Previous Names"]
	                previous_names = line_fields[column_index].strip()
	                if len(previous_names)>0:
                            previous_names = [ x.strip('"\s') for x in previous_names.split("\",\"") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "description", 
                                                                                value = x )) for x in previous_names ]

	                column_index = header_columns["Aliases"]
	                aliases_symbol = line_fields[column_index].strip()
	                if len(aliases_symbol)>0:
                            aliases_symbol = [ x.strip() for x in aliases_symbol.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "geneSymbol", 
                                                                                value = x,
                                                                                type = "alias" )) for x in aliases_symbol ]
	                	
	                column_index = header_columns["Accession Numbers"]
	                accession_numbers = line_fields[column_index].strip()
	                if len(accession_numbers)>0:
                            accession_numbers = [ x.strip() for x in accession_numbers.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "accessionNumber", 
                                                                                value = x,
                                                                                type = "cross-reference")) for x in accession_numbers ]


                        column_index = header_columns["Name Aliases"]
	                accession_numbers = line_fields[column_index].strip()
	                if len(accession_numbers)>0:
                            accession_numbers = [ x.strip() for x in accession_numbers.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "description", 
                                                                                value = x,
                                                                                type = "alias")) for x in accession_numbers ]

	                	
	                column_index = header_columns["Enzyme IDs"]
	                enzyme_IDs = line_fields[column_index].strip()
	                if len(enzyme_IDs)>0:
                            enzyme_IDs = [ x.strip() for x in enzyme_IDs.split(",") ]
                            new_enzyme_IDs = []
                            for id in enzyme_IDs:
                                m = re.match("\s*(.+\..+\..+\..+)", id)
                                if m:
                                    new_enzyme_IDs.append(m.group(1))
                            enzyme_IDs = new_enzyme_IDs

                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "EC", 
                                                                                value = x,
                                                                                type = "cross-reference" ) ) for x in enzyme_IDs ]

	                column_index = header_columns["Entrez Gene ID"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0:
                            geneIDs = [ x.strip() for x in column_value.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "geneID", 
                                                                                value = x,
                                                                                type = "cross-reference")) for x in geneIDs ]
	                	
	                column_index = header_columns["Mouse Genome Database ID"]
	                column_value = line_fields[column_index]
	                if len(column_value)>0:
                            MGD_IDs = [ x.lstrip("MGI:") for x in column_value.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "mgi", 
                                                                                value = x,
                                                                                type = "cross-reference")) for x in MGD_IDs ]
	                	
	                column_index = header_columns["RefSeq IDs"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0: 
                            refseqs = [ x.strip() for x in column_value.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "refseq", 
                                                                                value = x,
                                                                                type = "cross-reference")) for x in refseqs ]
	                	
	                column_index = header_columns["GDB ID (mapped data)"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0:
                            GDB_IDs = [ x.lstrip("GDB:") for x in column_value.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "gdb", 
                                                                                value = x,
                                                                                type = "cross-reference" )) for x in GDB_IDs ]
	                	
	                column_index = header_columns["Entrez Gene ID (mapped data supplied by NCBI)"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0:
                            mapped_geneIDs = [ x.strip() for x in column_value.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "geneID", 
                                                                                value = x,
                                                                                type = "cross-reference")) for x in mapped_geneIDs ]
	                	
	                column_index = header_columns["OMIM ID (mapped data supplied by NCBI)"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0:
                            omimIDs = [ x.strip() for x in column_value.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "mim",
                                                                                value = x, 
										type="cross-reference")) for x in omimIDs ]
	                	
	                column_index = header_columns["RefSeq (mapped data supplied by NCBI)"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0:
                            mapped_refseqs = [ x.strip() for x in column_value.split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "refseq", 
                                                                                value = x,
                                                                                type = "cross-reference")) for x in mapped_refseqs ]

	                column_index = header_columns["UniProt ID (mapped data supplied by UniProt)"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0:
                            uniprotIDs = [ x.strip() for x in column_value.strip().split(",") ]
                            [ hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "uniprotaccession", 
                                                                                value = x,
                                                                                type = "cross-reference")) for x in uniprotIDs ]

                        column_index = header_columns["Rat Genome Database ID (mapped data supplied by RGD)"]
	                column_value = line_fields[column_index].strip()
	                if len(column_value)>0:
                            RGD_IDs = [ x.lstrip("RGD:") for x in column_value.split(",") ]
                            for current_rgd_id in RGD_IDs:
                                if current_rgd_id.strip() != "":
                                    hgnc_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "rgd", 
                                                                                      value = current_rgd_id.strip(),
                                                                                      type = "cross-reference"))
	                	
	                # Save the object in the database            	
	                self.biana_access.insert_new_external_entity( externalEntity = hgnc_object )
	
	
                except:
                    traceback.print_exc()
                    sys.stderr.write("Error in parsing line %s\n" %(line_number))
                    raise Exception;


