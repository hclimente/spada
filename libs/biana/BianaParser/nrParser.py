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
File        : nr2piana.py
Author      : Ramon Aragues & Javier Garcia
Modified    : Javier Garcia October 2007
Creation    : 10.2004
Contents    : fills up tables in database piana with information from nr
Called from : 
=======================================================================================================

This file implements a program that fills up tables in database piana with information from nr

nr can be downloaded from ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

Before running, taxonomy table of piana must be populated (use taxonomy2piana for that)
"""


## STEP 1: IMPORT NECESSARY MODULES

#from Bio import Fasta  # needed to read the nr file (which is in fasta format) # Not used because error in Windows and many other systems because of Martel package
from bianaParser import *


class NrParser(BianaParser):
    """
    NCBI nr Parser Class
    """

    name = "nr"
    description = "This file implements a program that fills up tables in BIANA database with information from nr"
    external_entity_definition = "A external entity represents a protein"
    external_entity_relations = ""
    
    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "Non-Redundant NCBI database",
                             default_script_name = "nrParser.py",
                             default_script_description = NrParser.description)

        self.default_eE_attribute = "proteinSequence"
        self.initialize_input_file_descriptor()


    def parse_database(self):
        """
        Method that implements the specific operations of nr parser
        """

        # Read specific parameters

        verbose_all = 0

        self.initialize_input_file_descriptor()
        

        # Initialize protein number (to the output printing...)
        self.protein_number=0         # Counter of total of protein lines in file
        self.total_entries=0          # Counter of total of sequence-taxid entries in file
        self.total_inconsistencies=0  # counter of inconsistencies found (inconsistencies with identifiers cross-references)

        #nr_parser = Fasta.RecordParser()
        #nr_input_file = self.input_file_fd
        #nr_iterator = Fasta.Iterator(nr_input_file, nr_parser)


        #nr_record = nr_iterator.next()

        self.initialize_input_file_descriptor()


        if self.verbose:
            print "Processing file"

        # get dictionary for gis and tax ids
        if self.verbose:
            print "Reading gi/taxonomy information..."

        dict_name_tax = self.biana_access.get_taxonomy_names_taxID_dict()


        # Dictionary to save inconsistencies (when a proteinPiana does not exists but exists any external Identifier...)
        inconsistencies = {}

        # while record read is not None, parse the record and insert data into piana
        
        sequence = []
        protein_title_line = None

        for line in self.input_file_fd:

            if line[0]==">":
                if len(sequence)>0:
                    self.parse_nr_record(header_line = protein_title_line, sequence = "".join(sequence))
                protein_title_line = line[1:]
                sequence = []
            else:
                sequence.append(line.strip())
        
        if len(sequence)>0:
            self.parse_nr_record(header_line = protein_title_line, sequence = "".join(sequence))


    def parse_nr_record(self, header_line, sequence):

            self.protein_number += 1
            self.add_to_log("Number of proteins in input file")


            if self.time_control:
                if self.protein_number%20000==0:
                    sys.stderr.write("%s proteins\t%s seq-taxid entries\tin %s seconds\n" %(self.protein_number,self.total_entries,time.time()-self.initial_time))

            protein_title_line = header_line
            protein_sequence = sequence
            #protein_title_line = nr_record.title
            #protein_sequence = nr_record.sequence.strip()

            """
            Now, we have a proteinPiana that will be asigned to all other codes that we find in the title line
            """
            # title can look like this: gi|1346670|sp|P04175|NCPR_PIG NADPH--cytochrome P450 reductase (CPR) (P450R)
            #                           gi|5832586|dbj|BAA84019.1| maturase [Eucharis grandiflora]
            #                           gi|223154|prf||0601198A polymerase beta,RNA
            #                           gi|17538019|ref|NP_496216.1| surfeit 5 (18.1 kD) (2K591) [Caenorhabditis elegans]$gi|2497033|sp|Q23679|YWM3_CAEEL Hypothetical protein ZK970.3 in chromosome II$gi|3881897|emb|CAA88887.1| Hypothetical protein ZK970.3 [Caenorhabditis elegans]$gi|7511441|pir||T28127 hypothetical protein ZK970.3 - Caenorhabditis elegans
            #                           
            #    (character '>' from input file has been removed)
            #    (where character $ is actually octal character 001)
            #
            # in nr (as opposed to genpept) a title can have several gi identifiers for the same sequence
            # however, nr is a mess, so it is easier to get species and so on from genpept (that is why we parse it first)


            # So, given the messy format of nr, first of all split by octal character 001 to see if there is more than one gi
            title_entries = protein_title_line.split("\001")

            # for each of the title entries (normally just one) split it, retrieve information and insert codes, using same sequence for all)
            for title_entrie in title_entries:


                nr_object = ExternalEntity( source_database = self.database, type="protein" )

                nr_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "proteinsequence", 
                                                                 value = ProteinSequence(protein_sequence),
								 type = "unique" ))

                self.total_entries += 1

                gi = None

                if self.verbose:
                    sys.stderr.write("One entry of the title is: %s\n" %(title_entrie))

                title_atoms = title_entrie.split('|')

                # title_atom[1] (should) is always the gi code
                #if( title_atoms[0] == "gi" ):
                #    gi_id = int(title_atoms[1])

                #
                #  title_atom[2] can be (exhaustive list as of 10.2004):
                #
                #  - gb: for EMBL accessions
                #  - emb: for EMBL protein id
                #  - pir: for PIR accessions
                #  - sp: for swissprot
                #  - dbj: dna japan
                #  - prf: ???
                #  - ref: ref_seq
                #  - pdb: pdb code
                #  - tpg: ???
                #  - tpe: ???

                ### UPDATED: February 2007:

                #  - gb: for GenBank acession Number
                #  - emb: for EMBL accession Number
                #  - dbj: for DDBJ, DNA Database of Japan
                #  - pir: for NBRF PIR
                #  - prf: for Protein Research Foundation
                #  - sp:  for Swiss-Prot Entries
                #  - pdb: for Brookhaven Protein Data Bank
                #  - pat: for Patents
                #  - bbs: for GenInfo Backbone Id
                #  - gnl: for General database Identifier
                #  - ref: for NCBI Reference Sequence
                #  - lcl: for Local Sequence Identifier

                # title_atom[3] is the value (a protein id) indicated in title_atom[2]

                # In theory, only one has to exist... (because of the EMBL/GenBank/DDBJ nucleotide sequence database...
                #gb_accession = None
                #embl_accession = None
                #dbj_accession = None

                accessionNumber = None
                uniprot_accession = None  
                uniprot_entry = None
                ref_seq = None
                pdb_code = None
                pdb_chain = None
                pir = None
                description = None
                species_name = None

                # Read all the interesting cross-references:
                for i in xrange(len(title_atoms)):
                    if title_atoms[i] == "gi":
                        if gi is None:
                            gi = int(title_atoms[i+1])
                            nr_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "gi",
                                                                            value = gi,
                                                                            type = "unique" ))
                        else:
                            sys.stderr.write(" Error parsing... how an entry can have several gi values? \n")
                    elif title_atoms[i] == "gb" or title_atoms[i] == "emb" or title_atoms[i] == "dbj":
                        if accessionNumber is None:
                            accession = re.search("(\w+)\.(\d+)",title_atoms[i+1])
                            if accession:
                                accessionNumber = accession.group(1)
                                accessionNumber_version = accession.group(2)
                                nr_object.add_attribute(ExternalEntityAttribute(attribute_identifier="accessionNumber", 
                                                                                value = accessionNumber,
                                                                                version = accessionNumber_version,
                                                                                type = "unique" ))
                        else:
                            sys.stderr.write(" Error parsing... how an entry can have several accesionNumber.version values? \n")
                    elif title_atoms[i] == "pir":
                        if pir is None:
                            pirRe = re.search("(\w+)\.*",title_atoms[i+2])
                            if pirRe:
                                pir = pirRe.group(1)
                                nr_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "pir", 
                                                                                value  = pir,
                                                                                type = "cross-reference" ))
                        else:
                            sys.stderr.write(" Error parsing... How an entry can have several pir values? \n")

                    elif title_atoms[i] == "sp":
                        if uniprot_accession is None:
                            uniprot_accession = title_atoms[i+1][0:6]
                        else:
                            sys.stderr.write(" There are more than one uniprot accession for the same gi... \n")
                        if uniprot_entry is None:
                            # "Chapuza" because there is one entry that has no uniprotEntry name...and the parser does not work correctly
                            # This entry is: gi|20148613|gb|AAM10197.1| similar to high affinity sulfate transporter 2|sp|P53392 [Arabidopsis thaliana]
                            if len(title_atoms) > (i+2):
                                uniprot_entry = (title_atoms[i+2].split())[0].strip()
                                nr_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "uniprotentry", 
                                                                                 value = uniprot_entry,
                                                                                 type = "cross-reference"))
                        else:
                            sys.stderr.write(" There are more than one uniprot entry for the same gi... \n")
                    elif title_atoms[i] == "ref":
                        if ref_seq is None:
                            ref_seq = title_atoms[i+1].strip().split(".")
                            nr_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "refseq", 
                                                                            value = ref_seq[0],
                                                                            version = ref_seq[1],
                                                                            type = "cross-reference"))
                        else:
                            sys.stderr.write(" There are more than one uniprot RefSeq for the same gi... \n")
                    elif title_atoms[i] == "pdb":
                        if pdb_code is None:
                            pdb_code = title_atoms[i+1].strip()
                            pdb_chain = title_atoms[i+2][0].strip()
                            nr_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "pdb", 
                                                                            value = pdb_code,
                                                                            additional_fields = { "chain": pdb_chain },
                                                                            type = "cross-reference"))
                        else:
                            sys.stderr.write(" There are more than one PDB for the same gi... \n")

                # Take the description
#                has_description = re.search("\|([^\|]*)\[",title_entrie)
#                if has_description:
#                    description = has_description.group(1)
#                    if description != "":
#                        description = description.replace('"'," ").replace("\n", " ").replace("'"," ").replace('\\', " ").strip(),
#                        nr_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "description", 
#                                                                        value = description ))



                # TO CHECK:
                # Take the specie:
                #has_specie = re.search("\[(.*)\]$",title_entrie)
                #if has_specie:
                #    species_name = has_specie.group(1).replace('"','').replace("(","").strip()


                # Now, get the tax_ID associated with the specie
                # Process species name...

                starts = []
                ends = []

                for x in xrange(len(title_entrie)):
                    if title_entrie[x]=='[':
                        starts.append(x)
                    elif title_entrie[x]==']':
                        ends.append(x)

                starts.reverse()
                ends.reverse()

                try:
                    limit = starts[0]

                    for x in xrange(1,len(ends)):
                        if ends[x]>starts[x-1]:
                            limit = starts[x]
                            continue
                        limit = starts[x-1]
                        break

                    species_name = title_entrie[limit+1:-1]

                    try:

                        tax_id_value = dict_name_tax[species_name.lower()]
                        nr_object.add_attribute( ExternalEntityAttribute( attribute_identifier="taxid", value=tax_id_value, type = "unique") )

                    except:
                        # If not found, try to split by "," or by "="
                        try:
                            tax_id_value = dict_name_tax[species_name.split(",")[0].lower()]
                            nr_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID", value = tax_id_value, type = "unique"))
                        except:

                            try:
                                tax_id_value = dict_name_tax[species_name.split("=")[0].lower()]
                                nr_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID", value = tax_id_value, type = "unique"))

                            except:
                                if self.verbose:
                                    sys.stderr.write("%s not found\n" %(species_name))


                    nr_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "description",
                                                                     value = title_entrie[:limit].strip() ) )

                except:
                    nr_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "description",
                                                                     value = title_entrie.strip()))


                self.biana_access.insert_new_external_entity( externalEntity = nr_object )


            # reading next record
            if self.verbose:
                sys.stderr.write( "\n---------------reading next record---------------\n")
            #nr_record = nr_iterator.next()

        # END OF while nr_record is not None


