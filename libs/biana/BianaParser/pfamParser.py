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
File        : pfamParser.py
Author      : Javier Garcia
Creation    : 19 December 2007
Contents    : fills up tables in database biana with information from PFAM
Called from : 
=======================================================================================================

"""


## STEP 1: IMPORT NECESSARY MODULES

from bianaParser import *


class PFAMParser(BianaParser):
    """
    PFAM Parser Class
    """
    
    name = "pfam"
    description = "PFAM Parser. Collection of multiple sequence alignments and hidden markov models covering many common protein domains and families"
    external_entity_definition = "A external entity represents a pfam domain"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "Collection of multiple sequence alignments and hidden markov models covering many common protein domains and families",
                             default_script_name = "pfamParser.py",
                             default_script_description = PFAMParser.description,
                             additional_optional_arguments = [("include-sequence-mapping",True,"Includes sequence information and maps pfam to them"),
                                                              ("pfamA-file-name=","Pfam-A.full.gz","Name of the pfam A file. If None, it is not  used"),
                                                              ("pfamB-file-name=","Pfam-B.gz","Name of the pfam B file. If None, it is not used"),
                                                              ("pfamSeq-file-name=","pfamseq.gz","Name of the file with pfam sequences (only used if \"include-sequence-mapping\" option is selected")])
        self.default_eE_attribute = "pfam"
	#self.is_promiscuous = True


    def parse_database(self):
        """
        Method that implements the specific operations of PFAM parser
        """

        # Add the possibility to transfer pfam id using uniprotaccession
        self.biana_access._add_transfer_attribute( externalDatabaseID = self.database.get_id(), 
                                                   key_attribute = "uniprotaccession",
                                                   transfer_attribute="pfam" )

        self.biana_access._add_transfer_attribute( externalDatabaseID = self.database.get_id(), 
                                                   key_attribute = "uniprotentry",
                                                   transfer_attribute="pfam" )

        if not self.input_file.endswith(os.sep):
            self.input_file += os.sep

        self.uniprot_acc_seq_md5_dict = {} #! this dictionary is possibly causing memory problems

        self.insert_mapping = None

        # First, if "include_sequence_mapping" option is selected, gets the exact information of sequences used and inserts them to the database
        
        if self.arguments_dic["include-sequence-mapping"]:

            sequences_parsing_time = time.time()

            self.seq_database = ExternalDatabase( databaseName = self.sourcedb_name+"_sequences",
                                                  databaseVersion = self.sourcedb_version,
                                                  databaseFile = self.input_file.split(os.sep)[-1],
                                                  databaseDescription = self.database_description,
                                                  defaultExternalEntityAttribute = "proteinSequence")

            self.biana_access.insert_new_external_database( externalDatabase = self.seq_database )

            self.insert_mapping = 1            
            
            pfamseq_file = self.input_file+self.arguments_dic["pfamSeq-file-name"]

            if pfamseq_file.endswith(".gz"):
                pfamseq_input_file_fd = gzip.open(pfamseq_file,'r')
            else:
                pfamseq_input_file_fd = file(pfamseq_file, 'r')

            sequence = []
            protein_title_line = None

            protein_number = 0
            
            #self.title_regex = re.compile("^(\w+)\.(\d+)")

            uniprot_acc = None
            uniprot_acc_version = None
            uniprot_entry = None
            sequence = None

            self.title_regex = re.compile("^(\w+)\.(\d+)\s+(\S+)")

            for line in pfamseq_input_file_fd:

                if line[0]==">":
                    protein_number += 1

                    if self.time_control:
                        if protein_number%20000==0:
                            sys.stderr.write("%s proteins in %s seconds\n" %(protein_number,time.time()-self.initial_time))

                    if sequence is not None:
                        self.parse_pfam_seq_record(header_line = protein_title_line, sequence = "".join(sequence))

                    protein_title_line = line[1:]
                    sequence = []
                else:
                    sequence.append(line.strip())
        
            if len(sequence)>0:
                self.parse_pfam_seq_record(header_line = protein_title_line, sequence = "".join(sequence))

            # Updates the information that this external database has inserted
            self.seq_database.set_parsing_time( int(time.time() - sequences_parsing_time) )
            self.biana_access.update_external_database_external_entity_attributes( self.seq_database )

        pfamseq_input_file_fd.close()


        if self.arguments_dic["pfamA-file-name"].lower() != "none":
            self.parse_pfam_file(pfamType='A', pfamFile= self.input_file+self.arguments_dic["pfamA-file-name"])

        
        if self.arguments_dic["pfamB-file-name"].lower() != "none":
            self.parse_pfam_file(pfamType='B', pfamFile= self.input_file+self.arguments_dic["pfamB-file-name"])


    def parse_pfam_seq_record(self, header_line, sequence):

        pfamseq_object = ExternalEntity( source_database = self.seq_database, type="protein" )

        m = self.title_regex.match(header_line)

        if m:
            uniprot_acc = m.group(1)
            uniprot_acc_version = m.group(2)
            uniprot_entry = m.group(3)
        else:
            raise "ERROR in parsing line %s" %header_line

        pfamseq_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "uniprotentry", value = uniprot_entry,type = "cross-reference" ))
        pfamseq_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "uniprotaccession", value = uniprot_acc, version = uniprot_acc_version, type = "cross-reference"))

        sequenceObject = ProteinSequence(sequence.strip())

        pfamseq_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "proteinSequence", value = sequenceObject, type = "cross-reference"))

        self.uniprot_acc_seq_md5_dict["%s.%s" %(uniprot_acc, uniprot_acc_version)] = sequenceObject.get_sequence_MD5()

        self.biana_access.insert_new_external_entity( externalEntity = pfamseq_object )
            


    def parse_pfam_file(self, pfamType, pfamFile):

        # STOCKHOLM 1.0
        #=GF ID   14-3-3                  One word name for family
        #=GF AC   PF00244.11              Accession number in form PFxxxxx.version or PBxxxxxx.
        #=GF DE   14-3-3 protein          Short description of family
        
        #GF DR   PROSITE; PDOC00633;
        #GF DR   SMART; 14_3_3;
        #GF DR   PRINTS; PR00305;
        #=GF DR   SCOP; 1a4o; fa;
        #=GF DR   INTERPRO; IPR000308;
        #=GF DR   PFAMB; PB176422;
        #=GF DR   PRODOM; PD000197;
        #=GF DR   PFAMA; PF00297.13;

        #=GF SQ   560
        #=GS Q7M332_SHEEP/3-161      AC Q7M332.1
        #=GS 1433Z_BOVIN/3-236       DR PDB; 1a38 B; 3-228;

        if pfamFile.endswith(".gz"):
            pfam_input_file_fd = gzip.open(pfamFile,'r')
        else:
            pfam_input_file_fd = file(pfamFile, 'r')
            
        line_number=0

        # Create a new external entity object
        pfam_object = ExternalEntity( source_database = self.database, type="pattern")

        search_id = re.compile("^#=GF ID\s+(.+)$")
        search_description = re.compile("^#=GF DE\s+(.+)$")
        search_prosite = re.compile("^#=GF DR\s+PROSITE;\s+(\w+);")
        search_interpro = re.compile("^#=GF DR\s+INTERPRO;\s+(\w+);")
        search_scop = re.compile("^#=GF DR\s+SCOP;\s+(\w+);\s+(\w+);")
        search_prints = re.compile("^#=GF DR\s+PRINTS;\s+(\w+);")
        search_pfam_ac = re.compile("^#=GF AC\s+(.+)\.*(\d*)$")
        search_pfama = re.compile("^#=GF DR\s+PFAMA;\s+(.+)\.*(\d*);$")
        search_pfamb = re.compile("^#=GF DR\s+PFAMB;\s+(.+);$")
        search_prodom = re.compile("^#=GF DR\s+PRODOM;\s+(.+);$")
        search_end = re.compile("^\/\/$")


        #=GS Q28DR3_XENTR/4-241      AC Q28DR3.1
        search_uniprot = re.compile("^#=GS\s+([\w\_]+)\/(\d+\-\d+)\s+AC\s+(\w+)\.(\d+)")

        #=GS 1433T_HUMAN/3-236       DR PDB; 2btp A; 3-234;
        search_pdb_cross_ref = re.compile("^#=GS\s+[\w\_]+\/\d+\-\d+\s+DR\s+PDB;\s+(\w+)\s+(\w*);\s+(\d+)\-(\d+)")

        entries = 0
        
        for line in pfam_input_file_fd:

            pfamAC = None
            pfamDescription = None
            pfamType = pfamType
            pfamName = None

            line_number += 1

            try:
                line.strip()

                if( search_end.match(line) ):
                    # Insert current pfam_object
                    # Insert current pfam attribute properties

                    self.biana_access.insert_new_external_entity( externalEntity = pfam_object )

                    if self.time_control:
                        if entries%200==0:
                            sys.stderr.write("%s pfam entries in %s seconds\n" %(entries,time.time()-self.initial_time))
                    
                    pfam_object = ExternalEntity( source_database = self.database, type="pattern")
                    continue


                #ID
                m = search_id.match(line)
                if m:
                    entries += 1
                    pfamName = m.group(1)
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "name", value = pfamName, type = "unique"))
                    continue

                #DESCRIPTION
                m = search_description.match(line)
                if m:
                    pfamDescription = m.group(1)
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "description", value = pfamDescription))
                    continue
                
                #PRINTS
                m = search_prints.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "prints", value = m.group(1), type = "cross-reference"))


                #PRODOM
                m = search_prodom.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "prodom", value = m.group(1), type = "cross-reference"))

                #PFAMA
                m = search_pfama.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "pfam", value = m.group(1), type = "cross-reference"))
                
                #PFAMB
                m = search_pfamb.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "pfam", value = m.group(1), type = "cross-reference"))
                    continue
                    
                #INTERPRO
                m = search_interpro.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "interpro", value = m.group(1), type = "cross-reference"))
                    continue
                
                
                #AC
                m = search_pfam_ac.match(line)
                if m:
                    pfamAC = m.group(1)
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "pfam", value = pfamAC, type = "unique") )
                    continue

                #PROSITE
                m = search_prosite.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "prosite", value = m.group(1), type = "cross-reference"))
                    continue

                #UNIPROT
                m = search_uniprot.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "uniprotaccession", value = m.group(3), version = m.group(4), type = "cross-reference") )
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "uniprotentry", value = m.group(1), type = "cross-reference"))

                    if self.insert_mapping:
                        try:
                            sequenceMD5 = self.uniprot_acc_seq_md5_dict["%s.%s" %(m.group(3),m.group(4))]
                            pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "sequenceMap", value = sequenceMD5,
                                                                               additional_fields = {"seq_range": m.group(2)},
									       type="cross-reference" ) )
                        except:
                            sys.stderr.write("SequenceMD5 not found for %s.%s\n" %(m.group(3),m.group(4)))
                        
                                                                                               

                    continue

                #PDB
                m = search_pdb_cross_ref.match(line)
                if m:
                    pfam_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "pdb", value = m.group(1), type = "cross-reference",
                                                                       additional_fields = {"chain": m.group(2),
                                                                                            "pdb_range": "%s-%s" %(m.group(3),m.group(4))} ))
                    continue

            except:
                print line
                print traceback.print_exc()
                sys.stderr.write("Error in parsing line %s\n" %(line_number))
                raise ValueError("error in parsing line")


        if pfam_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = pfam_object )
        
        pfam_input_file_fd.close()
