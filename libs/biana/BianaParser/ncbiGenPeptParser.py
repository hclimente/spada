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


class NcbiGenPeptParserGenpeptParser(BianaParser):
    """
    Genpept Parser class
    """

    name = "ncbi_genpept"
    description = "Inserts NCBI Genpet database into BIANA"
    external_entity_definition = "A external entity represents a protein"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "NCBI genbank",
                             default_script_name = "ncbiGenPeptParser.py",
                             default_script_description = NcbiGenPeptParserGenpeptParser.description,
                             additional_optional_arguments = [])
        self.default_eE_attribute = "proteinSequence"
        self.initialize_input_file_descriptor()

    


    def parse_database(self):

        def insert_fsa_aa_files(arg,dirname,names):
            for name in names:
                if name.endswith('.fsa_aa'):
                    file_fd = open(os.path.join(dirname,name),'r')
                elif name.endswith('.fsa_aa.gz'):
                    file_fd = gzip.open(os.path.join(dirname,name),'r')
                else:
                    continue
                #try:
                self.parse_fasta_file(file_fd)
                #except:
                #    sys.stderr.write("Error parsing file %s\n" %name)
                #    traceback.print_exc()
                file_fd.close()
                

        self.protein_number=0

        self.dict_name_tax = self.biana_access.get_taxonomy_names_taxID_dict()

        if( len(self.dict_name_tax) == 0 ):
            raise ValueError("Taxonomies have not been inserted. You must previously insert information about taxonomy")

        # All these tags will be considered to be pointing to id type Accession
        self.accepted_coll_accessions = { "gb": None,
                                          "dbj": None,
                                          "emb": None}

        # Run all the path to insert all the hssp files of the path
        os.path.walk(self.input_file,insert_fsa_aa_files,None)



    def parse_fasta_file(self, fasta_file_fd):

        # Define the parser information to print in the help:

        # STEP 3: THE PARSER ITSELF. IT MUST USE THE METHODS OF PianaDBaccess

        

        if self.verbose:
            sys.stderr.write("Loading taxonomies\n")

        
        if self.verbose:
            sys.stderr.write( "Processing file\n")

        
        sequence = []
        protein_title_line = None

        for line in fasta_file_fd:

            if line[0]==">":
                
                if len(sequence)>0:
                    self.parse_genpept_entry(header_line = protein_title_line, sequence = "".join(sequence))

                protein_title_line = line[1:]
                sequence = []
            else:
                sequence.append(line.strip())
        
        if len(sequence)>0:
            self.parse_genpept_entry(header_line = protein_title_line, sequence = "".join(sequence))


    def parse_genpept_entry(self, header_line, sequence):

            genpept_object = ExternalEntity( source_database = self.database, type="protein" )

            self.protein_number += 1

            self.add_to_log("Number of proteins processed")

            if self.verbose:
                sys.stderr.write( "===================================\n")
                sys.stderr.write( "            NEW ENTRY\n")
                sys.stderr.write( "===================================\n")


            if self.time_control:
                if self.protein_number%20000==0:
                    sys.stderr.write("%s proteins done in %s seconds\n" %(self.protein_number,time.time()-self.initial_time))

            protein_title_line = header_line
            protein_sequence = sequence
            
            if protein_sequence:
                genpept_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "proteinsequence", 
                                                                       value = ProteinSequence(protein_sequence),
								       type = "unique") )

            if self.verbose:
                sys.stderr.write("title line is : %s\n" %(protein_title_line)) 

            # title looks like this: gi|gi_id|sourceDB|acc| protein name [species]
            #                        - gi|4105707|gb|AAD02507.1| carbamate kinase [Trichomonas vaginalis]
            #                        - gi|4105707|gb|AAD02507.1| carbamate kinase [Trichomonas vaginalis]
            #                        - gi|4105707|gb|AAD02507.1| carbamate kinase
            #                        - gi|4105707|dbj|BBAS02507.1| carbamate kinase 4243
            #                        
            #    (character '>' from input file has been removed)

            title_atoms = protein_title_line.split('|')

            gi_id = int(title_atoms[1])

            accessionNum = None
            accessionVersion = None
            
            if self.accepted_coll_accessions.has_key(title_atoms[2]):
                # we consider accessions from sources in accepted_coll_accessions
                coll_acc = title_atoms[3].strip()

                # Split Accession Number and Version
                try:
                    splited = re.split("\.",coll_acc)
                    accessionNum = splited[0]
                    accessionVersion = splited[1]
                except:
                    accessionNum = coll_acc

            temp = []
            for x in xrange(4,len(title_atoms)):
                temp.append(title_atoms[x])
            temp_protein_name = "".join(temp)


            # Process species name...

            starts = []
            ends = []

            for x in xrange(len(temp_protein_name)):
                if temp_protein_name[x]=='[':
                    starts.append(x)
                elif temp_protein_name[x]==']':
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


                species_name = temp_protein_name[limit+1:-1]

                try:
                
                    tax_id_value = self.dict_name_tax[species_name.lower()]
                    genpept_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID",
                                                                       value = tax_id_value,
								       type = "unique"))
                except:
                    # If not found, try to split by "," or by "="
                    try:
                        tax_id_value = self.dict_name_tax[species_name.split(",")[0].lower()]
                        genpept_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID",
                                                                               value = tax_id_value,
									       type = "unique"))
                    except:

                        try:
                            tax_id_value = self.dict_name_tax[species_name.split("=")[0].lower()]
                            genpept_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID",
                                                                                   value = tax_id_value,
										   type = "unique"))
                        
                        except:
                            if self.verbose:
                                sys.stderr.write("%s not found\n" %(species_name))


                genpept_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "description", 
                                                                      value = temp_protein_name[:limit].strip() ) )

                self.add_to_log("Number of protein descriptions inserted")

            except:
                genpept_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "description", 
                                                                      value = temp_protein_name.strip() ) )
                self.add_to_log("Number of protein descriptions inserted")


            """
            Now, we have a proteinPiana for the (sequence, tax_id), either the one already existing before or a new one.

            Just add protein information to database tables.

            TO DO!! this would be a good place to check for consistency between uniprot and genbank, using embl accessions
            """

            if gi_id:
                genpept_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "gi", 
                                                                       value = gi_id,
                                                                       type = "unique" ) )
                self.add_to_log("Number of gi codes inserted")

            if accessionNum != "":
                if accessionVersion is None:
                    genpept_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "accessionNumber",
                                                                          value = accessionNum,
                                                                          type = "unique" ) )
                else:
                    genpept_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "accessionNumber",
                                                                          value = accessionNum,
                                                                          version = accessionVersion,
                                                                          type = "unique"))
                self.add_to_log("Number of accessionNumber inserted")



            # reading next record
            if self.verbose:
                sys.stderr.write( "reading next record\n")


            self.biana_access.insert_new_external_entity( externalEntity = genpept_object )


