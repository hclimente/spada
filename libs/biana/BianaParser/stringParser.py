from bianaParser import *
from biana.BianaDB import database
from sets import *
import os, fnmatch, re, sys

class STRINGParser(BianaParser):
    """
    STRING Parser Class
    """
    N_MAX_ENTRY_AT_ONCE = 5000 #10000 #20000 #!
    name = "string"
    description = "This program fills up tables in database biana related with STRING"
    external_entity_definition = "A protein in STRING with sequence"
    external_entity_relations = "A link given between two proteins in STRING with score"

    datatype_to_biana_type = { "SWISSPROT_DE": "name",
                               "SWISSPROT_GN": "geneSymbol",
                               "SWISSPROT_ID": "uniprotEntry",
                               "SWISSPROT_AC": "uniprotAccession",
                               "SWISSPROT_DR_EMBL": "AccessionNumber",
                               "SWISSPROT_DR_PIR": "PIR",
                               "SWISSPROT_DR_SGD": "SGD",
                               "GenomeReviews": "GenomeReviews",
                               "Uniprot/SWISSPROT": "uniprotEntry",
                               "Uniprot/SPTREMBL": "uniprotEntry",
                               "Uniprot/SPTREMBL": "uniprotEntry",
                               "Ensembl_Uniprot/SWISS": "uniprotEntry",
                               "EMBL": "AccessionNumber", 
                               "RefSeq": "RefSeq",
                               "RefSeq_peptide": "RefSeq",
                               "RefSeq_dna": "RefSeq",
                               "EntrezGene": "GeneID",
                               "Ensembl": "Ensembl",
                               "UniGene": "Unigene",
                               "Ensembl_SGD": "SGD",
                               "SGD": "SGD", 
                               "PDB": "PDB", 
                               "IPI": "IPI",
                               "HGNC": "HGNC",
                               "HPRD": "HPRD",
                               "PUBMED": "Pubmed",
                               "HUGO": "GeneSymbol",
                               "gene": "Tair",
                               "RGD": "RGD",
                               "MIM": "MIM",
                               "wormbase_transcript": "WormBaseSequenceName", 
                               "Ensembl_HUGO_Approved_Name": "Name", # unique
                               "Ensembl_HUGO_Aliases": "Name", # alias
                               "Ensembl_HUGO_RefSeq_IDs": "RefSeq",
                               "Ensembl_HUGO_Accession_Numbers": "AccessionNumber",
                               "Ensembl_HUGO_HGNC_ID": "HGNC",
                               "HUGE": "Huge", # no entry w/ this source
                               "GO": "GO", # no entry
                               "PRINTS": "PRINTS", # no entry
                               "primary": None, # telling that x is x (primary), no idea why it is there..
                               "Extra": None,
                               "protein_id": None,
                               "wormpep_id": None, 
                               "CCDS": None,
                               "aniseed": None, 
                               "ciona_int_jgi": None,
                               "cint_jgi_v2": None,
                               "cint_jgi_v1": None,
                               "cint_aniseed_v1": None,
                               "Ensembl_HUGO_Previous_Symbols": None,
                               "Uniprot/SPTREMBL_predicted": None, 
                               "RefSeq_peptide_predicted": None,
                               "RefSeq_dna_predicted": None,
                               "protein_id_predicted": None,
                               "prediction_SPTREMBL": None,
                               "EMBL_predicted": None,
                               "Genoscope_annotated_gene": None,
                               "gadfly_translation_cgid": None, 
                               "gadfly_transcript_cgid": None,
                               "gadfly_gene_cgid": None,
                               "wormbase_gene": None, 
                               "wormbase_locus": None, 
                               "flybase_polypeptide_id": None,
                               "flybase_gene_id": None,
                               "flybase_transcript_id": None,
                               "flybase_annotation_id": None,
                               "FlyBaseName_translations": None,
                               "FlyBaseName_transcript": None,
                               "FlyBaseName_gene": None,
                               "Celera_Trans": None,
                               "Celera_Pep": "", 
                               "Celera_Gene": "", 
                               "MEDLINE": None,
                               "OTTT": None,
                               "Genoscope_pred_gene": None,
                               "Genoscope_pred_transcript": None,
                               "Tiffin": None,
                               "miRNA_Registry": None,
                               "AFFY_Mouse430_2": None,
                               "AFFY_C_elegans": None, 
                               "AFFY_Canine": None,
                               "AFFY_MG_U74Cv2": None,
                               "AFFY_Zebrafish": None,
                               "AFFY_HG_U95E": None,
                               "AFFY_HG_U95D": None,
                               "AFFY_HG_U95C": None,
                               "AFFY_HG_U95B": None,
                               "AFFY_HG_U133A": None,
                               "AFFY_HG_U133B": None,
                               "AFFY_HG_U133_PLUS_2": None,
                               "AFFY_HG_U95Av2": None,
                               "AFFY_MG_U74Av2": None,
                               "AFFY_HuGeneFL": None,
                               "AFFY_HC_G110 and CCDS": None,
                               "AFFY_HC_G110": None,
                               "AFFY_U133_X3P": None,
                               "AFFY_RG_U34B": None,
                               "AFFY_RG_U34C": None,
                               "AFFY_RG_U34A": None,
                               "AFFY_DrosGenome1": None,
                               "AFFY_Drosophila_2": None,
                               "AFFY_HG_U133A_2": None,
                               "AFFY_Chicken": None,
                               "AFFY_Mouse430A_2": None,
                               "AFFY_HG_Focus": None,
                               "AFFY_MG_U74Bv2": None,
                               "AFFY_Rat230_2": None,
                               "ZFIN_xpat": None,
                               "ZFIN_ID": None,
                               "Ens_Hs_translation": None,
                               "Uniprot/Varsplic": None,
                               "MarkerSymbol": None,
                               "Xenopus_Jamboree": None,
                               "AgilentCGH": None,
                               "Ens_Hs_transcript": None,
                               "AgilentProbe": None,
                               "Anopheles_symbol": None,
                              }

    swissprot_re = re.compile("SWISSPROT_[A-Z]{2}[.]*")

    def __init__(self):
        # Start with the default values
        BianaParser.__init__(self, default_db_description = "Search Tool for the Retrieval of Interacting Genes/Proteins",
                             default_script_name = "stringParser.py",
                             default_script_description = "This program fills up tables in database biana related to STRING",
                             additional_compulsory_arguments = [],
                             additional_optional_arguments = [("use-existing-temp-alias-table", 0, "Uses previously created temporary alias table in the database")])
        self.default_eE_attribute = "string"
        self.use_existing_temp = self.arguments_dic["use-existing-temp-alias-table"]
        self.string_protein_object_number = 0
        self.setUnknownDB = Set()
        self.alias_temp_table = database.TableDB( table_name = "temp_string_aliases",
                                                  table_fields = [ database.FieldDB(field_name = "id", data_type = "varchar(53)"),      
                                                  database.FieldDB(field_name = "alias", data_type = "varchar(130)"), 
                                                  database.FieldDB(field_name = "source_list", data_type = "varchar(190)") ],
                                                  primary_key=("id, alias") ) #, indices=[("alias")] )
        return 

    def parse_database(self):
        """
        Method that implements the specific operations of string parser
        """
        ## Get STRING protein sequence, alias and links files
        (sequences_file_fd, aliases_file, links_file_fd) = self._get_data_file_names()
        string_protein_object = None
        ## dictionary storing processed string ids and corresponding external ids assigned to them in Biana database
        processed_string_ids_to_external_ids = {}
        ## insert alias information to database as a temporary table to be used later
        if self.use_existing_temp == 0:
            self._insert_protein_information_from_alias_file_to_database(aliases_file)
        nSequence = 0
        ## insert the data of proteins one by one fetching sequence and protein from sequence file and then searching for the information of protein in alias table in database
        lineLastRed = sequences_file_fd.readline()
        while lineLastRed:
            ## get identifier and sequence information
            lineLastRed, dictIdToSequence = self._readGivenNumberOfSequencesToDictionary(sequences_file_fd, lineLastRed, self.N_MAX_ENTRY_AT_ONCE)
            ## fetch alias information for the current protein 
            setSequenced = Set(dictIdToSequence.keys())
            nSequence += len(setSequenced)
            list_information_tuple = self._fetch_given_list_protein_information_from_database(setSequenced)
            setAliased = Set()
            if list_information_tuple != ():
                id_word_prev = list_information_tuple[0][0]
                list_information_tuple_inner = []
            for id_word, alias, source_list in list_information_tuple:
                setAliased.add(id_word)
                if id_word_prev == id_word:
                    list_information_tuple_inner.append((alias, source_list))
                else:
                    ## Start & insert new entry
                    sequence = dictIdToSequence[id_word_prev] 
                    processed_string_ids_to_external_ids[id_word_prev] = self._insert_new_string_protein_object_into_database(id_word_prev, sequence, list_information_tuple_inner)
                    self.check_time()
                    list_information_tuple_inner = [(alias, source_list)]
                id_word_prev = id_word
                    
            if list_information_tuple != ():
                id_word_last = list_information_tuple[-1][0]
                sequence = dictIdToSequence[id_word_last] 
                processed_string_ids_to_external_ids[id_word_last] = self._insert_new_string_protein_object_into_database(id_word_last, sequence, list_information_tuple_inner)
                self.check_time()
            setToBeProcessed = setSequenced - setAliased
            if self.verbose and len(setToBeProcessed):
                print "Sequence entries without alias: ", len(setToBeProcessed)
            for id_word in setToBeProcessed:
                sequence = dictIdToSequence[id_word] 
                processed_string_ids_to_external_ids[id_word] = self._insert_new_string_protein_object_into_database(id_word, sequence, [])
                self.check_time()

        self._close_file_descriptor(sequences_file_fd)

        if self.verbose:
            print "Number of entries with sequence: ", nSequence
            print "Number of entries with alias: ", self.biana_access.db.select_db_content( sql_query= "SELECT COUNT(DISTINCT(id)) FROM %s" % self.alias_temp_table.get_table_name(), answer_mode = "single" )

        if self.verbose:
            print "Unknown databases:", self.setUnknownDB

	relations_inserted_before = set()
        ## Insert relations
        string_relation_object = None
        string_relation_object_number = 0
        for line in links_file_fd:
            ## get interaction & score information between two entries
	    if line.startswith("#"):
		continue
	    if line.startswith("protein"):
		continue
            line_fields = line.strip().split()   
            if len(line_fields) < 10:
                if self.verbose:
                    print "Warning: Format inconsistency - missing scores", line
                continue
            id_word1 = line_fields[0]
            id_word2 = line_fields[1]
            neighborhood = line_fields[2]
            fusion = line_fields[3]
            cooccurence = line_fields[4]
            coexpression = line_fields[5]
            experimental = line_fields[6]
            database = line_fields[7]
            textmining = line_fields[8]
            score = line_fields[9]
            ## if sequence information for these proteins was not available they were not inserted
            ## in these cases insert their alias information to the database 
            if not processed_string_ids_to_external_ids.has_key(id_word1):
                list_information_tuple = self._fetch_given_protein_information_from_database(id_word1)
                if list_information_tuple == ():
                    if self.verbose:
                        print "Warning: id not found in temp alias table:", id_word1
                processed_string_ids_to_external_ids[id_word1] = self._insert_new_string_protein_object_into_database(id_word1, None, list_information_tuple)

            if not processed_string_ids_to_external_ids.has_key(id_word2):
                list_information_tuple = self._fetch_given_protein_information_from_database(id_word2)
                if list_information_tuple == ():
                    if self.verbose:
                        print "Warning: id not found in temp alias table:", id_word2
                processed_string_ids_to_external_ids[id_word2] = self._insert_new_string_protein_object_into_database(id_word2, None, list_information_tuple)

	    if (id_word2, id_word1) in relations_inserted_before:
		relations_inserted_before.remove((id_word2, id_word1))
		continue
	    else:
		relations_inserted_before.add((id_word1, id_word2))
            ## Start & insert new entry relation
            string_relation_object = ExternalEntityRelation( source_database = self.database, relation_type="functional_association") #"interaction" )
            string_relation_object.add_participant( externalEntityID = processed_string_ids_to_external_ids[id_word1] )
            string_relation_object.add_participant( externalEntityID = processed_string_ids_to_external_ids[id_word2] )
            #string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore", 
            #    value = score, additional_fields = { "neighborhood": neighborhood, "fusion": fusion, "cooccurence": cooccurence, "coexpression": coexpression, "experimental": experimental, "db": database, "textmining": textmining })) # changed with below
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore", value = score))
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore_neighborhood", value = neighborhood))
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore_fusion", value = fusion))
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore_cooccurence", value = cooccurence))
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore_coexpression", value = coexpression))
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore_experimental", value = experimental))
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore_db", value = database))
            string_relation_object.add_attribute(ExternalEntityRelationAttribute( attribute_identifier = "STRINGScore_textmining", value = textmining))
            self.biana_access.insert_new_external_entity( externalEntity = string_relation_object )
            string_relation_object_number += 1
            if self.time_control:
                if string_relation_object_number%20000==0:
                    sys.stderr.write("%s relation entries done in %s seconds\n" %(string_relation_object_number,time.time()-self.initial_time))

        self._close_file_descriptor(links_file_fd)
        #print "Unknown databases:", self.setUnknownDB
        
        #self._remove_protein_information_from_database() #!

        return

    def check_time(self):
        self.string_protein_object_number += 1
        if self.time_control:
            if self.string_protein_object_number%20000==0:
                sys.stderr.write("%s entries done in %s seconds\n" %(self.string_protein_object_number,time.time()-self.initial_time))
        return

    def _readGivenNumberOfSequencesToDictionary(self, sequences_file_fd, lineLastRed, nSequence):
        i = 0
        dictIdToSequence = {}
        while lineLastRed and i< nSequence:
            lineLastRed, id_word, sequence = self._readNextSequenceInformationFromSequenceFile(sequences_file_fd, lineLastRed)
            i += 1
            dictIdToSequence[id_word] = sequence
        return lineLastRed, dictIdToSequence

    def _readNextSequenceInformationFromSequenceFile(self, sequences_file_fd, line):
        while line.startswith("#"):
            line = sequences_file_fd.readline()
            continue
        if not line.startswith(">"):
            if self.verbose:
                print "Warning: unexpected line in sequences file:", line
            return
        line_fields = line.strip().split()   
        id_word = line_fields[0][1:]
        line = sequences_file_fd.readline()
        if line.startswith(">"):
            if self.verbose:
                print "Warning: unexpected > in fasta file", line
            return
        sequence = ""
        while line and not line.startswith('>'):
            sequence += line.strip()
            line = sequences_file_fd.readline()
        return line, id_word, sequence

    def _insert_new_string_protein_object_into_database(self, id_word, sequence, list_information_tuple):

        index = id_word.find(".")
        id = id_word[index+1:]
        tax = id_word[:index]
        #print id_word, sequence, list_information_tuple 
        
        string_protein_object = ExternalEntity( source_database = self.database, type="protein" )

        value = self.biana_access._transform_attribute_value_data_type_to_biana_database_attribute_data_type( attribute_identifier="STRING", value=id )
        string_protein_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "STRING", value = id, type = "unique"))
        value = self.biana_access._transform_attribute_value_data_type_to_biana_database_attribute_data_type( attribute_identifier="taxID", value=tax )
        string_protein_object.add_attribute(ExternalEntityAttribute(attribute_identifier = "taxID", value = tax, type = "cross-reference")) 
        if sequence is not None:
            string_protein_object.add_attribute(ExternalEntityAttribute(attribute_identifier="proteinSequence", value = ProteinSequence("".join(sequence)), type = "cross-reference"))
        for alias, source_list_str in list_information_tuple:
            for source in source_list_str.split():
                search = self.swissprot_re.search(source)
                source_org = source
                if search:
                    source = source[search.start():]
                if not self.datatype_to_biana_type.has_key(source):
                    #sys.stderr.write("Warning: Unknown source db id - %s" % source+"\n")
                    self.setUnknownDB.add(source_org)
                    source_db = None
                else:
                    source_db = self.datatype_to_biana_type[source]
                if source_db is not None:
                    type = "cross-reference"
                    #if source_db.lower() == "ipi":
                    #    if alias.startswith("IPI"):
                    #        alias = alias[3:]
                    if source_db.lower() == "hgnc":
                        if alias.startswith("HGNC:"):
                            alias = alias[5:]
                    if source == "Ensembl_HUGO_Aliases":
                        type = "alias"
                    if source == "Ensembl_HUGO_Approved_Name":
                        type = "unique"
                    #if alias is None or len(alias) < 2:
                    #    if self.verbose:
                    #        print "Warning: None alias", id_word, alias, source_list_str
                    #    continue
		    if source_db.strip() == "":
			continue
		    value = self.biana_access._transform_attribute_value_data_type_to_biana_database_attribute_data_type( attribute_identifier=source_db, value=alias )
                    string_protein_object.add_attribute(ExternalEntityAttribute(attribute_identifier = source_db, value = alias, type = type))
           
        return self.biana_access.insert_new_external_entity( externalEntity = string_protein_object ) 

    def _insert_protein_information_from_alias_file_to_database(self, aliases_file):
        if self.verbose:
            print "Creating temporary database for aliases.."
        self.biana_access.db.insert_db_content( self.alias_temp_table.create_mysql_query(), answer_mode = None )
        aliases_file_fd = self._get_file_descriptor(aliases_file)
        #values = []
        #i = 0
        nMaxId = 0
        nMaxAlias = 0
        nMaxSourceList = 0
        for line in aliases_file_fd:
            #i+=1
            if line.startswith('#'):
                continue
            words = line.split('\t')
            tax = words[0]
            if self.verbose and tax is None:
                print "Warning: None taxId:", line
            id = words[1]
            id_word = "%s.%s" % (tax, id)
            alias = words[2].replace("\"", "") #.strip("\"")
            source_list = words[3].split()
            if self.verbose and alias is None:
                print "Warning: None alias:", line
            source_list_str = " ".join(source_list)
            #values.append(("\"%s\"" % id_word, "\"%s\"" % alias, "\"%s\"" % source_list_str))
            if len(id_word) > nMaxId: nMaxId = len(id_word)
            if len(alias) > nMaxAlias: nMaxAlias = len(alias)
            if len(source_list_str) > nMaxSourceList: nMaxSourceList = len(source_list_str)
            #if i > self.N_MAX_ENTRY_AT_ONCE:
            #    self.biana_access.db.insert_db_content( self.biana_access.db._get_multiple_insert_query(self.alias_temp_table.get_table_name(),('id', 'alias', 'source_list'),values), answer_mode = None )
            #    values = []
            #    i = 0
            self.biana_access.db.insert_db_content( self.biana_access.db._get_insert_sql_query(table = self.alias_temp_table.get_table_name(), column_values = [('id',id_word), 
                                                                                                                                                           ('alias', alias),
                                                                                                                                                           ('source_list', source_list_str)] ),                                                                                                                                            answer_mode = None )


        self.biana_access.db._empty_buffer()

        self._close_file_descriptor(aliases_file_fd)
        if self.verbose:
            print "Max id-alias-sourcelist: ", nMaxId, nMaxAlias, nMaxSourceList
        if self.verbose:
            print "Temporary database for aliases is created!"
        return

    def _fetch_given_list_protein_information_from_database(self, list_id):
        return self.biana_access.db.select_db_content( self.biana_access.db._get_select_sql_query(tables=[ self.alias_temp_table.get_table_name() ] , columns=['id', 'alias', 'source_list'], fixed_conditions=[('id', 'IN', "(\"%s\")" % "\",\"".join(list_id), None)] ), answer_mode = "raw") #, remove_duplicates="yes" )

    def _fetch_given_protein_information_from_database(self, id_word):
        return self.biana_access.db.select_db_content( self.biana_access.db._get_select_sql_query(tables=[ self.alias_temp_table.get_table_name() ] , columns=['alias', 'source_list'], fixed_conditions=[('id', '=', id_word)] ),                                                           answer_mode = "raw") #, remove_duplicates="yes" )

    def _remove_protein_information_from_database(self):
        self.biana_access.db.insert_db_content( "DELETE FROM %s" %self.alias_temp_table.get_table_name() )
        self.biana_access.db.insert_db_content( self.alias_temp_table.get_drop_query() )
        return

    def _get_data_file_names(self):
        print self.input_file
        (sequences_file_fd, aliases_file, links_file_fd) = (None, None, None)
        if( not self.input_file.endswith(os.sep) ):
            self.input_file += os.sep
        directoryData = os.path.dirname(self.input_file)+os.sep
        # find string data files -- not checking if different versions exists
        for file in os.listdir(directoryData):
            file = directoryData + file
            #print file
            if fnmatch.fnmatch(file, '*protein.links.detailed*%s*' % self.sourcedb_version):
                links_file_fd = self._get_file_descriptor(file)
            if fnmatch.fnmatch(file, '*protein.sequences*%s*' % self.sourcedb_version):
                sequences_file_fd = self._get_file_descriptor(file)
            if fnmatch.fnmatch(file, '*protein.aliases*%s*' % self.sourcedb_version):
                #aliases_file = file # self._get_file_descriptor(file)
                if file.endswith(".gz"):
                    os.system("gunzip %s" % file)
                    aliases_file = file[:-3] 
                else:
                    aliases_file = file # = self._get_file_descriptor(file)
        print links_file_fd, sequences_file_fd, aliases_file
        return (sequences_file_fd, aliases_file, links_file_fd)

    def _get_file_descriptor(self, file):
        print file
        if file.endswith(".gz"):
            return gzip.open(file,'r')
        else:
            return open(file, 'r')

    def _close_file_descriptor(self, fd):
        #if isinstance(fd, file):
        fd.close()
        #elif isinstance(fd, gzip):
        #gzip.close(fd)

