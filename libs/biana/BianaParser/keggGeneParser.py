"""
File        : keggGeneParser.py
Author      : Javier Garcia Garcia
Creation    : January 2008
Contents    : fills up tables in database biana with information from kegg gene database
Called from : 

=======================================================================================================

This file implements a program that fills up tables in database biana with information of kegg gene databases

"""

from bianaParser import *

class KeggGeneParser(BianaParser):
    """

    """

    name = "kegg_gene"
    description = "This file implements a program that fills up tables in database biana with information of kegg Gene Database"
    external_entity_definition = "A external entity represents a gene"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "KEGG GENE database",
                             default_script_name = "keggGeneParser.py",
                             default_script_description = KeggGeneParser.description )
        self.default_eE_attribute = "keggGene"
        self.initialize_input_file_descriptor()

    def parse_database(self):
        """
        """

        # General regex
        continue_field_regex = re.compile("^\s{3,}([^;]+);*$")
        field_regex = re.compile("^(\w+)\s+([^;]+);*$")
        pathway_regex = re.compile("PATH\:\s+(map|rn)(\d+)\s+(.+)$")
        ec_regex = re.compile("\[EC\:([\d\.])+\]")

        space_regex = re.compile("\s+")
        parenthesis_regex = re.compile("\(.+\)")  # used to eliminate extra information in sequence

        # In this case, the entry contains information about the specie
        #ENTRY       ZMO0001           CDS       Z.mobilis
        entry_regex = re.compile("ENTRY\s+(\w+)\s+([\w\_]+)\s+([\w\.]+)$")

        dblink_split_regex = re.compile("(\w+)\:")

        kegg_gene_object = None

        temp_value = []           # List used to store the information of those fields that can have more than a single line
        current_field = None

        number_of_entries = 0

        dict_name_tax = self.biana_access.get_taxonomy_names_taxID_dict()
        new_dict_name_tax = {}

        if len(dict_name_tax)==0:
            print "Taxonomy won't be inserted as Taxonomy database has not been previously inserted"
            
        # Transform species name
        for current_tax_name in dict_name_tax:
            splitted = current_tax_name.split(" ")
            if( len(splitted)==2 ):
                new_dict_name_tax[current_tax_name[0].upper()+"."+splitted[1]] = dict_name_tax[current_tax_name]

        del dict_name_tax
        dict_name_tax = new_dict_name_tax

        not_recognized_tax_id_names = set()

        for line in self.input_file_fd:

            m = entry_regex.search(line)

            if m:
                
                if kegg_gene_object is not None:
                    self.biana_access.insert_new_external_entity( externalEntity = kegg_gene_object )

                # It should be gene or protein... to check!!!
                if m.group(2) == "misc_RNA":
                    type = "RNA"
                elif m.group(2) == "tRNA":
                    type = "tRNA"
                elif m.group(2) == "rRNA":
                    type = "rRNA"
                elif m.group(2) == "mRNA":
                    type = "mRNA"
                elif m.group(2) == "CDS":
                    type = "CDS"
                elif m.group(2) == "snRNA":
                    type = "snRNA"
                elif m.group(2) == "snoRNA":
                    type = "snoRNA"
                elif m.group(2) == "gene":
                    type = "gene"
                else:
                    print "type %s not recognized..." %(m.group(2))


                kegg_gene_object = ExternalEntity( source_database = self.database, type = type )
                kegg_gene_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggGene", value = m.group(1), type =  "unique" ) )

                try:
                    kegg_gene_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "taxID", value = dict_name_tax[m.group(3)], type = "unique") )
                except:
                    not_recognized_tax_id_names.add(m.group(3))


                number_of_entries += 1
                if self.time_control:
                    if number_of_entries%20000==0:
                        sys.stderr.write("%s entries done in %s seconds\n" %(number_of_entries,time.time()-self.initial_time))
                
                
                continue


            new_field = field_regex.match(line)
            if new_field:
                if current_field == "DEFINITION":
                    kegg_gene_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "description", value = " ".join(temp_value) ) )
                    
                    ec_match = ec_regex.search("".join(temp_value))
                    if ec_match:
                        kegg_gene_object.add_attribute( ExternalEntityAttribute(attribute_identifier = "ec", value = ec_match.group(1), type="cross-reference" ) )

                    
                if current_field == "DBLINK":
                    all_db_links = " ".join(temp_value)
                    list_db_links = [ x.strip() for x in dblink_split_regex.split(all_db_links) ]

                    for actual_position in xrange(len(list_db_links)):
                        if list_db_links[actual_position] == "NCBI-GI":
                            [ kegg_gene_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "gi", value = x, type="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]
                            
                        elif list_db_links[actual_position] == "NCBI-GeneID":
                            [ kegg_gene_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "geneID", value=x, type="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]
                            
                        elif list_db_links[actual_position] == "UniProt":
                            [ kegg_gene_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "uniprotAccession", value=x, type ="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]
                        elif list_db_links[actual_position] == "TIGR":
                            [ kegg_gene_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "tigr", value=x, type ="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]

                # MOTIF not used because not standard nomenclature...
                elif current_field == "MOTIF":
                    all_db_links = " ".join(temp_value)
                    list_db_links = [ x.strip() for x in dblink_split_regex.split(all_db_links) ]
                    for actual_position in xrange(len(list_db_links)):
                        if list_db_links[actual_position] == "Pfam":
                            [ kegg_gene_object.add_attribute(ExternalEntityAttribute( attribute_identifier ="pfam", value=x,type="cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]
                        elif list_db_links[actual_position] == "PROSITE":
                            [ kegg_gene_object.add_attribute(ExternalEntityAttribute( attribute_identifier = "prosite", value=x, type = "cross-reference")) for x in list_db_links[actual_position+1].split(" ") ]

                elif current_field == "AASEQ":
                    aa_seq = "".join(temp_value[1:])
		    if len(aa_seq.strip()) > 0:
			kegg_gene_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "proteinSequence", value = ProteinSequence(aa_seq), type = "unique" ) )
                    
                elif current_field == "NTSEQ":
                    nn_seq = "".join(temp_value[1:])
                    kegg_gene_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "nucleotideSequence", value = DNASequence(nn_seq), type = "unique" ))
                    

                current_field = new_field.group(1)
                temp_value = [new_field.group(2)]
            else:
                cont_value = continue_field_regex.match(line)
                if cont_value:
                    temp_value.append(cont_value.group(1))        

                
        # Insert the last one
        if kegg_gene_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = kegg_gene_object )
            
        print "Not recognized specie names: \n%s" %"\n".join(not_recognized_tax_id_names)

        
