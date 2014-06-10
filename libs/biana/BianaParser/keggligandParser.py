"""
File        : keggligandParser.py
Author      : Javier Garcia Garcia
Creation    : January 2008
Contents    : fills up tables in database biana with information from kegg ligand database
Called from : 

=======================================================================================================

This file implements a program that fills up tables in database biana with information of kegg ligand databases

"""

from bianaParser import *
from biana.BianaObjects.Sequence import ProteinSequence

class KeggLigandParser(BianaParser):
    """
    Uniprot Parser Class
    """

    name = "kegg_ligand"
    description = "This file implements a program that fills up tables in database biana with information of kegg Ligand database"
    external_entity_definition = ""
    external_entity_relations = ""


    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "KEGG Ligand database",
                             default_script_name = "keggligandParser.py",
                             default_script_description = KeggLigandParser.description,
                             additional_compulsory_arguments = [], #[("kegg_ligand_path=",None,"Path where compound, drug, glycan, and enzyme files are")],
                             additional_optional_arguments = [])
        self.default_eE_attribute = "keggCode"


    def parse_database(self):
        """
        """

        kegg_ligand_path = self.input_file

        if kegg_ligand_path[-1] != os.sep:
            kegg_ligand_path += os.sep


        # General regex
        #continue_field_regex = re.compile("^\s{3,}([^;]+);*$")
        #continue_field_regex = re.compile("^\s{3,}(.+);$")
        continue_field_regex = re.compile("^\s{3,}(.+);*$")
        #field_regex = re.compile("^(\w+)\s+([^;]+);*$")
        field_regex = re.compile("^(\w+)\s+(.+);*$")
        pathway_regex = re.compile("PATH\:\s+(map|rn)(\d+)\s+(.+)$")

        space_regex = re.compile("\s+")
        parenthesis_regex = re.compile("\(.+\)")  # used to eliminate extra information in sequence
        

        pathway_dict_desc = {}         # This will store the pathway kegg code with its description as value
        pathway_dict_components = {}   # This will store in memory the pathway kegoo code as key with a list of its participants as values

        
        kegg_elements_dict = {}   # stores the uniqueID code from kegg and its correspondence with the external entity identifier
                                  # This is used later when inserting relations
        temp_code = None



        # PARSE COMPOUND FILE
        compound_f = file(kegg_ligand_path+"compound","r")

        entry_regex = re.compile("^ENTRY\s+(\w+)\s+.*\s+Compound")
        remark_regex = re.compile("^REMARK\s+Same\sas\:\s+(.+)$")
        formula_regex = re.compile("^FORMULA\s+(.+)$")
        comment_regex = re.compile("^COMMENT\s+(.+)$")

        peptide_regex = re.compile("^ENTRY.+Peptide.+Compound")
        sequence_regex = re.compile("^SEQUENCE\s+(.+)$")
            
        kegg_object = None

        temp_value = []           # List used to store the information of those fields that can have more than a single line
        current_field = None

        temp_pathway_codes = []

        is_peptide = None

        for line in compound_f:

            m = entry_regex.match(line)

            if m:
                if peptide_regex.match(line):
                    is_peptide = 1
                else:
                    is_peptide = None
                    
                if kegg_object is not None:
                    self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
                    kegg_elements_dict[temp_code] = kegg_object.get_id()
                    [ pathway_dict_components[actual_pathway_code].append(kegg_object.get_id()) for actual_pathway_code in temp_pathway_codes ] 

                id_type = None
                if is_peptide:
                    kegg_object = ExternalEntity( source_database = self.database, type="protein" )
                    #id_type = "keggCompound"
                    id_type = "keggCode"
                else:
                    kegg_object = ExternalEntity( source_database = self.database, type="compound" )
                    #id_type = "keggProtein"
                    id_type = "keggCode"
                    
                kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = id_type, value = m.group(1), type="unique" ) )
                temp_code = m.group(1)
                temp_pathway_codes = []

                continue

            m = pathway_regex.search(line)
            if m:
                temp_pathway_codes.append(m.group(2))
                if not pathway_dict_desc.has_key(m.group(2)):
                    pathway_dict_desc[m.group(2)] = m.group(3)
                    pathway_dict_components[m.group(2)] = []
                continue

            new_field = field_regex.match(line)
            
            if new_field:

                if current_field == "NAME":
                    [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value = x, type= "unique") ) for x in temp_value ]
                elif current_field == "FORMULA":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "formula", value = " ".join(temp_value), type="unique") )

                elif current_field == "COMMENT":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value = " ".join(temp_value) ) )

                elif current_field == "SEQUENCE":
                    #eliminate all between parenthesis and split by spaces
                    if is_peptide:
                        sequence_list = space_regex.split( parenthesis_regex.sub('', " ".join(temp_value)).strip() )
                        if "(Disulfide" in sequence_list:
                            print 
                            print parenthesis_regex.sub('', " ".join(temp_value))
                        sequence = [ ProteinSequence.get_aminoacid_code_3to1( code = actual_residue.replace("-NH2","").replace("Acetyl-","").replace("6-Bromo-","").replace("N-Formyl-Met","") ) for actual_residue in sequence_list ]
                        kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "ProteinSequence", value = ProteinSequence("".join(sequence)), type = "unique" ))
                    
                    
                    #kegg_object

                elif current_field == "REMARK":
                    for current_remark_line in temp_value:
                        m = remark_regex.match(current_field+" "+current_field)
                        if m:
                            #print m.group(1)
                            [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value=x,type="cross-reference")) for x in m.group(1).split(" ") ]

                current_field = new_field.group(1)
                temp_value = [new_field.group(2).strip()]
            else:
                cont_value = continue_field_regex.match(line)
                if cont_value:
                    temp_value.append(cont_value.group(1).strip())
                
                

        # Insert the last one
        if kegg_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
            kegg_elements_dict[temp_code] = kegg_object.get_id()
            
        compound_f.close()


        #########################################################################################

        # PARSE DRUG FILE
        drug_f = file(kegg_ligand_path+"drug","r")

        entry_regex = re.compile("^ENTRY\s+(\w+)\s+.*\s+Drug")
        remark_regex = re.compile("^REMARK\s+Same\sas\:\s+(.+)$")
        formula_regex = re.compile("^FORMULA\s+(.+)$")

        kegg_object = None

        temp_value = []           # List used to store the information of those fields that can have more than a single line
        current_field = None

        temp_pathway_codes = []


        for line in drug_f:


            m = entry_regex.match(line)

            if m:
                if kegg_object is not None:
                    self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
                    kegg_elements_dict[temp_code] = kegg_object.get_id()
                    [ pathway_dict_components[actual_pathway_code].append(kegg_object.get_id()) for actual_pathway_code in temp_pathway_codes ]
                
                kegg_object = ExternalEntity( source_database = self.database, type="drug" )
                kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value= m.group(1), type="unique") )
                temp_code = m.group(1)
                temp_pathway_codes = []

                continue


            m = pathway_regex.search(line)
            if m:
                temp_pathway_codes.append(m.group(2))
                if not pathway_dict_desc.has_key(m.group(2)):
                    pathway_dict_desc[m.group(2)] = m.group(3)
                    pathway_dict_components[m.group(2)] = []
                continue
            

            new_field = field_regex.match(line)
            if new_field:
                if current_field == "NAME":
                    [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value=x, type="unique") ) for x in temp_value ]
                elif current_field == "FORMULA":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "formula", value=" ".join(temp_value), type="unique") )

                elif current_field == "COMMENT":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value=" ".join(temp_value)) )

                elif current_field == "REMARK":
                    for current_remark_line in temp_value:
                        m = remark_regex.match(current_field+" "+current_field)
                        if m:
                            #print m.group(1)
                            [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value=x,type="cross-reference")) for x in m.group(1).split(" ") ]

                current_field = new_field.group(1)
                temp_value = [new_field.group(2).strip()]
            else:
                cont_value = continue_field_regex.match(line)
                if cont_value:
                    temp_value.append(cont_value.group(1).strip())        

                
        # Insert the last one
        if kegg_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
            kegg_elements_dict[temp_code] = kegg_object.get_id()
            
        drug_f.close()

        
        #########################################################################################

        # PARSE GLYCAN FILE
        glycan_f = file(kegg_ligand_path+"glycan","r")

        entry_regex = re.compile("^ENTRY\s+(\w+)\s+.*\s+Glycan")
        formula_regex = re.compile("^COMPOSITION\s+(.+)$")
        remark_regex = re.compile("^REMARK\s+Same\sas\:\s+(.+)$")


        kegg_object = None

        temp_value = []           # List used to store the information of those fields that can have more than a single line
        current_field = None

        temp_pathway_codes = []
        
        for line in glycan_f:


            m = entry_regex.match(line)

            if m:
                if kegg_object is not None:
                    self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
                    kegg_elements_dict[temp_code] = kegg_object.get_id()
                    [ pathway_dict_components[actual_pathway_code].append(kegg_object.get_id()) for actual_pathway_code in temp_pathway_codes ]
                
                kegg_object = ExternalEntity( source_database = self.database, type="glycan" )
                kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value=m.group(1),type="unique") )
                temp_code = m.group(1)
                temp_pathway_codes = []

                continue

            m = pathway_regex.search(line)
            if m:
                temp_pathway_codes.append(m.group(2))
                if not pathway_dict_desc.has_key(m.group(2)):
                    pathway_dict_desc[m.group(2)] = m.group(3)
                    pathway_dict_components[m.group(2)] = []
                continue

            new_field = field_regex.match(line)
            if new_field:
                    
                if current_field == "NAME":
                    [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value = x,type="unique") ) for x in temp_value ]
                elif current_field == "COMPOSITION":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "formula", value = " ".join(temp_value), type="unique") )

                elif current_field == "COMMENT":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value = " ".join(temp_value)) )

                elif current_field == "REMARK":
                    for current_remark_line in temp_value:
                        m = remark_regex.match(current_field+" "+current_field)
                        if m:
                            #print m.group(1)
                            [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value=x,type="cross-reference")) for x in m.group(1).split(" ") ]

                current_field = new_field.group(1)
                temp_value = [new_field.group(2).strip()]
            else:
                cont_value = continue_field_regex.match(line)
                if cont_value:
                    temp_value.append(cont_value.group(1).strip())     


        # Insert the last one
        if kegg_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
            kegg_elements_dict[temp_code] = kegg_object.get_id()
            
        glycan_f.close()
        



        #########################################################################################

        # PARSE ENZIME FILE
        enzyme_f = file(kegg_ligand_path+"enzyme","r")

        entry_regex = re.compile("^ENTRY\s+EC\s*([\d\.]+)\s+.*\s+Enzyme")

        kegg_object = None

        temp_value = []           # List used to store the information of those fields that can have more than a single line
        current_field = None

        temp_pathway_codes = []

        sysname_regex = re.compile("^SYSNAME\s+(.+)$")
        structure_regex = re.compile("PDB\:\s+(.+)$")

        for line in enzyme_f:

            m = entry_regex.match(line)

            if m:
                if kegg_object is not None:
                    self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
                    kegg_elements_dict[temp_code] = kegg_object.get_id()
                    [ pathway_dict_components[actual_pathway_code].append(kegg_object.get_id()) for actual_pathway_code in temp_pathway_codes ]
                
                kegg_object = ExternalEntity( source_database = self.database, type="enzyme" )
                kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "EC", value=m.group(1),type="unique") )
                temp_code = m.group(1)
                temp_pathway_codes = []

                continue


            m = pathway_regex.search(line)
            if m:
                temp_pathway_codes.append(m.group(2))
                if not pathway_dict_desc.has_key(m.group(2)):
                    pathway_dict_desc[m.group(2)] = m.group(3)
                    pathway_dict_components[m.group(2)] = []
                continue

            m = sysname_regex.match(line)
            if m:
                kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value=m.group(1), type="unique") )


            new_field = field_regex.match(line)
            if new_field:

                if current_field == "NAME":
                    [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value=x, type="unique") ) for x in temp_value ]

                elif current_field == "COMMENT":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value = " ".join(temp_value)) )

                elif current_field == "STRUCTURES":
                    all_str = " ".join(temp_value).strip()
                    m = structure_regex.search(all_str)
                    if m:
                        for actual_pdb in space_regex.split(m.group(1)):
                            kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "pdb", value = actual_pdb, type="cross-reference") )
                
                current_field = new_field.group(1)
                temp_value = [new_field.group(2).strip()]
            else:
                cont_value = continue_field_regex.match(line)
                if cont_value:
                    temp_value.append(cont_value.group(1).strip())     


        # Insert the last one
        if kegg_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
            kegg_elements_dict[temp_code] = kegg_object.get_id()

        enzyme_f.close()



        #########################################################################################

        # PARSE REACTION FILE
        reaction_f = file(kegg_ligand_path+"reaction","r")

        entry_regex = re.compile("^ENTRY\s+(\w+)\s+.*\s+Reaction")
        enzyme_regex = re.compile("^ENZYME\s+([\d\.\s]+)$")
        equation_regex = re.compile("^EQUATION\s+(.+)\s*\<\=\>\s+(.+)\s*$")

        # special case for dna
        parenthesis_regex = re.compile("\(.+\)")

        kegg_object = None

        temp_value = []           # List used to store the information of those fields that can have more than a single line
        current_field = None

        temp_pathway_codes = []

        for line in reaction_f:

            m = entry_regex.match(line)

            if m:
                if kegg_object is not None:
                    self.biana_access.insert_new_external_entity( externalEntity = kegg_object )
                    kegg_elements_dict[temp_code] = kegg_object.get_id()
                    [ pathway_dict_components[actual_pathway_code].append(kegg_object.get_id()) for actual_pathway_code in temp_pathway_codes ]
                
                kegg_object = ExternalEntityRelation( source_database = self.database, relation_type="reaction" )
                kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value=m.group(1),type="unique"))
                temp_code = m.group(1)
                temp_pathway_codes = []

                continue


            m = pathway_regex.search(line)
            if m:
                temp_pathway_codes.append(m.group(2))
                if not pathway_dict_desc.has_key(m.group(2)):
                    pathway_dict_desc[m.group(2)] = m.group(3)
                    pathway_dict_components[m.group(2)] = []
                continue

            new_field = field_regex.match(line)
            if new_field:
                if current_field == "NAME":
                    [ kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "name", value=x,type="unique") ) for x in temp_value ]

                elif current_field == "COMMENT":
                    kegg_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value=" ".join(temp_value)) )

                elif current_field == "ENZYME":
                    m = enzyme_regex.match(current_field+" "+" ".join(temp_value).strip())
                    if m:
                        #get the externalEntityID for this enzymeID (stored in memory) and add as a participant
                        for actual_enzyme in space_regex.split(m.group(1)):
                            kegg_object.add_participant( externalEntityID = kegg_elements_dict[actual_enzyme] )
                            kegg_object.add_participant_attribute( externalEntityID = kegg_elements_dict[actual_enzyme],
                                                                   participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", 
                                                                                                                                      value = "catalyst" ) )
                            
                elif current_field == "EQUATION":
                    m = equation_regex.match(current_field+" "+" ".join(temp_value))
                    if m:
                        substrates = m.group(1)
                        products = m.group(2)

                        # Add substrates
                        for actual_substrat in [ x.strip() for x in substrates.split(" + ") ]:
                            splitted = actual_substrat.split(" ")
                            if len(splitted)==1:
                                num = 1
                                #code = splitted[0]
                                code = parenthesis_regex.sub('', splitted[0])
                                
                            elif len(splitted)==2:
                                num = splitted[0].replace('n','').replace('(','').replace(')','')
                                #code = splitted[1]
                                code = parenthesis_regex.sub('', splitted[1])
                            else:
                                raise ValueError("How is possible to have more than 2 elements?")

                            try:
                                kegg_object.add_participant( externalEntityID = kegg_elements_dict[code] )
                                kegg_object.add_participant_attribute( externalEntityID = kegg_elements_dict[code],
                                                                       participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", 
                                                                                                                                       value = "substrate" ) )
                                
                                if num != '':
                                    if int(num)>1:
                                        kegg_object.add_participant_attribute( externalEntityID = kegg_elements_dict[code],
                                                                               participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "cardinality", 
                                                                                                                                                   value = num ) )
                            except:
                                sys.stderr.write("Kegg element %s is not defined in kegg database\n" %code)

                        # Add products
                        for actual_product in [ x.strip() for x in products.split(" + ") ]:
                            splitted = actual_product.split(" ")
                            if len(splitted)==1:
                                num = 1
                                #code = splitted[0]
                                code = parenthesis_regex.sub('', splitted[0])
                            elif len(splitted)==2:
                                num = splitted[0].replace('n','').replace('(','').replace(')','')
                                #code = splitted[1]
                                code = parenthesis_regex.sub('', splitted[1])
                            else:
                                raise ValueError("How is possible to have more than 2 elements? [ %s ]\nPRODUCTS: %s" %(actual_product,products))

                            try:
                                kegg_object.add_participant( externalEntityID = kegg_elements_dict[code] )
                                kegg_object.add_participant_attribute( externalEntityID = kegg_elements_dict[code],
                                                                       participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "role", 
                                                                                                                                           value = "product" ) )
                                if num != '':
                                    if int(num)>1:
                                        kegg_object.add_participant_attribute( externalEntityID = kegg_elements_dict[code],
                                                                               participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "cardinality", 
                                                                                                                                                   value = num ) )
                            except:
                                sys.stderr.write("Kegg element %s is not defined in kegg database\n" %code)

                current_field = new_field.group(1)
                temp_value = [new_field.group(2).strip()]
            else:
                cont_value = continue_field_regex.match(line)
                if cont_value:
                    temp_value.append(cont_value.group(1).strip())     


        # Insert the last one
        if kegg_object is not None:
            self.biana_access.insert_new_external_entity( externalEntity = kegg_object )

        reaction_f.close()




        ###############################

        # INSERT PATHWAYS INFORMATION

        # Each distinct pathway is inserted as an external entity
        # All the elements mapped to the pathway are assigned to them
        # Maybe it would be sufficient to map the reactions, as the components are mapped to the reactions...

        for actual_pathway_code in pathway_dict_desc.keys():
            kegg_pathway_object = ExternalEntityRelation( source_database = self.database, relation_type="pathway" )
            kegg_pathway_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "keggCode", value=actual_pathway_code, type="unique") )
            kegg_pathway_object.add_attribute( ExternalEntityAttribute( attribute_identifier = "description", value=pathway_dict_desc[actual_pathway_code] ) )

            [ kegg_pathway_object.add_participant( externalEntityID = actual_participant_external_entity_id )
              for actual_participant_external_entity_id in pathway_dict_components[actual_pathway_code] ]

            self.biana_access.insert_new_external_entity( externalEntity = kegg_pathway_object )


            

