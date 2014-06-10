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
from obo_index import obo_name_to_MI

class iRefIndexParser(BianaParser):
    """
    Parser for iRefIndex MITAB 2.5 PPI files
    """

    name = "iRefIndex"
    description = "This file implements the parser for iRefIndex"
    external_entity_definition = "A protein"
    external_entity_relations = "physical interaction"


    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "iRefIndex Database",
                             default_script_name = "iRefIndexParser.py",
                             default_script_description = iRefIndexParser.description,
                             additional_optional_arguments = [])
        self.default_eE_attribute = "iRefIndex_ROGID"
        self.initialize_input_file_descriptor()
        return


    def parse_database(self):
        """
        Parses the iRefIndex file

        iRefIndex format explained at http://irefindex.uio.no/wiki/README_iRefIndex_MITAB_4.0
        """

        #irefindex:NeQ0QcVrpNzyqWlmp/HvG5FIRHA9606       irefindex:xjgJ54UxwS5DwnuMgn6A8XxjvsY9606       uniprotkb:P24593|refseq:NP_000590|entrezgene/locuslink:3488     uniprotkb:P35858|refseq:NP_004961|entrezgene/locuslink:3483  uniprotkb:IBP5_HUMAN|entrezgene/locuslink:IGFBP5        uniprotkb:ALS_HUMAN|entrezgene/locuslink:IGFALS MI:0000(-)      -       pubmed:-        taxid:9606      taxid:9606  MI:0000(aggregation)     MI:0923(irefindex)|MI:0000:(ophid)      irefindex:++1ALo0GitoAXTTz5lyLFlYtoO8|ophid:-   lpr:-|hpr:-|np:-        3488    3483    MI:0326(protein)        MI:0326(protein)    ++1ALo0GitoAXTTz5lyLFlYtoO8      X       2       3196230 4869527
        
        fields_dict = {}

        external_entities_dict = {}   # Dictionary protein_id -> external Entity ID

        external_entity_relations_dict = {}

        mi_re = re.compile("MI:\d+\((.+)\)")

        for line in self.input_file_fd:

            # PROCESS HEADER
            #uidA   uidB    altA    altB    aliasA  aliasB  method  author  pmids   taxa    taxb    interactionType sourcedb        interactionIdentifiers  confidence      entrezGeneA     entrezGeneB     atype
            #btype   rigid   edgetype        numParticipants ROGA    ROGB
            if line[0]=="#":
                fields = line[1:].strip().split("\t")
                for x in xrange(0, len(fields)):
                    fields_dict[fields[x].lower()] = x
            
                continue

            fields = line.strip().split("\t")

            eE1_id = None
            eE2_id = None

            # Create or get external entities
            if "0326" in fields[fields_dict["atype"]]:

                if fields[0] not in external_entities_dict:

                    eE1 = ExternalEntity( source_database = self.database, type="protein" )

                    # Primary ids
                    primary_id_t = fields[fields_dict["uida"]]
                    t = primary_id_t.split(":")
                    eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "irefindex_ROGID", value=t[1].replace("irefindex:",""), type="unique") )


                    alternative_ids = fields[fields_dict["alta"]].split("|")
                    for current_id in alternative_ids:
                        if current_id == "-":
                            continue
                        t = current_id.split(":")
                        idtype = t[0].lower()
                        idvalue = t[1]
                        if idvalue=="-":
                            continue
                        if idtype=="uniprotkb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value=idvalue, type="cross-reference") )
                        elif idtype=="refseq":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "refseq", value=idvalue, type="cross-reference") )
                        elif idtype=="entrezgene/locuslink":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "geneID", value=idvalue, type="cross-reference") )
                        elif idtype=="pdb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "pdb", value=idvalue[0:4], type="cross-reference", additional_fields = {"chain": idvalue[-1]} ) )
                        elif idtype=="gb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="dbj":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="pir":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "pir", value = idvalue, type="cross-reference" ) )
                        elif idtype=="kegg":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "keggCode", value = idvalue, type="cross-reference" ) )
                        elif idtype=="emb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="uniprot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="swiss-prot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="intact":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "intact", value = idvalue.replace("EBI-",""), type="cross-reference" ) )
                        elif idtype=="genbank_protein_gi":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "GI", value = idvalue, type="cross-reference" ) )
                        elif idtype=="tigr":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "tigr", value = idvalue, type="cross-reference" ) )
                        elif idtype=="prf" or idtype=="pubmed" or idtype=="uniparc":
                            pass
                        else:
                            pass
                            #sys.stderr.write("Alternative id type %s not recognized\n" %idtype)

                    aliases_ids = fields[fields_dict["aliasa"]].split("|")
                    for current_id in aliases_ids:
                        if current_id == "-":
                            continue
                        t = current_id.split(":")
                        idtype = t[0].lower()
                        idvalue = t[1]
                        if idvalue=="-":
                            continue
                        if idtype=="uniprotkb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotentry", value=idvalue, type="cross-reference") )
                        elif idtype=="refseq":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "refseq", value=idvalue, type="cross-reference") )
                        elif idtype=="entrezgene/locuslink":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "geneSymbol", value=idvalue, type="cross-reference") )
                        elif idtype=="pdb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "pdb", value=idvalue[0:4], type="cross-reference", additional_fields = {"chain": idvalue[-1]} ) )
                        elif idtype=="gb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="dbj":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="pir":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "pir", value = idvalue, type="cross-reference" ) )
                        elif idtype=="kegg":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "keggCode", value = idvalue, type="cross-reference" ) )
                        elif idtype=="emb":
                            eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="uniprot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="swiss-prot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="intact":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "intact", value = idvalue.replace("EBI-",""), type="cross-reference" ) )
                        elif idtype=="genbank_protein_gi":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "GI", value = idvalue, type="cross-reference" ) )
                        elif idtype=="tigr":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "tigr", value = idvalue, type="cross-reference" ) )
                        elif idtype=="prf" or idtype=="pubmed" or idtype=="uniparc":
                            pass
                        else:
                            #sys.stderr.write("Alternative id type %s not recognized\n" %idtype)
                            pass


                    taxID = fields[fields_dict["taxa"]].replace("taxid:","")
                    if taxID!="-":
                        eE1.add_attribute( ExternalEntityAttribute( attribute_identifier= "taxID", value=taxID, type = "cross-reference") )

                    external_entities_dict[fields[0]] = self.biana_access.insert_new_external_entity( externalEntity = eE1 )

                eE1_id = external_entities_dict[fields[0]]
                

            if "0326" in fields[fields_dict["btype"]]:

                if fields[1] not in external_entities_dict:

                    eE2 = ExternalEntity( source_database = self.database, type="protein" )

                    # Primary ids
                    primary_id_t = fields[fields_dict["uidb"]]
                    t = primary_id_t.split(":")
                    eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "irefindex_ROGID", value=t[1].replace("irefindex:",""), type="cross-reference") )

                    alternative_ids = fields[fields_dict["altb"]].split("|")
                    for current_id in alternative_ids:
                        if current_id == "-":
                            continue
                        t = current_id.split(":")
                        idtype = t[0].lower()
                        idvalue = t[1]
                        if idvalue=="-":
                            continue
                        if idtype=="uniprotkb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value=idvalue, type="cross-reference") )
                        elif idtype=="refseq":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "refseq", value=idvalue, type="cross-reference") )
                        elif idtype=="entrezgene/locuslink":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "geneID", value=idvalue, type="cross-reference") )
                        elif idtype=="pdb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "pdb", value = idvalue[0:4], type="cross-reference", additional_fields = {"chain": idvalue[-1]} ) )
                        elif idtype=="gb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="dbj":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="pir":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "pir", value = idvalue, type="cross-reference" ) )
                        elif idtype=="kegg":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "keggCode", value = idvalue, type="cross-reference" ) )
                        elif idtype=="emb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="uniprot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="swiss-prot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="intact":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "intact", value = idvalue.replace("EBI-",""), type="cross-reference" ) ) 
                        elif idtype=="genbank_protein_gi":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "GI", value = idvalue, type="cross-reference" ) )
                        elif idtype=="tigr":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "tigr", value = idvalue, type="cross-reference" ) )
                        elif idtype=="prf" or idtype=="pubmed" or idtype=="uniparc":
                            pass
                        else:
                            #sys.stderr.write("Alias id type %s not recognized\n" %idtype)
                            pass

                    aliases_ids = fields[fields_dict["aliasb"]].split("|")
                    for current_id in aliases_ids:
                        if current_id == "-":
                            continue
                        t = current_id.split(":")
                        idtype = t[0].lower()
                        idvalue = t[1]
                        if idvalue=="-":
                            continue
                        if idtype=="uniprotkb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotentry", value=idvalue, type="cross-reference") )
                        elif idtype=="refseq":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "refseq", value=idvalue, type="cross-reference") )
                        elif idtype=="entrezgene/locuslink":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "geneSymbol", value=idvalue, type="cross-reference") )
                        elif idtype=="pdb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "pdb", value=idvalue[0:4], type="cross-reference", additional_fields = {"chain": idvalue[-1]} ) )
                        elif idtype=="gb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="dbj":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="pir":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "pir", value = idvalue, type="cross-reference" ) )
                        elif idtype=="kegg":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "keggCode", value = idvalue, type="cross-reference" ) )
                        elif idtype=="emb":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "accessionnumber", value = idvalue, type="cross-reference" ) )
                        elif idtype=="uniprot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="swiss-prot":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "uniprotaccession", value = idvalue, type="cross-reference" ) )
                        elif idtype=="intact":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "intact", value = idvalue.replace("EBI-",""), type="cross-reference" ) )
                        elif idtype=="genbank_protein_gi":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "GI", value = idvalue, type="cross-reference" ) )
                        elif idtype=="tigr":
                            eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "tigr", value = idvalue, type="cross-reference" ) )
                        elif idtype=="prf" or idtype=="pubmed" or idtype=="uniparc":
                            pass
                        else:
                            #sys.stderr.write("Alias id type %s not recognized\n" %idtype)
                            pass

                    taxID = fields[fields_dict["taxb"]].replace("taxid:","")
                    if taxID!="-":
                        eE2.add_attribute( ExternalEntityAttribute( attribute_identifier= "taxID", value=taxID, type = "cross-reference") )

                    external_entities_dict[fields[1]] = self.biana_access.insert_new_external_entity( externalEntity = eE2 )

                eE2_id = external_entities_dict[fields[1]]


            ######################################
            ## INTERACTION SPECIFIC INFORMATION ##
            ######################################

            rigid_id = fields[fields_dict["rigid"]]
            

            add_common_attributes = True

            # Does the edge represent a binary interaction (X), complex membership (C), or a multimer (Y)?
            edgetype = fields[fields_dict["edgetype"]]
            if edgetype=="X":
                eEr = ExternalEntityRelation( source_database = self.database, relation_type = "interaction" )
            elif edgetype=="C":
                if rigid_id in external_entity_relations_dict:
                    add_common_attributes = False
                eEr = external_entity_relations_dict.setdefault(rigid_id, ExternalEntityRelation( source_database = self.database, relation_type = "complex" ) )
            elif edgetype=="Y":
                eEr = ExternalEntityRelation( source_database = self.database, relation_type = "interaction" )

                
            if "0326" in fields[fields_dict["atype"]]:
                eEr.add_participant( externalEntityID = eE1_id )
            
            if "0326" in fields[fields_dict["btype"]]:
                eEr.add_participant( externalEntityID = eE2_id )


            if edgetype=="Y":
                eEr.add_participant_attribute(externalEntityID = eE1_id,
                                              participantAttribute = ExternalEntityRelationParticipantAttribute( attribute_identifier = "cardinality", 
                                                                                                                 value = fields[fields_dict["numparticipants"]] ))
            else:
                if "0326" in fields[fields_dict["btype"]]:
                    eEr.add_participant( externalEntityID = eE2_id )
                                              

            
            if add_common_attributes:

                eEr.add_attribute( ExternalEntityAttribute( attribute_identifier = "iRefIndex_RIGID", value = rigid_id, type="unique" ) )

                # pubmed:9199353|pubmed:10413469|pubmed:14759368|pubmed:11805826
                pubmeds = fields[fields_dict["pmids"]].split("|")
                for current_pubmed in pubmeds:
                    pubmed_id = current_pubmed.replace("pubmed:","").strip()
                    if pubmed_id != "-" and not pubmed_id.startswith("unassigned") and pubmed_id!="-2" and pubmed_id!="null":
                        eEr.add_attribute( ExternalEntityAttribute( attribute_identifier= "pubmed", value=current_pubmed.replace("pubmed:",""), type="cross-reference") )


                # METHOD
                methods = fields[fields_dict["method"]].split("|")
                # MI:0000(two-hybrid-test)|MI:0000(2h fragment pooling)|MI:0000(two hybrid pooling)
                for current_method in methods:
                    m = mi_re.match(current_method)
                    if m:
                        try:
                            method_name = m.group(1).lower()
                            if method_name=="-" or method_name=="other" or method_name=="not-specified" or method_name=="na":
                                continue
                            method_MI = obo_name_to_MI[method_name]
                            if method_MI is not None:
                                eEr.add_attribute( ExternalEntityAttribute( attribute_identifier= "method_id", value=method_MI, type="cross-reference" ) )
                        except:
                            sys.stderr.write("Method MI not found: %s\n" %m.group(1))

                
                # CONFIDENCE
                confidences = fields[fields_dict["confidence"]].split("|")
                for current_confidence in confidences:
                    if current_confidence.startswith("lpr"):
                        eEr.add_attribute( ExternalEntityAttribute( attribute_identifier= "iRefIndex_lpr", value=current_confidence[3:] ) )
                    elif current_confidence.startswith("hpr"):
                        eEr.add_attribute( ExternalEntityAttribute( attribute_identifier= "iRefIndex_hpr", value=current_confidence[3:] ) )
                    elif current_confidence.startswith("np"):
                        eEr.add_attribute( ExternalEntityAttribute( attribute_identifier= "iRefIndex_np", value=current_confidence[2:] ) )

            if edgetype != "C":
                self.biana_access.insert_new_external_entity( externalEntity = eEr )


        # Insert all complexes
        for current_complex_eEr in external_entity_relations_dict.values():
            self.biana_access.insert_new_external_entity( externalEntity = current_complex_eEr )



            


                    

                    
            
            

