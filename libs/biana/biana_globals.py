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


# This file contains specific default user parameters used in biana

# DATABASE CONNECTION DEFAULT PARAMETERS
DBNAME = 'biana'
DBUSER = 'root'
DBPASS = 'root'
DBHOST = 'localhost'
DBPORT = None
DBSOCKET = '/home/jgarcia/local/mysql/var/mysql.sock'

EXTERNAL_ENTITY_TYPES = ["protein",
                         "DNA",
                         "RNA",
                         "mRNA",
                         "tRNA",
                         "rRNA",
                         "CDS",
                         "gene",
                         "sRNA",
                         "snRNA",
                         "snoRNA",
                         "structure",
                         "pattern",
                         "compound",
                         "drug",
                         "glycan",
                         "enzyme",
                         "relation",
                         "ontology",
                         "SCOPElement",
                         "taxonomyElement",
                         "PsiMiOboOntologyElement",
                         "GOElement"
                         ]





EXTERNAL_ENTITY_RELATION_TYPES = [ "interaction",
                                   "no_interaction",
                                   "reaction",
                                   "functional_association",
                                   "cluster",
                                   "homology",
                                   "pathway",
                                   "alignment",
                                   "complex",
                                   "regulation",
                                   "cooperation",
                                   "forward_reaction",
                                   "backward_reaction" 
                                   ]

# EXTERNAL ENTITY ATTRIBUTE TYPES
EXTERNAL_ENTITY_IDENTIFIER_ATTRIBUTES = [ ("CHEBI", "integer unsigned"),
                                          ("COG", "varchar(10)"),
                                          ("CYGD", "varchar(15)"), # normally 7 (YDR172w) but sometimes 9 (YLR312w-a) (in mips there are some errors... because of that, we increase it to 15
                                          ("DIP", "varchar(6)"), # DIP:216N (~17000 entries)
                                          ("EC", "varchar(30)"),
                                          ("Encode", "varchar(14)"),
                                          ("Ensembl", "varchar(40)"),
                                          ("FlyBase", "varchar(13)"),
                                          ("GDB", "integer(3) unsigned"),
                                          ("GeneID", "integer(4) unsigned"),
                                          ("GeneSymbol", "varchar(255)"),
                                          ("GenomeReviews", "varchar(15)"),
                                          ("GI", "integer(4) unsigned"),
                                          ("GO", "integer(3) unsigned"),
                                          ("HGNC", "integer(2) unsigned"),
                                          ("Homologene", "integer(3) unsigned"),
                                          ("HPRD", "integer(3) unsigned"),
                                          ("Huge", "smallint unsigned"),
                                          ("IMGT", "varchar(10)"),
                                          ("IntAct", "integer(3) unsigned"),
                                          ("IntEnz", "varchar(10)"),
                                          ("InterPro", "varchar(12)"),
                                          #("IPI", "varchar(20)"), # Moved to versionable attributes
                                          ("KeggCode", "char(6)"),
                                          ("KeggGene", "varchar(155)"),
                                          ("Method_id", "integer(2) unsigned"),    #psi_mi obo mi code
                                          ("MGI", "integer(3) unsigned"),
                                          ("MIM", "integer(3) unsigned"),
                                          ("MINT", "integer(3) unsigned"),
                                          ("MIPS", "integer(2) unsigned"),
                                          ("OrderedLocusName", "varchar(255)"),
                                          ("ORFName", "varchar(255)"), # Actually at most 7: YAL213W: Yeast (Y) 1st (A) chromosome's left (L) at 213th (213) position on Watson (W) strand
                                          ("PFAM", "varchar(255)"),
                                          ("PIR", "varchar(8)"),
                                          ("PRINTS", "varchar(15)"),
                                          ("PRODOM", "varchar(15)"),
                                          ("Prosite", "varchar(255)"),
                                          ("psimi_name", "varchar(255)"),
                                          ("PubChemCompound", "integer(3) unsigned"),
                                          ("Ratmap", "integer(3) unsigned"),
                                          ("Reactome", "integer unsigned"),
                                          ("RGD", "integer unsigned"),
                                          ("SCOP", "integer(3) unsigned"), 
                                          ("SGD", "varchar(15)"),
                                          ("STRING", "varchar(25)"), # gives ordered locus names, so called ensembl codes and many more
                                          ("Tair", "varchar(100)"),
                                          ("TaxID", "integer(3) unsigned"),
                                          ("Unigene", "varchar(10)"),
                                          ("UniParc", "binary(10)"),
                                          ("UniprotEntry", "varchar(15)"),
                                          ("WormBaseGeneID", "integer(3) unsigned"),
                                          ("WormBaseSequenceName", "varchar(255)"),
                                          ("YPD", "varchar(15)"),
                                          ("iRefIndex_ROGID", "varchar(255)"),
                                          ("iRefIndex_RIGID", "varchar(255)"),
                                          ]


EXTERNAL_ENTITY_GENERAL_ATTRIBUTES = []


PROMISCUOUS_EXTERNAL_ENTITY_TYPES_DICT = [ ("SCOPElement", "PDB") ]


VALID_IDENTIFIER_REFERENCE_TYPES = ["unique", "previous", "alias", "cross-reference", "synonym","short-name", "exact_synonym", "related_synonym"]

CROSSABLE_ATTRIBUTES = set(["sequence","taxid","ipi","uniprotentry","uniprotaccession","genesymbol","geneid","refseq","ec"])


EXTERNAL_ENTITY_VERSIONABLE_IDENTIFIER_ATTRIBUTE_TYPES = [("AccessionNumber", "varchar(15)"),
                                                          ("RefSeq", "varchar(15)"),
                                                          ("TIGR", "varchar(255)"),
                                                          ("UniprotAccession", "varchar(9)"),
							  ("IPI", "varchar(20)"),
                                                          ]


EXTERNAL_ENTITY_DESCRIPTIVE_SEARCHABLE_ATTRIBUTE_TYPES = [("Disease", "text(2000)"),
                                                          ("Function", "text"),
                                                          ("Keyword", "varchar(50)"),
                                                          ("Description", "text"),
                                                          ("SubcellularLocation", "text(400)"),
                                                          ("Name", "varchar(255)")
                                                          ]

EXTERNAL_ENTITY_DESCRIPTIVE_ATTRIBUTE_TYPES = [("Pubmed", "integer(3) unsigned"),
                                               ("Formula", "varchar(255)")
                                               ]


EXTERNAL_ENTITY_NUMERIC_ATTRIBUTE_TYPES = [("Pvalue", "double"),
                                           ("Score", "double"),
                                           ("iRefIndex_lpr", "integer unsigned"),
                                           ("iRefIndex_hpr", "integer unsigned"),
                                           ("iRefIndex_np", "integer unsigned"),
                                           ("STRINGScore", "int(2)"),
                                           ("STRINGScore_neighborhood","int(2)"), 
                                           ("STRINGScore_fusion","int(2)"), 
                                           ("STRINGScore_cooccurence","int(2)"), 
                                           ("STRINGScore_coexpression","int(2)"),
                                           ("STRINGScore_experimental","int(2)"), 
                                           ("STRINGScore_db", "int(2)"), 
                                           ("STRINGScore_textmining","int(2)")]


EXTERNAL_ENTITY_SPECIAL_ATTRIBUTE_TYPES = { "PDB": {"fields": [ ("value","char(4)"),
                                                                ("chain","varchar(4)",True),
                                                                ("pdb_range","varchar(255)",True) ], 
                                                    "indices": ("value","chain","pdb_range")},
                                            
                                            "ProteinSequence": { "fields": [ ("value","binary(16)"), 
                                                                             ("sequenceType","ENUM(\"peptide\")",False) ],
                                                                 "indices": ("value",) },
                                            
                                            "NucleotideSequence": { "fields": [ ("value","binary(16)"),
                                                                                ("sequenceType","ENUM(\"dna\",\"rna\")",False)],
                                                                    "indices": ("value",)},
                                            
                                            "SequenceMap": { "fields": [ ("value","binary(16)"),
                                                                         ("seq_range","varchar(255)",False) ], 
                                                             "indices": ()},
                                            
                                            "Pattern": { "fields": [ ("value","varchar(255)"),
                                                                     ("patternExpression","varchar(255)",False)], 
                                                         "indices": ("value",)}, # Stores a regex
                                            
                                            #"STRINGScore": { "fields": [ ("value","int(2)"), # moved to regular attributes as seperate score attributes
                                            #                             ("neighborhood","int(2)",True), 
                                            #                             ("fusion","int(2)",True), 
                                            #                             ("cooccurence","int(2)",True), 
                                            #                             ("coexpression","int(2)",True),
                                            #                             ("experimental","int(2)",True), 
                                            #                             ("db", "int(2)", True), 
                                            #                             ("textmining","int(2)",True)],
                                            #                 "indices": () }
                                            }


# EXTERNAL ENTITY RELATION PARTICIPANT ATTRIBUTE TYPES 
EXTERNAL_ENTITY_RELATION_PARTICIPANT_ATTRIBUTE_TYPES = [ ("cardinality", "smallint unsigned"),
                                                         ("detection_method", "smallint unsigned"),
                                                         ("GO", "integer(3) unsigned"),
                                                         ("KeggCode", "varchar(6)"),
                                                         ("role", "ENUM(\"batch\",\"product\",\"substrate\",\"catalyst\",\"prey\",\"bait\",\"neutral\",\"acceptor\",\"donor\",\"self\",\"ancillary\",\"enzyme\",\"enzyme target\",\"inhibitor\",\"cofactor\",\"stimulator\",\"activates\",\"inhibits\",\"allosteric_inhibition\",\"competitive_inhibition\",\"irreversible_inhibition\",\"non_competitive_inhibition\",\"uncompetitive_inhibition\",\"allosteric_activation\",\"nonallosteric_activation\",\"transcription_factor\",\"regulated_DNA\",\"onward_effect\",\"reverse_effect\")")
                                                         ]





# EXTERNAL SOFTWARE EXECUTABLES
CLUSTALW_EXEC = '/soft/bio/sequence/clustalw' # '/usr/local/bin/clustalw'
FORMATDB_EXEC = '/soft/bio/sequence/blast-2.2.17/bin/formatdb'
BLASTALL_EXEC = '/soft/bio/sequence/blast-2.2.17/bin/blastall'
BL2SEQ_EXEC = '/soft/bio/sequence/blast-2.2.17/bin/bl2seq'
TCOFFEE_EXEC = '/soft/bio/sequence/t-coffee-6.92' #'/Users/javigx2/phD/external_software/tcoffee/T-COFFEE_distribution_Version_7.04/bin/macosx/t_coffee'
#DSSP_EXEC = '/Users/javigx2/phD/external_software/dssp/dssp/dsspcmbi'
DSSP_EXEC = '/soft/bio/structure/dssp/dsspcmbi'
#CD_HIT_PATH = '/home/emre/lib/cd-hit/'
CD_HIT_PATH = '/soft/bio/sequence/cd-hit/'
#CD_HIT_PATH = '/home/jgarcia/programs/CD-HIT/cd-hit/'
# DEFAULT TEMPORAL DATA PATHS
TEMPORAL_PATH = None

