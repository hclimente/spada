#!/soft/devel/python-2.7/bin/python

from biana import *
from biana.utilities import identifier_utilities
from rpy2.robjects import r
import sys

top = 10
if(len(sys.argv) == 1):
	top = sys.argv[1]

# Start Biana Session
bianaSession = create_new_session(
									sessionID = 'SmartAS', 
									dbname = 'BIANA_MARCH_2013',
									dbhost = 'ben-yehuda',
									dbuser = 'biana_user',
									dbpassword = 'biana_password', 
									unification_protocol = 'uniprot_geneID_seqtax'
								  )

# Create A List With All The Seed Identifiers

list_input_identifiers = []

with open("Results/candidateList.lst", "r") as candidates:
	for line in candidates:
		elements = line.split("\t")
		list_input_identifiers.append(("ensembl", elements[1]))
		list_input_identifiers.append(("ensembl", elements[2]))

#list_input_identifiers = identifier_utilities.read_identifier_list_from_file(file_name = "ENSTranscripts.lst", id_type = "ensembl")
list_input_restriction_identifiers = []
list_input_negative_restriction_identifiers = []

# Create The Set
user_entity_set_object = bianaSession.create_new_user_entity_set( 
																	identifier_description_list = list_input_identifiers, 
																	attribute_restriction_list = list_input_restriction_identifiers,
																	negative_attribute_restriction_list = list_input_negative_restriction_identifiers,					
																	id_type = 'embedded', 
																	new_user_entity_set_id = "SmartAS_entitySet"
																)

# Create Network Commands
bianaSession.create_network(
								user_entity_set_id = 'SmartAS_entitySet', 
								level = 1, 
								include_relations_last_level = False, 
								relation_type_list = ["interaction"], 
								relation_attribute_restriction_list = [], 
								use_self_relations = True,
								expansion_attribute_list = [],
								expansion_relation_type_list = [], 
								expansion_level = 1, 
								attribute_network_attribute_list = [], 
								group_relation_type_list = []
							)

# Output Commands
bianaSession.output_user_entity_set_details(
												user_entity_set_id = 'SmartAS_entitySet', 
												out_method = open('Results/candidatesInteractions.tsv','w').write, 
												attributes = ["ensembl","uniprotaccession","uniprotentry"], 
												include_level_info = True,
												include_degree_info = True,
												level = None,
												only_selected = True, 
												output_format = 'tabulated', 
												include_tags_info = False,
												include_tags_linkage_degree_info = [], 
												output_1_value_per_attribute = False
											)

bianaSession.output_user_entity_set_details(	
												user_entity_set_id = 'SmartAS_entitySet', 
												out_method = open('Results/candidatesInteractions_extended.tsv','w').write, 
												attributes = ["ensembl","uniprotaccession","uniprotentry"], 
												include_level_info = True,
												include_degree_info = True,
												level = None,
												only_selected = False, 
												output_format = 'tabulated', 
												include_tags_info = False,
												include_tags_linkage_degree_info = [], 
												output_1_value_per_attribute = False
											)

bianaSession.output_user_entity_set_network(  
											  user_entity_set_id = 'SmartAS_entitySet', 
										 	  out_method = open('Results/allInteractions.tsv','w').write, 
											  node_attributes = ["ensembl","method_id"],
											  participant_attributes = [],
											  relation_attributes = ['psimi_name', 'Pubmed'],
											  allowed_relation_types = 'all',
											  include_participant_tags = False,
											  include_relation_tags = False,
											  include_relation_ids = True,
											  include_participant_ids = True,
											  include_relation_type = True,
											  include_relation_sources = True,
											  output_1_value_per_attribute = True,
											  output_format = 'tabulated',
											  only_selected = False
											)

r('load("SmartAS.RData")')
r('nodeDetails <- read.delim("~/SmartAS/nodeDetails_onlyUnique.tsv", quote="")')
r('write.table( nodeDetails[ order(-nodeDetails$Degree) ], paste(wd, "/Results/candidatesInteractions.sorted.tsv", sep=""), sep="\t", row.names=F, col.names=F)')

iLoopsPairs = open("Results/candidateList.top.lst", "w")
with open("Results/candidateList.lst", "r") as candidates:
	for line in candidates:
		elements = line.split("\t")
		candidate1 = elements[1]
		candidate2 = elements[2]

		hubCandidate = False
	
		with open("Results/candidatesInteractions.sorted.tsv", "r") as nodes:
			for count in range(10):
				line = nodes.readline()
				if line.find(candidate1) != -1 or line.find(candidate2) != -1:
					hubCandidate = True
					break

		if hubCandidate:
			iLoopsPairs.write(candidate1 + "\t" + candidate2 + "\n")

iLoopsPairs.close()