#!/soft/devel/python-2.7/bin/python

from biana import *
from biana.utilities import identifier_utilities
#from rpy2.robjects import r
#import sets

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
list_input_identifiers = identifier_utilities.read_identifier_list_from_file(file_name = "ENSTranscripts.lst", id_type = "ensembl")
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
												#out_method = open('tumor_ints_All_nodes_details.txt','w').write, 
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
												#out_method = open('tumor_ints_All_nodes_details_only_unique.txt','w').write, 
												attributes = ["ensembl","uniprotaccession","uniprotentry"], 
												include_level_info = True,
												include_degree_info = True,
												level = None,
												only_selected = False, 
												output_format = 'tabulated', 
												include_tags_info = False,
												include_tags_linkage_degree_info = [], 
												#output_only_unique_values = True,
												output_1_value_per_attribute = False
											)

bianaSession.output_user_entity_set_network(  
											  user_entity_set_id = 'SmartAS_entitySet', 
										 	  out_method = open('Results/allInteractions.tsv','w').write, 
										 	  #out_method = open('tumor_ints_All_nodes_network.txt','w').write, 
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

#r('load("SmartAS.RData")')
#r('nodeDetails <- read.delim("~/SmartAS/nodeDetails_onlyUnique.tsv", quote="")')
#r('mask <- nodeDetails$Level == 0')
#r('write.table( nodeDetails[mask,], paste(wd, "/Results/nodeDetails.tsv", sep=""), sep="\t", row.names=F)')

candidates = open("Results/candidateList.lst", "r")
for line in candidates:
	elements = line.split("\t", )
	transcript1 = elements[1]
	transcript2 = elements[2]

	nodes = open("Results/nodeDetails.tsv", "r")
	for node in nodes:
		if node.find(transcript1) != -1:
			print(transcript1 + " " + node)
			break
		if node.find(transcript2) != -1:
			print(transcript2 + " " + node)
			break