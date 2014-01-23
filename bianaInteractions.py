from biana import *
import sets

from biana.utilities import identifier_utilities



# START BIANA SESSION

biana_session_object = create_new_session(sessionID='session_ID', 
				dbname='BIANA_MARCH_2013',
				dbhost='ben-yehuda',
				dbuser='biana_user',
				dbpassword='biana_password', 
				unification_protocol='uniprot_geneID_seqtax')


# CREATE A LIST WITH ALL THE SEED IDENTIFIERS

list_input_identifiers = []
list_input_identifiers.extend( identifier_utilities.read_identifier_list_from_file(file_name = "ensgenes.txt", id_type="ensembl")  )
list_input_restriction_identifiers = []
list_input_negative_restriction_identifiers = []


# CREATE THE SET

user_entity_set_object = biana_session_object.create_new_user_entity_set( identifier_description_list = list_input_identifiers, 
							attribute_restriction_list = list_input_restriction_identifiers,
							negative_attribute_restriction_list = list_input_negative_restriction_identifiers,					
							id_type='embedded', 
							new_user_entity_set_id="my_user_entity_set")


# CREATE NETWORK COMMANDS

biana_session_object.create_network(user_entity_set_id = 'my_user_entity_set', 
					level = 1, 
					include_relations_last_level = False, 
					relation_type_list = ["interaction"], 
					relation_attribute_restriction_list = [], 
					use_self_relations = True,
					expansion_attribute_list = [],
					expansion_relation_type_list = [], 
					expansion_level = 1, 
					attribute_network_attribute_list = [], 
					group_relation_type_list = [])


# OUTPUT COMMANDS

biana_session_object.output_user_entity_set_details(user_entity_set_id = 'my_user_entity_set', 
						out_method = open('tumor_ints_All_nodes_details.txt','w').write, 
						attributes = ["ensembl","uniprotaccession","uniprotentry"], 
						include_level_info = True,
						include_degree_info = True,
						level = None,
						only_selected = False, 
						output_format = 'tabulated', 
						include_tags_info = False,
						include_tags_linkage_degree_info = [], 
						output_1_value_per_attribute = False)
biana_session_object.output_user_entity_set_details(user_entity_set_id = 'my_user_entity_set', 
						out_method = open('tumor_ints_All_nodes_details_only_unique.txt','w').write, 
						attributes = ["ensembl","uniprotaccession","uniprotentry"], 
						include_level_info = True,
						include_degree_info = True,
						level = None,
						only_selected = False, 
						output_format = 'tabulated', 
						include_tags_info = False,
						include_tags_linkage_degree_info = [], 
						output_only_unique_values = True,
						output_1_value_per_attribute = False)
biana_session_object.output_user_entity_set_network(user_entity_set_id = 'my_user_entity_set', 
		             	        	  out_method = open('tumor_ints_All_nodes_network.txt','w').write, 
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
						  only_selected = False)name => 
type => 
tmp_name => 
error => 4
size => 0

name => ensgenes.txt
type => text/plain
tmp_name => /tmp/phpBNIKil
error => 0
size => 1391