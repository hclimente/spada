#!/soft/devel/python-2.7/bin/python

from biana import *
from biana.utilities import identifier_utilities
from rpy2.robjects import r
import sys
import csv

top = 10
if(len(sys.argv) == 1):
	top = int(sys.argv[1])

print("\t* Prioritizing by known interactions of the genes (BIANA).")

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

with open("Results/candidateList.tsv", "r") as candidates:
	for line in candidates:
		elements = line.split("\t")
		list_input_identifiers.append(("ensembl", elements[2]))
		list_input_identifiers.append(("ensembl", elements[3].strip()))

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
								expansion_level = 2, 
								attribute_network_attribute_list = [], 
								group_relation_type_list = []
							)

# Output Commands

#bianaSession.output_user_entity_set_details(	
#												user_entity_set_id = 'SmartAS_entitySet', 
#												out_method = open('Results/candidateInteractions_extended.tsv','w').write, 
#												attributes = ["ensembl","uniprotaccession","uniprotentry"], 
#												include_level_info = True,
#												include_degree_info = True,
#												level = None,
#												only_selected = False, 
#												output_format = 'tabulated', 
#												include_tags_info = False,
#												include_tags_linkage_degree_info = [], 
#												output_1_value_per_attribute = False
#											)
#
#bianaSession.output_user_entity_set_network(  
#											  user_entity_set_id = 'SmartAS_entitySet', 
#										 	  out_method = open('Results/allInteractions.tsv','w').write, 
#											  node_attributes = ["ensembl","method_id"],
#											  participant_attributes = [],
#											  relation_attributes = ['psimi_name', 'Pubmed'],
#											  allowed_relation_types = 'all',
#											  include_participant_tags = False,
#											  include_relation_tags = False,
#											  include_relation_ids = True,
#											  include_participant_ids = True,
#											  include_relation_type = True,
#											  include_relation_sources = True,
#											  output_1_value_per_attribute = True,
#											  output_format = 'tabulated',
#											  only_selected = False
#											)

user_entities_to_print = set(user_entity_set_object.get_user_entity_ids(level=0))
bianaSession.select_user_entities_from_user_entity_set(
														user_entity_set_id = 'SmartAS_entitySet', 
                                                		user_entity_id_list = user_entities_to_print, 
                                                		clear_previous_selection = True
                                                	  )
bianaSession.output_user_entity_set_details(
												user_entity_set_id = 'SmartAS_entitySet', 
												out_method = open('Results/candidateInteractions.tsv','w').write, 
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


r('load("SmartAS.RData")')
r('nodeDetails <- read.delim("Results/candidateInteractions.tsv", quote="")')
r('write.table( nodeDetails[ order(-nodeDetails$Degree), ][1:10,], paste(wd, "Results/candidateInteractions.sorted.tsv", sep=""), sep="\t", row.names=F, col.names=F)')

print("\t* Prioritizing by known driver genes (IntOGen).")
intogenDrivers = set()
with open('Data/Intogen.tsv','r') as Intogen:
	Intogen = csv.reader(Intogen, delimiter='\t')
	for row in Intogen:
		intogenDrivers.add(row[0])

with open("Results/candidateList.top.tsv", "w") as iLoopsPairs:
	with open("Results/candidateList.tsv", "r") as candidates:
		for line in candidates:
			elements = line.split("\t")
			name = elements[0]
			gene = elements[1]
			candidate1 = elements[2]
			candidate2 = elements[3].strip()
			
			with open("Results/candidateInteractions.sorted.tsv", "r") as nodes:
				for line in nodes:
					if line.find(candidate1) != -1 or line.find(candidate2) != -1:
						iLoopsPairs.write(name + "\t" + gene + "\t" + candidate1 + "\t" + candidate2 + "\t" + "Hub" + "\n")
						break
			
			if gene in intogenDrivers:	
				iLoopsPairs.write(name + "\t" + gene + "\t" + candidate1 + "\t" + candidate2 + "\t" + "Driver gene" + "\n")