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


intogenDrivers = set()
articleCompilation = {}

with open('Data/Intogen.tsv','r') as Intogen_r:
	print("\t* Reading known driver genes from IntOGen.")
	Intogen = csv.reader(Intogen_r, delimiter='\t')
	for row in Intogen:
		intogenDrivers.add(row[0])
with open("Data/compilationTable.tsv", "r") as compilationTable:
	print("\t* Reading information available in bibliography.")
	for line in compilationTable:
		splitted = line.strip().split("\t")
		if line.find("yes") != -1:
			articleCompilation[splitted[1]] = splitted[3:15]
		else:
			articleCompilation[splitted[1]] = []
			for i in range(1,12):
				articleCompilation[splitted[1]].append("")
			articleCompilation[splitted[1]].append(splitted[14])


candidateList = []
with open("Results/candidateList.tsv", "r") as candidates:
	print("\t* Checking coincidences.")
	for line in candidates:
		aCandidate = {}
		aCandidate["Basic"] = line.strip()
		aCandidate["Hub"] = "na"
		aCandidate["IntOGen"] = "no"
		aCandidate["Articles"] = []
		
		splitted = aCandidate["Basic"].split("\t")
		name = splitted[0]
		gene = splitted[1]
		candidate1 = splitted[2]
		candidate2 = splitted[3]
		
		with open("Results/candidateInteractions.tsv", "r") as nodes:
			for line in nodes:
				if line.find(candidate1) != -1 or line.find(candidate2) != -1:
					aCandidate["Hub"] = (line.split("\t"))[2]
					break
		
		if gene in intogenDrivers:
			aCandidate["IntOGen"] = "yes"
		if name in articleCompilation.keys():
			for info in articleCompilation[name]:
				aCandidate["Articles"].append(info)

		candidateList.append(aCandidate)

with open("Results/candidateList.top.tsv", "w") as candidates:
	candidates.write("hugo_id\tENSG\tENST_normal\tENST_tumor\tKnown PPI\tIntOGen\t")
	candidates.write("baltz_a\tcastello_a\tkwon_a\tgonzalez_a\tbrosseau_a\tvogelstein_a\than_a\tjuan_ap\tjuan_pr\tbiomart_a\tcosmic_a\treactome\n")
	for aCandidate in candidateList:
		candidates.write(aCandidate["Basic"] + "\t" + aCandidate["Hub"] + "\t" + aCandidate["IntOGen"])
		for info in aCandidate["Articles"]:
			candidates.write("\t" + info)
		candidates.write("\n")