#!/soft/devel/python-2.7/bin/python

from include.biana import *
from include.biana.utilities import identifier_utilities
import csv
import sys

out = sys.argv[1]
inputType = ""
if sys.argv[2] == "GENCODE":
	inputType = "ensembl"
elif sys.argv[2] == "TCGA":
	inputType = "geneid"

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

with open("Results/" + out + "/candidateList.tsv", "r") as candidates:
	for line in candidates:
		fields = line.strip().split("\t")
		if inputType == "ensembl":
			list_input_identifiers.append(("ensembl", fields[2]))
			list_input_identifiers.append(("ensembl", fields[3]))
		elif inputType == "geneid":
			nameComponents = fields[1].split("|")
			if len(nameComponents) == 2:
				list_input_identifiers.append(("geneid", nameComponents[1] ))

list_input_restriction_identifiers = []
list_input_negative_restriction_identifiers = []

# Create The Set for Y2H-detected interactions
userEntity_Y2H = bianaSession.create_new_user_entity_set( 
																		identifier_description_list = list_input_identifiers, 
																		attribute_restriction_list = list_input_restriction_identifiers,
																		negative_attribute_restriction_list = list_input_negative_restriction_identifiers,					
																		id_type = 'embedded', 
																		new_user_entity_set_id = "entitySet_Y2H"
																	)

# Create Network Commands
bianaSession.create_network(
								user_entity_set_id = 'entitySet_Y2H', 
								level = 1, 
								include_relations_last_level = False, 
								relation_type_list = ["interaction"], 
								relation_attribute_restriction_list = [("method_id", 18)], 
								use_self_relations = True,
								expansion_attribute_list = [],
								expansion_relation_type_list = [], 
								expansion_level = 2, 
								attribute_network_attribute_list = [], 
								group_relation_type_list = []
							)

# Create The Set for Lumier-detected interactions
userEntity_Lumier = bianaSession.create_new_user_entity_set( 
																			identifier_description_list = list_input_identifiers, 
																			attribute_restriction_list = list_input_restriction_identifiers,
																			negative_attribute_restriction_list = list_input_negative_restriction_identifiers,					
																			id_type = 'embedded', 
																			new_user_entity_set_id = "entitySet_Lumier"
																		)

# Create Network Commands
bianaSession.create_network(
								user_entity_set_id = 'entitySet_Lumier', 
								level = 1, 
								include_relations_last_level = False, 
								relation_type_list = ["interaction"], 
								relation_attribute_restriction_list = [("method_id", 696)], 
								use_self_relations = True,
								expansion_attribute_list = [],
								expansion_relation_type_list = [], 
								expansion_level = 2, 
								attribute_network_attribute_list = [], 
								group_relation_type_list = []
							)

# Make union of the previous
userEntity_union = bianaSession.get_union_of_user_entity_set_list(
																			[userEntity_Y2H, userEntity_Lumier], 
																			include_relations=True, 
																			new_user_entity_set_id="SmartAS_entitySet"
																		)

# Output Commands
user_entities_to_print = set(userEntity_union.get_user_entity_ids(level=0))
bianaSession.select_user_entities_from_user_entity_set(
														user_entity_set_id = 'SmartAS_entitySet', 
                                                		user_entity_id_list = user_entities_to_print, 
                                                		clear_previous_selection = True
                                                	  )
bianaSession.output_user_entity_set_details(
												user_entity_set_id = 'SmartAS_entitySet', 
												out_method = open("Results/" + out + "/candidateInteractions.tsv",'w').write, 
												attributes = [inputType,"uniprotaccession","uniprotentry"], 
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
txInfo = {}

with open('Data/Databases/Intogen.tsv','r') as Intogen_r:
	print("\t* Reading known driver genes from IntOGen.")
	Intogen = csv.reader(Intogen_r, delimiter='\t')
	for row in Intogen:
		if inputType == "ensembl":
			intogenDrivers.add(row[0])
		elif inputType == "geneid":
			intogenDrivers.add(row[1])
with open("Data/Databases/compilationTable.tsv", "r") as compilationTable:
	print("\t* Reading information available in bibliography.")
	for line in compilationTable:
		splitted = line.strip().split("\t")
		if line.find("yes") != -1:
			articleCompilation[splitted[1]] = splitted[3:15]
		else:
			articleCompilation[splitted[1]] = []
			for i in range(1,12):
				articleCompilation[splitted[1]].append("-")
			articleCompilation[splitted[1]].append(splitted[14])

with open("Data/TCGA/knownGene.txt", "r") as TX_INFO:
	print("\t* Calculating transcript differences.")
	for line in TX_INFO:
		elements = line.strip().split("\t")
		name		= elements[0] 
		chrom		= elements[1] 
		strand		= elements[2] 
		txStart		= int(elements[3])
		txEnd		= int(elements[4])
		cdsStart	= int(elements[5])
		cdsEnd		= int(elements[6])
		exonCount	= int(elements[7])
		exonStarts	= map(int, filter(None, elements[8].split(",") ) )
		exonEnds	= map(int, filter(None, elements[9].split(",") ) )
		proteinID	= elements[10]
		alignID		= elements[11]

		txInfo[name] = {}
		txInfo[name]["txStart"] 	= txStart
		txInfo[name]["txEnd"] 		= txEnd
		txInfo[name]["cdsStart"] 	= cdsStart
		txInfo[name]["cdsEnd"] 		= cdsEnd
		txInfo[name]["exonStarts"] 	= exonStarts
		txInfo[name]["exonEnds"] 	= exonEnds

candidateList = []
with open("Results/" + out + "/candidateList.tsv", "r") as candidates:
	print("\t* Checking coincidences.")
	for line in candidates:
		aCandidate = {}
		aCandidate["Basic"] = line.strip()

		#Data from external sources doesn't have number of replicates
		if len(aCandidate["Basic"].split("\t")) == 5:
			aCandidate["Basic"] += "\t-1"
		
		aCandidate["UTR Change"] = "No"
		aCandidate["CDS Change"] = "No"
		aCandidate["CDS?"] = "Yes"
		aCandidate["Hub"] = "na"
		aCandidate["IntOGen"] = "no"
		aCandidate["Articles"] = []
		
		splitted = aCandidate["Basic"].split("\t")
		name = splitted[0]
		gene = splitted[1]
		candidate1 = splitted[2]
		candidate2 = splitted[3]
		
		#BIANA info
		with open("Results/" + out + "/candidateInteractions.tsv", "r") as nodes:
			for line in nodes:
				ids = set(line.split("\t")[3].split())
				if inputType == "ensembl" and (candidate1 in ids or candidate2 in ids): 
					aCandidate["Hub"] = (line.split("\t"))[2]
					break
				elif inputType == "geneid":
					nameComponents = gene.split("|")
					if len(nameComponents) == 2 and nameComponents[1] in ids:
						aCandidate["Hub"] = (line.split("\t"))[2]
						break
		
		#Intogen info
		if (inputType == "ensembl" and gene in intogenDrivers) or (inputType == "geneid" and name in intogenDrivers):
			aCandidate["IntOGen"] = "yes"
		
		#Journals info
		if name in articleCompilation.keys():
			for info in articleCompilation[name]:
				aCandidate["Articles"].append(info)
		else:
			aCandidate["Articles"] = ["-","-","-","-","-","-","-","-","-","-","-","NA"]

		#CDS and UTR change
		if txInfo[candidate1]["cdsStart"] == txInfo[candidate1]["cdsEnd"] or txInfo[candidate2]["cdsStart"] == txInfo[candidate2]["cdsEnd"]:
			aCandidate["CDS?"] = "No"

		if txInfo[candidate1]["exonStarts"] != txInfo[candidate2]["exonStarts"] or txInfo[candidate1]["exonEnds"] != txInfo[candidate2]["exonEnds"]:
			if txInfo[candidate1]["cdsStart"] != txInfo[candidate2]["cdsStart"] or txInfo[candidate1]["cdsEnd"] != txInfo[candidate2]["cdsEnd"]:
				aCandidate["CDS Change"] = "Yes"
			if txInfo[candidate1]["txStart"] != txInfo[candidate2]["txStart"] or txInfo[candidate1]["txEnd"] != txInfo[candidate2]["txEnd"]:
				aCandidate["UTR Change"] = "Yes"

		#Check if each of the exons at the CDS are the same.
		for anExonStart_1, anExonEnd_1 in zip(txInfo[candidate1]["exonStarts"], txInfo[candidate1]["exonEnds"] ):
			exonMatch = False
			
			for anExonStart_2, anExonEnd_2 in zip(txInfo[candidate2]["exonStarts"], txInfo[candidate2]["exonEnds"] ):
				if anExonStart_1 == anExonStart_2 and anExonEnd_1 == anExonEnd_2:
					exonMatch = True
					break
			
			if not exonMatch:
				if anExonStart_1 < txInfo[candidate1]["cdsStart"] or anExonStart_1 > txInfo[candidate1]["cdsEnd"]:
					aCandidate["UTR Change"] = "Yes"
				else:
					aCandidate["CDS Change"] = "Yes"

		candidateList.append(aCandidate)

with open("Results/" + out + "/candidateList.top.tsv", "w") as CANDIDATES_TOP:
	CANDIDATES_TOP.write("hugo_id\tGene\tTranscript_normal\tTranscript_tumor\tReplicates\tKnown PPI\tCDS?\tUTR Change\tCDS Change\tIntOGen\t")
	CANDIDATES_TOP.write("baltz_a\tcastello_a\tkwon_a\tgonzalez_a\tbrosseau_a\tvogelstein_a\than_a\tjuan_ap\tjuan_pr\tbiomart_a\tcosmic_a\treactome\n")
	for aCandidate in candidateList:
		CANDIDATES_TOP.write(aCandidate["Basic"] + "\t" + aCandidate["Hub"] + "\t" + aCandidate["CDS?"] + "\t")
		CANDIDATES_TOP.write(aCandidate["UTR Change"] + "\t" + aCandidate["CDS Change"] + "\t" + aCandidate["IntOGen"])
		
		for info in aCandidate["Articles"]:
			CANDIDATES_TOP.write("\t" + info)
		CANDIDATES_TOP.write("\n")