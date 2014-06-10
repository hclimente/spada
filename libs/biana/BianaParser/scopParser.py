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
from sets import *
import os
from biana.BianaObjects.PDB import PDBFragment

class SCOPParser(BianaParser):
    """
    SCOP Parser Class
    """

    name = "scop"
    description = "This program fills up tables in database biana related with SCOP"
    external_entity_definition = "A external entity represents a SCOP entity (fold, class, domain,...)"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "Structural Clasification Of Proteins",
                             default_script_name = "scopParser.py",
                             default_script_description = "This program fills up tables in database biana related to SCOP",
                             additional_compulsory_arguments = [])
        self.default_eE_attribute = "scop"
	#self.is_promiscuous = True


    def parse_database(self):
        """
        Method that implements the specific operations of scop parser
        """

        self.biana_access.add_valid_external_entity_attribute_type( name = "SCOP_Category",
                                                                    data_type = "ENUM(\"class\",\"fold\",\"superfamily\",\"family\",\"domain\")",
                                                                    category = "eE attribute" )


        # IMPORTANT: As we have added new types and attributes that are not in the default BIANA distribution, we must execute the follwing command:
        self.biana_access.refresh_database_information()


        def new_list():
            return []

        categories = {"cl":"class",
                      "cf":"fold",
                      "sf":"superfamily",
                      "fa":"family",
                      "dm":"domain"}

        number_of_lines = 0

        #
        # Reading external DB "scop" and inserting its data into Biana DB
        # 
        
        #d1uvya_ 1uvy    A:      a.1.1.1 100068  cl=46456,cf=46457,sf=46458,fa=46459,dm=46460,sp=46461,px=100068
	cl_regex = re.compile("cl=(\d+)")
	cf_regex = re.compile("cf=(\d+)")
	sf_regex = re.compile("sf=(\d+)")
	fa_regex = re.compile("fa=(\d+)")
	dm_regex = re.compile("dm=(\d+)")
	sp_regex = re.compile("sp=(\d+)")
	px_regex = re.compile("px=(\d+)")
	range_regex = re.compile("(\w+):(\S*)")
        tax_regex = re.compile("\[TaxId:\s(\d+)\]")

	domains_dict = {}
	#domains_to_id_dict = {}
	#hierarchy_dict = {"cl":{},"cf":{},"sf":{},"fa":{},"dm":{},"sp":{}}
	hierarchy_dict = {"cf":{},"sf":{},"fa":{},"dm":{}}
	descriptions_dict = {"cl":{},"cf":{},"sf":{},"fa":{},"dm":{},"sp":{},"px":{}}
        sp_dict = {}

        scop_entry_to_eE_id = {}

        if not self.input_file.endswith(os.sep):
            self.input_file += os.sep

        scop_dir_cla_fd = file(self.input_file+"dir.cla.scop.txt_"+self.sourcedb_version.replace("\"",""),'r')

        for line in scop_dir_cla_fd:

	    if line.startswith("#"):
		continue
            
	    line_fields = line.strip().split()   # line_fields[0] is complete pdb entry
                                       		 # line_fields[1] is pdb code
	                                         # line_fields[2] is pdb chain follow by : XX-YY (optional)
        	                                 # line_fields[3] is ???
                	                         # line_fields[4] is ???
                        	                 # line_fields[5] is comma-separated codes

            if len(line_fields) != 6:
                # skip incomplete lines
                print "skipping..."
                continue

	    pdb_code = line_fields[1]

	    #problem: a domain can be in different chains...
	    #m = range_regex.search(line_fields[2])
	    #if m:
	    #	chain = m.group(1)
	    #	range = m.group(2)
	    range = line_fields[2]

	    cl = cl_regex.search(line_fields[5]).group(1)
            cf = cf_regex.search(line_fields[5]).group(1)
            sf = sf_regex.search(line_fields[5]).group(1)
            fa = fa_regex.search(line_fields[5]).group(1)
            dm = dm_regex.search(line_fields[5]).group(1)
            sp = sp_regex.search(line_fields[5]).group(1)
                                        
	    hierarchy_dict["cf"][cf] = cl
            hierarchy_dict["sf"][sf] = cf
	    hierarchy_dict["fa"][fa] = sf
            hierarchy_dict["dm"][dm] = fa

            #sp_dict.setdefault(dm,new_list()).append(sp)
            sp_dict.setdefault(dm,Set(new_list())).add(sp)
            domains_dict.setdefault(dm,new_list()).append((pdb_code,range))
            #domains_to_id_dict[dm] = line_fields[3]

        scop_dir_cla_fd.close()

        #print len(hierarchy_dict), hierarchy_dict
        #print len(sp_dict), sp_dict
        #print len(domains_dict), domains_dict

        scop_des_fd = file(self.input_file+"dir.des.scop.txt_"+self.sourcedb_version.replace("\"",""),'r')

        for line in scop_des_fd:

            if line.startswith("#"):
                continue

	    line_fields = line.strip().split("\t")
            descriptions_dict[line_fields[1]][line_fields[0]] = line_fields[4]

        scop_des_fd.close()

        #print descriptions_dict

        for current_category in descriptions_dict:
            if current_category!="px" and current_category!="sp":
                for current_scop_entry in descriptions_dict[current_category]:
                    eE = ExternalEntity( source_database = self.database, type = "SCOPElement" )
                    eE.add_attribute( ExternalEntityAttribute(attribute_identifier="SCOP", value = current_scop_entry, type="unique") )
                    eE.add_attribute(ExternalEntityAttribute( attribute_identifier="SCOP_Category", value = categories[current_category], type="unique" ))
                    if current_category == "dm":
                        for current_pdb in domains_dict[current_scop_entry]:

                            fragments = PDBFragment.fragment_parser( fragment_str = current_pdb[1], separator = "," )

                            chain = fragments[0].chain  # Only takes the chain of the first fragment (all fragments should be on the same chain!). 
                                                        # If different fragments belong to different chains, it is not taken into account

                            additional_fields =  { "pdb_range": current_pdb[1] }
                            
                            if chain is not None:
                                additional_fields["chain"] = chain 

                            eE.add_attribute(ExternalEntityAttribute(attribute_identifier = "pdb", value=current_pdb[0],
                                                                     additional_fields = additional_fields,
								     type="cross-reference"))
                            eE.add_attribute(ExternalEntityAttribute(attribute_identifier="description", value = descriptions_dict[current_category][current_scop_entry]))
                            
                        #[ eE.add_attribute(ExternalEntityAttribute(attribute_identifier="taxid", value = current_sp_id)) for current_sp_id in sp_dict[current_scop_entry] ]

                        
                        for current_sp_id in sp_dict[current_scop_entry]:
                            
                            m = tax_regex.search(descriptions_dict["sp"][current_sp_id])

                            if m:
                                eE.add_attribute( ExternalEntityAttribute( attribute_identifier="taxid", value = m.group(1), type = "cross-reference" ) )

                            else:
                                print current_sp_id
                        
                        #[ eE.add_attribute(ExternalEntityAttribute(attribute_identifier="taxid", value = tax_regex.search(descriptions_dict["sp"][current_sp_id]).group(1))) for current_sp_id in sp_dict[current_scop_entry] ]

                    self.biana_access.insert_new_external_entity(eE)

                    scop_entry_to_eE_id[current_scop_entry] = eE.get_id()

        
        ontology = Ontology( source_database = self.database, linkedAttribute="scop", name="scop", descriptionAttribute="description", levelAttribute="SCOPCategory" )
        

        for current_dm_scop_element in hierarchy_dict["dm"]:
            ontology.add_element( ontologyElementID = scop_entry_to_eE_id[current_dm_scop_element],
                                 isA = [scop_entry_to_eE_id[hierarchy_dict["dm"][current_dm_scop_element]]] )
        for current_fa_scop_element in hierarchy_dict["fa"]:
            ontology.add_element( ontologyElementID = scop_entry_to_eE_id[current_fa_scop_element],
                                 isA = [scop_entry_to_eE_id[hierarchy_dict["fa"][current_fa_scop_element]]] )
        for current_sf_scop_element in hierarchy_dict["sf"]:
            ontology.add_element( ontologyElementID = scop_entry_to_eE_id[current_sf_scop_element],
                                 isA = [scop_entry_to_eE_id[hierarchy_dict["sf"][current_sf_scop_element]]] )
        for current_cf_scop_element in hierarchy_dict["cf"]:
            ontology.add_element( ontologyElementID = scop_entry_to_eE_id[current_cf_scop_element],
                                 isA = [scop_entry_to_eE_id[hierarchy_dict["cf"][current_cf_scop_element]]] )
        for current_cl in descriptions_dict["cl"]:
            ontology.add_element( ontologyElementID = scop_entry_to_eE_id[current_cl] )
                            
        self.biana_access.insert_new_external_entity( ontology )
