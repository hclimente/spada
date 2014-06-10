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

def read_identifier_list_from_file(file_name, id_type=None):
    """
    File either contains an identifier per line
    or (attribute_name, identifier) couple on each line

    - Designated for use from command line scripts
    ------
    file_name: name of the file identifiers will be read from
    identifier_description_list: list of identifiers (or (id_type, id_string) tuples in case id_type is "embedded")
    id_type: type of the identifiers in the file, if "embedded" then file contains "attribute_name" "identifier" (seperated by tab) instead of just identifier at each line, if "auto" certain regular expressions will be tried to matched for automatic detection
    """
    fileIdentifierList = open(file_name)
    line = fileIdentifierList.readline()
    identifier_description_list = []

    ## read input file
    while line:
	words = line[:-1].split("\t")
	if len(words) == 2 and id_type=="embedded":
	    identifierType = words[0].strip()
	    identifierString = words[1].strip() #.replace("'","")
	    identifier_description_list.append((identifierType, identifierString))
	elif len(words) == 1:
	    identifierString = words[0].strip() #.replace("'","")
	    if id_type is None:
		identifier_description_list.append(identifierString)
	    elif id_type == "auto":
		detected_id_type = detect_identifier_type(identifierString)
		if detected_id_type is None:
		    print "Warning: Can not auto-detect type of %s" % identifierString
		else:
		    identifier_description_list.append((detected_id_type, identifierString))
	    else:
		identifier_description_list.append((id_type, identifierString))
	else:
	    raise Exception("Unrecognized id_type: %s" % id_type)

	line = fileIdentifierList.readline()

    return identifier_description_list



def detect_identifier_type(id):
    """
    Return auto-detected type of the given id

    Works for ids with strict regular format (such as Q7Z4I7-1, RAW1_HUMAN, ENSP00000215832, AC1423234.1, IPI00157836, EC1.2.3.12, PF681234, ...)
    """
    import re

    re_uniprot_entry = re.compile("\w+_[A-Z]+$")
    re_uniprot_accession = re.compile("[A-Z]\d\w{3}\d(-\d{1,3}){0,1}$")
    re_ensembl = re.compile("ENS[A-Z]+\d{11}$")
    re_accession_number = re.compile("([A-Z]){1,2}\d{3,7}(.\d{1,3}){0,1}$")
    #re_ipi = re.compile("IPI\d{8}$")
    #re_ec = re.compile("EC\d{1,2}")
    #re_ordered_locus_name = re.compile("[A-Z]{1,3}\d{3,4}[W|C|](.\d{1,3}){0,1}$")

    detected_type = None
    if re_uniprot_accession.match(id):
	detected_type = "uniprotaccession"
    elif re_uniprot_entry.match(id):
	detected_type = "uniprotentry"
    #elif id.startswith("ENS"):
    elif re_ensembl.match(id):
	detected_type = "ensembl"
    #elif id.startswith("AC"):
    elif re_accession_number.match(id):
	detected_type = "accessionnumber"
    #elif re_orf_name.search(id):
    #    detected_type = "orfname"
    elif id.startswith("IPI"):
	detected_type = "ipi"
    elif id.startswith("EC"):
	detected_type = "ec"
    elif id.startswith("PF"):
	detected_type = "pfam"
    #else: # handled in the caller
    #    print "Warning: Undetermined type: %s" % id

    return detected_type


