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

#########################################################################
# BIANA output modification Utility Library 
# Methods 
#   to get primary uniprot accession ids from tsv formatted node output file
#   to map uniprots to genes (vice verca) in tsv node output
#   to filter interaction network based on method type
#   to select/convert sequence from BIANA tsv output to fasta
#
# eg 22/06/2009
#########################################################################

from biana.utilities import TsvReader
from biana.utilities import FastaWriter

def get_primary_uniprot_accessions(input_file):
    """
	Returns a list of primary uniprot accession ids from tsv formatted node output file for each user entity
    """
    import re
    swissprot_exp = re.compile("\w\d\d\d\d\d")
    reader = TsvReader.TsvReader(input_file)
    columns, id_to_vals = reader.read(fields_to_include = ["User Entity ID", "uniprotaccession"], overwrite_keys = True)
    #print columns, id_to_vals.items()[0]
    ids = set()
    for ueid, vals in id_to_vals.iteritems():
	words = vals[columns["uniprotaccession"]].split(",")
	local_ids = []
	containedP = None
	containedO = None
	for id in words:
	    id = id.strip()
	    if re.match(swissprot_exp, id):
		#ids.add(id)
		if id.startswith("P"):
		    containedP = id
		if id.startswith("O"):
		    containedO = id
	    local_ids.append(id)
	if len(local_ids) > 1:
	    if containedP is not None:
		ids.add(containedP)
	    elif containedO is not None:
		ids.add(containedO)
	    else:
		#print local_ids
		ids.add(local_ids[0])
	#elif len(local_ids) == 0:
	#    print "no swissprot related primary accession found:", ueid
	else:
	    #print local_ids
	    ids.add(local_ids[0])
    return ids


def convert_sequence_in_tsv_to_fasta(tsv_file_name, fasta_file_name, select_representative_sequence=False):
    reader = TsvReader.TsvReader(tsv_file_name)
    columns, id_to_vals = reader.read(fields_to_include = ["user entity id", "proteinsequence"]) # If overwrite_keys = False stores list of list as values
    #print columns
    #print id_to_vals
    file_out = open(fasta_file_name, "w")
    fwriter = FastaWriter.FastaWriter(out_method=file_out.write, one_line_per_sequence=True) 
    for id, vals in id_to_vals.iteritems():
    #for id, vals_list in id_to_vals.iteritems():
	#for vals in vals_list:
	#file_out.write(">%s\n%s\n" % (id, vals[columns["proteinsequence"]])) 
	if select_representative_sequence:
	    seqToCount = {}
	    for val in vals[columns["proteinsequence"]].split(","):
		val = val.strip()
		seqToCount[val] = seqToCount.setdefault(val, 0) + 1
	    seqCounts = seqToCount.items()
	    seqCounts.sort(lambda x,y: cmp(x[1],y[1]))
	    seqCounts.reverse()
	    #seq = seqCounts[0][0]
	    #if id == "38347":
	    #	print seqCounts
	    nCount = seqCounts[0][1]
	    seq_selected = seqCounts[0][0]
	    seq_len = 0
	    for seq, n in seqCounts:
		if n < nCount:
		    break
		if len(seq) > seq_len:
		    seq_selected = seq
		    seq_len = len(seq)
		    #print id,
	    #print
	else:
	    seq_selected = vals[columns["proteinsequence"]]
	if seq_selected == "-":
	    continue
	fwriter.output_sequence(id, seq_selected)
    file_out.close()
    #file_out = open(fasta_file_name+"2", "w")
    #reader.process(file_out.write, fields_to_include=["user entity id", "proteinsequence"])
    #file_out.close()
    return

def get_attribute_to_attribute_mapping(tsv_file_name, attribute, to_attribute, in_value_separator = ",", out_file_name=None, keys_to_include=None):
    reader = TsvReader.TsvReader(tsv_file_name)
    columns, id_to_vals = reader.read(fields_to_include = [attribute, to_attribute], keys_to_include = keys_to_include) 
    #print columns
    #print id_to_vals

    if out_file_name is not None:
	file_out = open(out_file_name, "w")
	reader.process(file_out.write, fields_to_include=[attribute, to_attribute]) 
	file_out.close()
	return
    attribute_to_nodes = {}
    node_to_attributes = {}
    for id, vals in id_to_vals.iteritems():
    #for id, vals_list in id_to_vals.iteritems():
	#for vals in vals_list:
	#file_out.write("%s\t%s\n" % (id, vals[columns["genesymbol"]])) 
	words = vals[columns[to_attribute]].split(in_value_separator)
	[ attribute_to_nodes.setdefault(val.strip(), set()).add(id) for val in words ] 
	[ node_to_attributes.setdefault(id, set()).add(val.strip()) for val in words ]  
    #file_out.close()
    return node_to_attributes, attribute_to_nodes

def get_user_entity_attribute_mapping(tsv_file_name, attribute, out_file_name=None):
    return get_attribute_to_attribute_mapping(tsv_file_name, "user entity id", attribute, out_file_name)

def filter_network_by_interaction_type(network_attribute_file_name, network_out_file_name, interaction_type="y2h", reverse_selection=False):
    """
	interaction_type: "y2h" | "tap" 
	reverse_selection: True => no tap / no y2h
    """
    y2h = set([18,397,727,728,437,398,399])
    tap = set([4,96,676,729,19,6,7,858,59,109])
    other = set([114, 441, 492, 493, 802]) # xray, sga, in vitro, in vivo, enhancement
    if interaction_type=="y2h":
	valid_ids = y2h
    elif interaction_type=="tap":
	valid_ids = tap
    else:
	print "Unsupported interaction type", interaction_type
	return
    f = open(network_attribute_file_name)
    valid_interactions = set()
    for line in f:
	words = line[:-1].split()
	if len(words) == 1 and line.startswith("method_id"):
	    continue
	if words[3] != "=":
	    print "format error", line
	    continue
	if reverse_selection:
	    if int(words[4]) not in valid_ids:
		valid_interactions.add(" ".join([words[0], words[1], words[2]]))
	else:
	    if int(words[4]) in valid_ids:
		valid_interactions.add(" ".join([words[0], words[1], words[2]]))
    f.close()
    f_out = open(network_out_file_name, 'w')
    for i in valid_interactions:
	f_out.write(i+"\n")
    f_out.close()
    return

def get_tap_method_id_querry():
    return "select * from ExternalEntityOntology_isA A, externalEntitypsimi_name N, key_attribute_3 K, key_attribute_3 K2 where K.externalEntityID=A.is_a and A.externalEntityID=N.externalEntityID and N.externalEntityID=K2.externalEntityID and K.value IN (4,96,676,19,6,7,729,858,59,109)"
    # 59 (gst pull down)
    # 109 (tap tag coimmunoprecipitation)
    # considering 25 (copurification, obsolute because not specific, can be tap related or not..)
    # 46, 63 (interaction prediction) also 58, 101, 105

def get_y2h_method_id_querry():
    return "select * from ExternalEntityOntology_isA A, externalEntitypsimi_name N, key_attribute_3 K, key_attribute_3 K2 where K.externalEntityID=A.is_a and A.externalEntityID=N.externalEntityID and N.externalEntityID=K2.externalEntityID and K.value IN (18,397,727,728,437,398,399)"


if __name__ == "__main__":
    node_to_genes, gene_to_nodes = get_user_entity_gene_mapping(tsv_file_name="/home/emre/arastirma/data/human_interactome_biana/test_nodes.tsv", out_file_name="test.txt")
    print node_to_genes, gene_to_nodes

