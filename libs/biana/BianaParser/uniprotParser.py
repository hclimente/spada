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


class UniprotParser(BianaParser):
    """
    Uniprot Parser Class
    """

    name = "uniprot"
    description = "This file implements a program that fills up tables in database biana with information of uniprot databases"
    external_entity_definition = "A external entity represents a protein"
    external_entity_relations = ""

    def __init__(self):

        # Start with the default values

        BianaParser.__init__(self, default_db_description = "Uniprot database",
                             default_script_name = "uniprotParser.py",
                             default_script_description = UniprotParser.description,
                             additional_optional_arguments = [])
        self.default_eE_attribute = "uniprotaccession"
        self.initialize_input_file_descriptor()

    def parse_database(self):
        """
        Method that implements the specific operations of uniprot parser

        If executing in self.mode "tables", it is necessary to insert the tables here
        """
        
        protein_number=0


        # New entry regex
        new_regex = re.compile("^\/\/\s*$")

        # General regex
        id_regex = re.compile("^ID\s+(\S+)\s*")
        #ac_regex = re.compile("^AC\s+(\S+)\;\s*$")
        ac_regex = re.compile("^AC\s+(.+);\s*$")
        ac_version_regex = re.compile("sequence version (\d+)")
        de_regex = re.compile("^DE\s+(.+)\s*$")
        taxID_regex = re.compile("^OX\s+NCBI_TaxID=(\d+);")
        keyword_regex = re.compile("^KW\s+(.+);$")
        

        # GeneName regex
        geneName_regex = re.compile("^GN")
        gene_name_regex = re.compile("Name=([^;]+);")
        gene_orf_name_regex = re.compile("ORFNames=([^;]+);")
        gene_synonyms_regex = re.compile("Synonyms=([^;]+);")
        gene_orderedLocusNames = re.compile("OrderedLocusNames=([^;]+);")
        
        #Cross-references regular expressions
        cross_regex = re.compile("^DR")
        
        pfam_regex = re.compile("^DR\s+Pfam;\s*(\S+);")
        kegg_regex = re.compile("^DR\s+KEGG;\s*(\S+);")
        interpro_regex = re.compile("^DR\s+InterPro;\s*(\S+);")
        prosite_regex = re.compile("^DR\s+PROSITE;\s*(\S+);")
        prodom_regex = re.compile("^DR\s+ProDom;\s*(\S+);")
        mim_regex = re.compile("^DR\s+MIM;\s*(\S+);")
        pir_regex = re.compile("^DR\s+PIR;\s*(\S+);")
        prints_regex = re.compile("^DR\s+PRINTS;\s*(\S+);")
        #ensembl_regex = re.compile("^DR\s+Ensembl;\s*(\S+);")
        ensembl_regex = re.compile("^DR\s+Ensembl;((\s*\S+;)+)+")
        #embl_regex = re.compile("^DR\s+EMBL;\s*(\S+);")
	embl_regex = re.compile("^DR\s+EMBL;((\s*\S+;)+)+")
        geneID_regex = re.compile("^DR\s+GeneID;\s*(\S+);")
        go_regex = re.compile("^DR\s+GO;\s*GO\:(\d+);")
        refseq_regex = re.compile("^DR\s+RefSeq;\s*(\S+);")
        unigene_regex = re.compile("^DR\s+UniGene;\s*(\S+);")
        hgnc_regex = re.compile("^DR\s+HGNC;\s*HGNC\:(\d+);")
        pdb_regex = re.compile("^DR\s+PDB;\s*(\S+);.+;.+;(.+).")
        flybase_regex = re.compile("^DR\s+FlyBase;\s*(\S+);")
        mgi_regex = re.compile("^DR\s+MGI;\s*MGI:(\d+);")
        reactome_regex = re.compile("^DR\s*Reactome;\s*REACT_(\d+);")
        sgd_regex = re.compile("^DR\s+SGD;\s*(\w+);")
        

        tigr_regex = re.compile("^DR\s+TIGR\;\s+(.+)\;")
        #intact_regex = re.compile("^DR\s+IntAct") # It's the same as UniprotAccession?
        dip_regex = re.compile("^DR\s+DIP\;\s+DIP\:(.+)\;")
        cygd_regex = re.compile("^DR\s+CYGD\;\s+(.+)\;")
        #arrayexpress_regex = re.compile("")  # It's the same as UniprotAccession?
        WormPep_regex = re.compile("^DR\s+WormPep\;\s+(.+)\;\s*CE(\d+)\.\s*$")
        WormBase_regex = re.compile("^DR\s+WormBase\;\s*WBGene(\d+)\;\s*(.+)\.\s*$")
        rgd_regex = re.compile("^DR\s+RGD\;\s+(\d+)\;")
        
        #Sequence
        sequence_regex = re.compile("^\s+(.+)$")

        #Comments
        new_comment_regex = re.compile("^CC\s+\-\!\-")
        general_comment_regex = re.compile("^CC\s+(.+)$")
        subcellular_location_regex = re.compile("SUBCELLULAR LOCATION:\s*(.*)$")
        function_regex = re.compile("FUNCTION:\s*(.*)$")
        disease_regex = re.compile("DISEASE:\s*(.*)$")
        
        #Start first uniprotObject
        uniprotObject = ExternalEntity( source_database = self.database, type="protein" )


        description = []
        sequence = []
        comments = { "SubcellularLocation": [],
                     "Disease": [],
                     "Function": [] }

        actual_comment = None

        self.initialize_input_file_descriptor()

        uniprot_accession_list = []

        # START PARSING
        for line in self.input_file_fd:

            # New entry
            if new_regex.match(line):

                if uniprotObject is not None:
                    #add sequence
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="proteinSequence", value=ProteinSequence("".join(sequence)), type = "unique"))

                    #add description
                    if len(description)>0:
                        desc_str = " ".join(description)
                        uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="description", value=desc_str))

                    # Detect EC code in description
                        #if desc_str  != "":
                        #    enzymes = re.findall("\(EC[\s\=]*\d+\.\d+\.\d+\.\d+\)", desc_str)
                        #    for enzyme in enzymes:
                        #        enzyme = re.sub("[(^\(EC[\s\=]*)(\)$)]", "", enzyme).strip()
                        #        if enzyme != "":
                        #            uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="EC", value=enzyme, type="cross-reference"))
                        
                        if desc_str != "":
                            enzymes = re.findall("EC=([\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+)\;", desc_str)
                            for enzyme in enzymes:
                                uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="EC", value=enzyme, type="cross-reference"))

                    #add comments
                    if len(comments["Function"])>0:
                        uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="function", value= " ".join(comments["Function"]), type="cross-reference"))

                    if len(comments["Disease"])>0:
                        uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="disease", value = " ".join(comments["Disease"]), type="cross-reference"))

                    if len(comments["SubcellularLocation"])>0:
                        uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="subcellularLocation", value = " ".join(comments["SubcellularLocation"]), type="cross-reference"))

                    # restart variables
                    description = []
                    sequence = []
                    comments = { "SubcellularLocation": [],
                                 "Disease": [],
                                 "Function": [] }
                    actual_comment = None

                    # Insert
                    self.biana_access.insert_new_external_entity( externalEntity = uniprotObject )


                # Start new object
                uniprotObject = ExternalEntity( source_database = self.database, type="protein" )
                protein_number += 1

                sequence = []
                
                if self.time_control:
                    if protein_number%20000==0:
                        sys.stderr.write("%s proteins done in %s seconds\n" %(protein_number,time.time()-self.initial_time))


            m = id_regex.match(line)
            if m:
                uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="uniprotentry", value=m.group(1), type="unique"))
                continue

            m = ac_regex.match(line)
            if m:
                uniprot_accession_list.extend([ x.strip() for x in m.group(1).split(";") ])
                continue

            m = ac_version_regex.search(line)
            if m:
                #print uniprot_accession_list
                #[ uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="uniprotaccession", value=x, version=m.group(1), type="unique")) for x in uniprot_accession_list ]
		# First one is the primary accession the followings are previous accessions
		for i in xrange(len(uniprot_accession_list)):
		    x = uniprot_accession_list[i]
		    if i == 0:
			uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="uniprotaccession", value=x, version=m.group(1), type="unique")) 
		    else:
			uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="uniprotaccession", value=x, version=m.group(1), type="previous")) 
                uniprot_accession_list = []
                continue

            m = de_regex.match(line)
            if m:
                description.append( m.group(1) )
                continue

            m = taxID_regex.match(line)
            if m:
                uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="taxID", value=m.group(1), type = "unique"))
                continue

            m = keyword_regex.match(line)
            if m:
                [ uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="keyword", value=x, type="cross-reference")) for x in m.group(1).split(";") ]
                continue


            m = sequence_regex.match(line)
            if m:
                sequence.append( m.group(1).replace(" ","")  )

            # Gene
            m = geneName_regex.match(line)

            if m:
                m = gene_name_regex.search(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="geneSymbol", value=m.group(1),type="unique"))

                m = gene_orf_name_regex.search(line)
                if m:
                    [ uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="ORFName", value=x,type="alias")) for x in m.group(1).split(",") ]
                    
                m = gene_synonyms_regex.search(line)
                if m:
                    [ uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="geneSymbol", value=x, type="synonym")) for x in m.group(1).split(",") ]

                m = gene_orderedLocusNames.search(line)
                if m:
                    [ uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="OrderedLocusName", value=x, type="alias")) for x in m.group(1).split(",") ]


                continue


            # COMMENTS
            m = general_comment_regex.match(line)
            if m:

                if( new_comment_regex.match(line)):
                    actual_comment = None
                    m = subcellular_location_regex.search(line)
                    if m:
                        actual_comment = "SubcellularLocation"
                    else:
                        m = function_regex.search(line)
                        if m:
                            actual_comment = "Function"
                        else:
                            m = disease_regex.search(line)
                            if m:
                                actual_comment = "Disease"

                    if actual_comment is not None: 
                        comments[actual_comment].append(m.group(1))

                else:
                    if actual_comment is not None:
                        comments[actual_comment].append(m.group(1))


            # CROSS-REFERENCES
            if cross_regex.match(line):

                m = WormBase_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="WormBaseGeneID", 
                                                                        value=m.group(1),type="cross-reference"))
                    continue

                m = WormPep_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="WormBaseSequenceName",
                                                                        value=m.group(1),type="cross-reference"))
                    continue

                m = dip_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="DIP",
                                                                        value=m.group(1),type="cross-reference"))
                    continue

                m = tigr_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="tigr",
                                                                        value=m.group(1),type="cross-reference"))

                    continue

                m = cygd_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="cygd",
                                                                        value=m.group(1),type="cross-reference"))

                    continue

                
                m = rgd_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="rgd",
                                                                        value=m.group(1),type="cross-reference"))
                    continue
            
                m = pfam_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="pfam", value=m.group(1),type="cross-reference"))
                    continue

                m = kegg_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="kegggene", value=m.group(1),type="cross-reference"))
                    continue

                m = interpro_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="interpro", value=m.group(1),type="cross-reference"))
                    continue

                m = prosite_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="prosite", value=m.group(1),type="cross-reference"))
                    continue

                m = prodom_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="prodom", value=m.group(1), type="cross-reference"))
                    continue

                m = mim_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="mim", value=m.group(1), type="cross-reference"))
                    continue

                m = pir_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="pir", value=m.group(1), type="cross-reference"))
                    continue

                m = prints_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="prints", value=m.group(1), type="cross-reference"))
                    continue

                m = ensembl_regex.match(line)
                if m:
                    #uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="ensembl", value=m.group(1), type="cross-reference"))
		    words = m.group(1).split(";")
		    for w in words:
			w=w.strip()
			if w != "-" and w != "":
			    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="ensembl", value=w, type="cross-reference"))
                    continue

                m = embl_regex.match(line)
                if m:
                    #uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="accessionNumber", value=m.group(1), type="cross-reference"))
		    words = m.group(1).split(";")
		    for w in words:
			w=w.strip()
			if w != "-" and w != "":
			    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="accessionNumber", value=w, type="cross-reference"))
                    continue

                m = geneID_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="geneID", value=m.group(1), type="cross-reference"))
                    continue

                m = go_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="go", value=m.group(1), type="cross-reference"))
                    continue

                m = refseq_regex.match(line)
                if m:
                    rs = m.group(1).split('.')
                    if len(rs)==2:
                        uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="refseq", value=rs[0], version=rs[1], type="cross-reference"))
                    else:
                        print "Refseq %s has no version?" %m.group(1)
                    continue
            
                m = unigene_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="unigene", value = m.group(1), type="cross-reference"))
                    continue

                m = hgnc_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="hgnc", value=m.group(1),type="cross-reference"))
                    continue

                m = pdb_regex.match(line)
                if m:
                    pdb_code = m.group(1)

                    fragments = m.group(2).split(",")

                    for actual_frag in fragments:
                        m = re.search("\s*(.+)=(.+)\s*",actual_frag)
                        if m:
                            chains = m.group(1).split("/")
                            m = re.search("(\d+)-(\d+)",m.group(2))
                            if m:
                                range = "%s-%s" %(m.group(1),m.group(2))

                                [ uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="pdb", value=pdb_code, type = "cross-reference",
                                                                                      additional_fields = {"chain": x,
                                                                                                           "pdb_range": range })) for x in chains ]
                            else:
                                [ uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="pdb", value=pdb_code, type="cross-reference",
                                                                                      additional_fields = {"chain": x})) for x in chains ]
                                
                    continue

                m = flybase_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="flybase", value=m.group(1), type = "cross-reference"))
                    continue

                m = mgi_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="MGI", value=m.group(1),type="cross-reference"))
                    continue


                m = reactome_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="reactome", value = m.group(1), type="cross-reference"))
                    continue

                m = sgd_regex.match(line)
                if m:
                    uniprotObject.add_attribute(ExternalEntityAttribute(attribute_identifier="SGD", value=m.group(1), type="cross-reference"))

                    continue
