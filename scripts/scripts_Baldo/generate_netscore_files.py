import sys
import argparse
import os
import re
import biana
try: from biana import *
except: sys.exit(10)


def main():

    options = parse_user_arguments()
    extract(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Generate bPPI network with seed genes",
        epilog      = "@oliva's lab 2013")
    parser.add_argument('-iseed','--seeds_input_file',dest='seed',action = 'store',default='input_seed',
                        help = 'Seeds Input file (default is input_seed)')
    parser.add_argument('-rseed','--restrict_network_to_seeds',dest='restricted_to_seeds',action = 'store_true',
                        help = 'Flag to restrict the network only to the list of seed proteins')
    parser.add_argument('-radius','--radius_of_subnetwork_around_seeds',dest='radius',default=0,action = 'store',type=int,
                        help = 'Network is built in a radius of connections around the seed proteins')
    parser.add_argument('-iloc','--input_localization_file',dest='input_of_localization',action = 'store',default='input_localization',
                        help = 'Input file with sub-cellular (or tissue-cell, or cancer) localization (default is input_localization)')
    parser.add_argument('-rloc','--restricted_to_known_localization',dest='restricted_to_localization',action = 'store_true',
                        help = 'Flag to use only nodes with known localization')
    parser.add_argument('-iarr','--Annotation_input_file',dest='input_array_annotation',action = 'store',default='input_annotation',
                        help = 'Annotation Input file of the array ESTs(default is input_annotation)')
    parser.add_argument('-oarr','--Annotation_output_file',dest='output_array_annotation',action = 'store',default='output_annotation',
                        help = 'Annotation Output file of the array ESTs(default is output_annotation)')
    parser.add_argument('-taxid','--TaxID',dest='taxid',action = 'store',default='9606',
                        help = 'Tax ID (i.e. human=9606 is default)')
    parser.add_argument('-stype','--seed_type',dest='stype',action = 'store',default='genesymbol',
                        help = 'Type of identifier for seeds (default is genesymbol)')
    parser.add_argument('-ttype','--translation_type',dest='ttype',action = 'store',default='accessionnumber',
                        help = '''Type of identifier for the output translation of codes (default is accessionnumber)
                        Using "proteinsequence" provides with the longest sequence of all codes''')
    parser.add_argument('-trans','--translation_of_nodes_file',dest='translation_file',action = 'store',default='translation_nodes.txt',
                        help = 'File with the translation of codes from BIANA to the selected type for all nodes')
    parser.add_argument('-node','--node_file',dest='node',action = 'store', default='biana_nodes',
                        help = 'Output file with nodes(default is biana_nodes)')
    parser.add_argument('-loc','--localization_file',dest='localization',action = 'store', default='biana_localization',
                        help = 'Output file with sub-cellular/tissue/cancer localization of nodes(default is biana_localization)')
    parser.add_argument('-edge','--edge_file',dest='edge',action = 'store', default='biana_edges',
                        help = 'Output file with edges(default is biana_edges)')
    parser.add_argument('-score','--Non-seed_score',dest='score',action = 'store',default='0.1',type=float,
                        help = 'Score of non-seed nodes (default is 0.1)')
    parser.add_argument('-nmeth','--number_of_methods',dest='nmeth',action = 'store',default='1',type=int,
                        help = 'Minimum number of methdos describing an interaction (default is 1)')
    parser.add_argument('-ndb','--number_of_databases',dest='ndb',action = 'store',default='1',type=int,
                        help = 'Minimum number of databases describing an interaction (default is 1)')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Flag to use verbose mode')
    parser.add_argument('-rAFF','--restricted_to_TAP',dest='restricted_to_TAP',action = 'store_true',
                        help = 'Flag to use interactions at least described by affinity methods (i.e. Tandem Affinity Purification)')
    parser.add_argument('-rY2H','--restricted_to_Y2H',dest='restricted_to_Y2H',action = 'store_true',
                        help = 'Flag to use interactions at least described by yeast two hybrid methods (Y2H)')
    parser.add_argument('-rUSER','--restricted_to_user',dest='restricted_to_user',action = 'store',default='restricted_methods',
                        help = 'File to use interactions described by the user selected methods')
    parser.add_argument('-eAFF','--except_TAP',dest='except_TAP',action = 'store_true',
                        help = 'Flag to use all interactions except those described by affinity methods (i.e. Tandem Affinity Purification)')
    parser.add_argument('-eY2H','--except_Y2H',dest='except_Y2H',action = 'store_true',
                        help = 'Flag to use all interactions except those described by yeast two hybrid methods (Y2H)')
    parser.add_argument('-eUSER','--except_user',dest='except_user',action = 'store',default='restricted_methods',
                        help = 'File to reject interactions described by the user selected methods')
    parser.add_argument('-format','--output_format',dest='format',action = 'store',default='guild',
                        help = 'Format files of nodes and edges:\tguild (default) or \tnetscore')
    options=parser.parse_args()

    return options

def fileExist (file):               #Checks if a file exists AND is a file

    return os.path.exists(file) and os.path.isfile(file)

def extract(options):

    affinity_dict={
    '492':'in vivo',
    '493':'in vitro',
    '0':'molecular interaction',
    '4':'affinity chromatography technology',
    '6':'anti bait coimmunoprecipitation',
    '7':'anti tag coimmunoprecipitation',
    '8':'array technology',
    '9':'bacterial display',
    '19':'coimmunoprecipitation',
    '28':'cosedimentation in solution',
    '29':'cosedimentation through density gradient',
    '30':'cross-linking',
    '34':'display technology',
    '47':'far western blotting',
    '48':'filamentous phage display',
    '49':'filter binding',
    '66':'lambda phage display',
    '71':'molecular sieving',
    '73':'mrna display',
    '81':'peptide array',
    '84':'phage display',
    '89':'protein array',
    '92':'protein in situ array',
    '95':'proteinchip(r) on a surface-enhanced laser desorption/ionization',
    '96':'pull down',
    '98':'ribosome display',
    '108':'t7 phage display',
    '115':'yeast display',
    '225':'chromatin immunoprecipitation array',
    '400':'affinity technology',
    '402':'chromatin immunoprecipitation assay',
    '405':'competition binding',
    '411':'enzyme linked immunosorbent assay',
    '412':'electrophoretic mobility supershift assay',
    '413':'electrophoretic mobility shift assay',
    '440':'saturation binding',
    '657':'systematic evolution of ligands by exponential enrichment',
    '676':'tandem affinity purification',
    '678':'antibody array',
    '695':'sandwich immunoassay',
    '729':'luminescence based mammalian interactome mapping',
    '813':'proximity enzyme linked immunosorbent assay',
    '858':'immunodepleted coimmunoprecipitation',
    '892':'solid phase assay',
    '899':'p3 filamentous phage display',
    '900':'p8 filamentous phage display',
    '921':'surface plasmon resonance array',
    '946':'ping',
    '947':'bead aggregation assay',
    '963':'interactome parallel affinity capture',
    '1017':'rna immunoprecipitation',
    '1028':'modified chromatin immunoprecipitation',
    '1029':'proteomics of isolated chromatin segments',
    '1031':'protein folding/unfolding',
    '1087':'monoclonal antibody blockade'
    }
    affinity=set(affinity_dict.keys())

    complementation_dict={
    '492':'in vivo',
    '493':'in vitro',
    '0':'molecular interaction',
    '10':	'beta galactosidase complementation',
    '11':	'beta lactamase complementation',
    '14':	'adenylate cyclase complementation',
    '18':	'two hybrid',
    '90':	'protein complementation assay',
    '97':	'reverse ras recruitment system',
    '111':	'dihydrofolate reductase reconstruction',
    '112':	'ubiquitin reconstruction',
    '228':	'cytoplasmic complementation assay',
    '229':	'green fluorescence protein complementation assay',
    '230':	'membrane bound complementation assay',
    '231':	'mammalian protein protein interaction trap',
    '232':	'transcriptional complementation assay',
    '369':	'lex-a dimerization assay',
    '370':	'tox-r dimerization assay',
    '397':	'two hybrid array',
    '398':	'two hybrid pooling approach',
    '399':	'two hybrid fragment pooling approach',
    '432':	'one hybrid',
    '437':	'protein tri hybrid',
    '438':	'rna tri hybrid',
    '588':	'3 hybrid method',
    '655':	'lambda repressor two hybrid',
    '726':	'reverse two hybrid',
    '727':	'lexa b52 complementation',
    '728':	'gal4 vp16 complementation',
    '809':	'bimolecular fluorescence complementation',
    '895':	'protein kinase A complementation',
    '916':	'lexa vp16 complementation',
    '1037':	'Split renilla luciferase complementation',
    }
    complementation=set(complementation_dict.keys())

    if not fileExist(options.restricted_to_user):
        print "No restriction on methods selected by the user"
        user_selection=False
    else:
        use_methods=[]
        input_method=open(options.restricted_to_user)
        for line in input_method:
            fields = line.strip().split("\t")
            use_methods.append(fields[0])
        input_method.close()
        user_selection=True
        #print "Input to use only Methods:",repr(use_methods)
    if not fileExist(options.except_user):
        print "No rejection of methods selected by the user"
        user_rejection=False
    else:
        no_methods=[]
        input_method=open(options.except_user)
        for line in input_method:
            fields = line.strip().split("\t")
            no_methods.append(fields[0])
        input_method.close()
        user_rejection=True
        #print "Input of rejected Methods:",repr(no_methods)

    session = create_new_session(sessionID="biana_session",dbname="BIANA_MARCH_2013",dbhost="ben-yehuda",
                                 dbuser="biana_user",dbpassword="biana_password",
                                 unification_protocol="uniprot_geneID_seqtax")

    if options.restricted_to_seeds or options.radius>0:
        if not fileExist(options.seed):
            print "File with seeds is missing or not found"
            sys.exit(10)
        else:
            level=options.radius
            seed_list=[]
            input_seed = open(options.seed,'r')
            for line in input_seed:
                fields = line.strip().split("\t")
                seed_list.append(fields[0])
            input_seed.close()
            proteome = session.create_new_user_entity_set(  identifier_description_list =seed_list,
                                                          attribute_restriction_list=[("taxid",options.taxid)],
                                                          id_type=options.stype,new_user_entity_set_id="proteome",
                                                          negative_attribute_restriction_list=[] )

    else:
        level=0
        proteome = session.create_new_user_entity_set( identifier_description_list = [("taxid",options.taxid)],
                                                      attribute_restriction_list=[], id_type="embedded",
                                                      new_user_entity_set_id="proteome",
                                                      negative_attribute_restriction_list=[] )

#Select interactions
    if options.restricted_to_TAP:
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               relation_attribute_restriction_list = [("Method_id",400)],
                                #relation_attribute_restriction_list  = [("psimi_name","affinity technology")],
                               include_relations_last_level = (not options.restricted_to_seeds) , use_self_relations = False)
    elif options.restricted_to_Y2H:
	#print "restricting to y2h"
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               relation_attribute_restriction_list = [("Method_id",18)],
                               #relation_attribute_restriction_list = [("psimi_name","y2h2")],
                               include_relations_last_level = (not options.restricted_to_seeds) , use_self_relations = False)
    else:
        session.create_network( user_entity_set_id = "proteome" , level = level, relation_type_list=["interaction"] ,
                               include_relations_last_level = (not options.restricted_to_seeds) , use_self_relations = False)



#Summary of interactions

    out_network = open(options.edge,'w')
    all_interactions = proteome.getRelations()
    print "Num interactions:", len(all_interactions)

    nodes=set()

    skip_interactions=0
    for (uE_id1, uE_id2) in all_interactions:

	#self.dbAccess.get_external_entities_dict( externalEntityIdsList = [external_entity_relation_id] )

        eErIDs_list = proteome.get_external_entity_relation_ids(uE_id1, uE_id2)

        method_names = set()
        method_ids = set()
        source_databases = set()
        use_method_ids=set()


        relationObj_dict = session.dbAccess.get_external_entities_dict(
                                externalEntityIdsList = eErIDs_list, attribute_list = [],
                                relation_attribute_list = ["method_id","psimi_name"], participant_attribute_list = [] )

        num_methods=0
        for current_eErID in eErIDs_list:
            relationObj = relationObj_dict[current_eErID]
            if options.verbose:
                print "Interaction: (",uE_id1,",",uE_id2,")"
                print relationObj

            #if relationObj.get_attribute(attribute_identifier="psimi_name") is not None:
            #    print "\t".join([ x.value for x in relationObj.get_attribute(attribute_identifier="psimi_name") ])
            #if relationObj.get_attribute(attribute_identifier="method_id") is not None:
            #print "\t".join([ x.value for x in relationObj.get_attribute(attribute_identifier="method_id") ])
            #print relationObj.get_attributes_dict()
            #print [ x.value for x in relationObj.get_attributes_dict()["psimi_name"] ]
            #print [ x.value for x in relationObj.get_attributes_dict()["method_id"] ]

            if "psimi_name" in relationObj.get_attributes_dict():
                method_names.update([ str(x.value) for x in relationObj.get_attributes_dict()["psimi_name"] ])
            if "method_id" in relationObj.get_attributes_dict():
                method_ids.update([ x.value for x in relationObj.get_attributes_dict()["method_id"]])
            source_databases.add(str(session.dbAccess.get_external_database(
                                    database_id = relationObj.get_source_database()) ))
            if options.except_TAP:
                #print "check TAP"
                for m in method_ids:
                    if m not in affinity:
                        use_method_ids.add(m)
                        #print "Add", m
            elif options.except_Y2H:
                #print "check Y2H"
                for m in method_ids:
                    if m not in complementation:
                        use_method_ids.add(m)
                        #print "Add", m
            elif user_rejection:
                for m in method_ids:
                    if m not in no_methods:
                        use_method_ids.add(m)
            elif user_selection:
                for m in method_ids:
                    #print "Check",repr(use_methods)
                    if m in set(use_methods):
                        use_method_ids.add(m)
                    elif options.verbose:
                        print "Not among selected methods ",m
            else:
                use_method_ids.update(method_ids)

        info_sources=";".join([str(x) for x in source_databases])
        info_methods=";".join([str(x) for x in method_names])
        info_methods_ids=";".join([str(x) for x in use_method_ids])
        num_databases=len(source_databases)
        num_methods=len(use_method_ids)

        if options.verbose:
            print "Methods",num_methods,info_methods,"\tSelected:",info_methods_ids
            print "Databases",num_databases,info_sources

        if num_methods >= options.nmeth:
            use=True
        else:
            use=False

        if use and num_databases >= options.ndb:
            use=True
        else:
            use=False

        if not use:
            skip_interactions=skip_interactions+1
        #print method_names, method_ids, source_databases
        if use:
            #print uE_id1, uE_id/2
            nodes.add(uE_id1)
            nodes.add(uE_id2)
            #print "Attribute ",(uE_id1,uE_id2).get_attribute(
            if options.verbose:
                #out_network.write("\t%s\t%s\t%f\t\t# Methods:%s # Databases:%s \n" %(uE_id1,uE_id2,1,info_methods,info_sources))
                if options.format == 'netscore' :
                    out_network.write("\t{0}\t{1}\t{2:f}\t\t# Methods:{3} # Databases:{4} \n".
                                  format(uE_id1,uE_id2,1.,info_methods,info_sources))
                else:
                    out_network.write("{0} {1:f} {2}\t\t# Methods:{3} # Databases:{4} \n".
                                  format(uE_id1,1.,uE_id2,info_methods,info_sources))
            else:
                if options.format == 'netscore' :
                #out_network.write("\t%s\t%s\t%f\n" %(uE_id1,uE_id2,1))
                    out_network.write("\t{0}\t{1}\t{2:f}\n".format(uE_id1,uE_id2,1.))
                else:
                    out_network.write("{0} {1:f} {2}\n".format(uE_id1,1.,uE_id2))

    print "Num neglected interactions:", skip_interactions
    out_network.close()

    out_annotation= file(options.output_array_annotation,'w')
    if not fileExist(options.input_array_annotation):
         for protein in nodes:
             out_annotation.write("{0}\t{1}\n".format(protein,protein))
    else:
         fd=open(options.input_array_annotation,'r')
         accnum2microarray = {}
         genesymbol2microarray = {}
         translations = {}
         for line in fd:
             fields = line.strip().split("\t")
             accnum2microarray.setdefault(fields[2].lower(), []).append(fields[0])
             genesymbol2microarray.setdefault(fields[1].lower(), []).append(fields[0])
         fd.close()
         for protein in nodes:
             for microarray_codes in translations.setdefault(protein, get_translation(protein,session,accnum2microarray,genesymbol2microarray)):
                 out_annotation.write("{0}\t{1}\n".format(protein, microarray_codes))
    out_annotation.close()

    if not fileExist(options.seed):
        out_proteins = open(options.node,'w')
        for protein in nodes:
            if options.format=='netscore':
                #out_proteins.write("%s\t%10.5f\t%10.5f\t%10.5f\n" %(protein,1,1,options.score))
                out_proteins.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(protein,1.,1.,options.score))
            else:
                out_proteins.write("{0} {1:10.5f}\n".format(protein,score))
        out_proteins.close()
        out_translation = open(options.translation_file,'w')
        for protein in nodes:
            uE = session.get_user_entity(protein)
            if option.ttype == "proteinsequence":
                maxlen=0;
                for current_id in uE.get_attribute(attribute_identifier=options.ttype):
                    if maxlen < len(current_id.value.get_sequence().upper()):
                        maxlen=len(current_id.value.get_sequence().upper())
                translation=",".join([str(current_id.value.get_sequence().upper()) for current_id in uE.get_attribute(attribute_identifier=options.ttype) if len(str(current_id.value.get_sequence().upper())) == maxlen ] )
                #print "Translation",protein,translation
                #print("{0}\t'{1}'\n".format(protein,translation))
            else:
                translation=",".join([str(current_id) for current_id in uE.get_attribute(attribute_identifier=options.ttype)])
            #out_translation.write("%s\t%s\n" % (protein,translation))
            out_translation.write("{0}\t'{1}'\n".format(protein,translation))
        out_translation.close()
        out_localization = open(options.localization,'w')
        if not fileExist(options.input_of_localization):
            for protein in nodes:
                #out_localization.write("%25s\t%10d\t%10d\t%10d\n" %(protein,1,0,0))
                out_localization.write("{0:>25s}\t{1:10d}\t{2:10d}\t{3:10d}\n".format(str(protein),1,0,0))
        else:
            localization=getLocalization(options.input_of_localization)
            places=set()
            for g,u in localization.items():
                for l in u:
                    places.add(l)
            maxloc=len(places)
            for protein in nodes:
                uE = session.get_user_entity(protein)
                value=set()
                for current_id in uE.get_attribute(attribute_identifier="ensembl"):
                    value.update(set(localization.get(current_id)))
                listed=list(value)
                listed.sort()
                if len(value)>0:
                    #subcel="\t".join([str(x) for x in listed])
                    subcel="\t".join(["{0:10d}".format(x) for x in listed])
                elif not options.restricted_to_localization:
                    #subcel="\t".join([str(x) for x in range(1,maxloc+1)])
                    subcel="\t".join(["{0:10d}".format(x) for x in range(1,maxloc+1)])
                #print "Protein ",protein," Subcel ",subcel
                out_localization.write("{0:>25s}\t{1}\n".format(str(protein),subcel))
        out_localization.close()
    else:
        seeds=set()
        input_seed = open(options.seed,'r')
        for line in input_seed:
            fields = line.strip().split("\t")
            seeds.add(fields[0].lower())
        input_seed.close()
        out_proteins = open(options.node,'w')
        translate={}
        for protein in nodes:
            score=options.score
            uE = session.get_user_entity(protein)
            for current_id in uE.get_attribute(attribute_identifier=options.stype):
                if current_id.value.lower() in seeds:
                    translate.setdefault(current_id.value.lower(),[])
                    translate[current_id.value.lower()].append(protein)
                    score=1.0
            if options.format=='netscore':
                #out_proteins.write("%s\t%10.5f\t%10.5f\t%10.5f\n" %(protein,1,1,score)) #obsolete format
                #out_proteins.write("{0:s}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(str(protein),1.,1.,score)) #overepresentation of format
                out_proteins.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(protein,1.,1.,score))
            else:
                out_proteins.write("{0} {1:10.5f}\n".format(protein,score))

        out_proteins.close()


# Get the IDS of single nodes that were not previously found in the network

        single=set()
        for uE_id in proteome.get_unconnected_nodes():
            single.add(uE_id)
        for protein in single:
            uE = session.get_user_entity(protein)
            for current_id in uE.get_attribute(attribute_identifier=options.stype):
                if current_id.value.lower() in seeds:
                    translate.setdefault(current_id.value.lower(),[])
                    translate[current_id.value.lower()].append(protein)

# Get all IDS of SEEDS, defined as "proteome", and check missing codes to be
# added for translation

        allseed=set()
        for uE_id in proteome.get_user_entity_ids():
            allseed.add(uE_id)
        for protein in allseed:
            if protein not in single and protein not in nodes:
                uE = session.get_user_entity(protein)
                for current_id in uE.get_attribute(attribute_identifier=options.stype):
                    if current_id.value.lower() in seeds:
                        translate.setdefault(current_id.value.lower(),[])
                        translate[current_id.value.lower()].append(protein)



        out_translation = open("translation_seeds_to_BIANA_codes.txt",'w')
        for s in seeds:
            if s == '': continue
            if s in translate:
                codes=set(translate[s])
                translation="','".join([str(x) for x in codes])
                #out_translation.write("%s\t'%s'\n" % (s.upper(),translation))
                out_translation.write("{0}\t'{1}'\n".format(s.upper(),translation))
            else:
                out_translation.write("{0}\t'Unknown'\n".format(s.upper()))
        out_translation.close()
        out_translation = open(options.translation_file,'w')
        for protein in nodes:
            uE = session.get_user_entity(protein)
            translate=set()
            if options.ttype == "proteinsequence":
                maxlen=0;
                for current_id in uE.get_attribute(attribute_identifier=options.ttype):
                    if maxlen < len(current_id.value.get_sequence().upper()):
                        maxlen=len(current_id.value.get_sequence().upper())
                translation=",".join([str(current_id.value.get_sequence().upper()) for current_id in uE.get_attribute(attribute_identifier=options.ttype) if len(str(current_id.value.get_sequence().upper())) == maxlen ] )
                #print "Translation",protein,translation
                #print("{0}\t'{1}'\n".format(protein,translation))
            else:
                for current_id in uE.get_attribute(attribute_identifier=options.ttype):
                    translate.add(current_id.value.upper())
                translation="','".join(["{0}".format(x) for x in translate])
            out_translation.write("{0}\t'{1}'\n".format(protein,translation))
        out_translation.close()
        out_localization = open(options.localization,'w')
        if not fileExist(options.input_of_localization):
            for protein in nodes:
                out_localization.write("{0:>25s}\t{1:10d}\t{2:10d}\t{3:10d}\n".format(str(protein),1,0,0))
        else:
            localization=getLocalization(options.input_of_localization)
            places=set()
            for g,u in localization.items():
                for l in u:
                    places.add(l)
            maxloc=len(places)
            for protein in nodes:
                uE = session.get_user_entity(protein)
                value=set()
                #print "UnitEntity ", repr(uE)
                #print "UnitEntity Attribute ", repr(uE.get_attribute(attribute_identifier="ensembl"))
                for current_id in uE.get_attribute(attribute_identifier="ensembl"):
                    #print "Ensembl ",current_id.value.lower()
                    if current_id.value.lower() in localization:
                        value.update(set(localization.get(current_id.value.lower())))
                listed=list(value)
                listed.sort()
                if len(value)>0:
                    #subcel="\t".join([str(x) for x in listed])
                    subcel="\t".join(["{0:10d}".format(x) for x in listed])
                elif not options.restricted_to_localization:
                    #subcel="\t".join([str(x) for x in range(1,maxloc+1)])
                    subcel="\t".join(["{0:10d}".format(x) for x in range(1,maxloc+1)])
                #print "Protein ",protein," Subcel ",subcel
                out_localization.write("{0:>25s}\t{1}\n".format(str(protein),subcel))
        out_localization.close()

    print "\nDone\n"
    sys.exit(10)

def getLocalization(file_localization):
    input_localization = open(file_localization,'r')
    title=True
    region=set()
    localname={}
    localization={}
    for line in input_localization:
        #fields = line.strip().split(",")
        #field = re.compile(r'"(\.+)"')
        field = line.strip().split(",")
        fraction=[]
        if  title:
            title=False
            subcellular=False
            tissue=False
            cancer=False
            for word in field:
                fraction.append(word.strip('"'))
            if fraction[1]=="Main location":
                subcellular=True
            elif fraction[1]=="Tissue":
                tissue=True
            elif fraction[1]=="Tumor":
                cancer=True
            else:
                print "The file with Localization is unknown\n"
                sys.exit(10)
        else:
            if subcellular:
                for word in field:
                    fraction.append(word.strip('"'))
                gene=fraction[0].lower()
                local=fraction[1].lower()+";"+fraction[2].lower()
                for ubiq in local.split(";"):
                    if ubiq == '': continue
                    region.add(ubiq)
                    localname.setdefault(gene, [])
                    localname[gene].append(ubiq)
                    if ubiq == "nucleus but not nucleoli":
                        localname[gene].append("nucleus")
                    if ubiq == "cytoskeleton (microtubules)":
                        localname[gene].append("cytoskeleton")
                    if ubiq == "cytoskeleton (intermediate filaments)":
                        localname[gene].append("cytoskeleton")
                    if ubiq == "cytoskeleton (actin filaments)":
                        localname[gene].append("cytoskeleton")
                    if ubiq == "cytoskeleton (cytokinetic bridge)":
                        localname[gene].append("cytoskeleton")
            if tissue:
                for word in field:
                    fraction.append(word.strip('"'))
                gene=fraction[0].lower()
                local=fraction[2].lower()+" in "+fraction[1].lower()
                region.add(local)
                localname.setdefault(gene, [])
                localname[gene].append(local)
            if cancer:
                for word in field:
                    fraction.append(word.strip('"'))
                gene=fraction[0].lower()
                local=fraction[1].lower()
                region.add(local)
                localname.setdefault(gene, [])
                localname[gene].append(local)

    output_localization=open("localization_definition.txt",'w')
    dictionary_region={}
    n=1
    for u in region:
        if u == '': continue
        #print "N ",n,"Localization ",repr(u)
        dictionary_region.setdefault(u,n)
        output_localization.write("%d\t%s\n" %(n,u))
        n=n+1
    output_localization.close()
    input_localization.close()
    for gene,ubiq in localname.items():
        localization.setdefault(gene,[])
        for sub in ubiq:
            if sub == '':continue
            #print "Gene ",gene, "Localization ", sub, "Value ",dictionary_region[sub]
            localization[gene].append(dictionary_region[sub])
    return localization

def get_translation(uE_id,session,accnum2microarray,genesymbol2microarray):

	uE = session.get_user_entity(uE_id)

	microarray_codes = set()

	for current_acc in uE.get_attribute(attribute_identifier="accessionnumber"):
		if current_acc.value.lower() in accnum2microarray:
			microarray_codes.update(accnum2microarray[current_acc.value.lower()])
	# If not found, the use genesymbol
	if len(microarray_codes)==0:
		for current_gs in uE.get_attribute(attribute_identifier="genesymbol"):
			if current_gs.value.lower() in genesymbol2microarray:
				microarray_codes.update(genesymbol2microarray[current_gs.value.lower()])

	return microarray_codes




if  __name__ == "__main__":
    main()



