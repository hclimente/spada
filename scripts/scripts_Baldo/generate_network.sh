
fgrep -v locus nodes_scored.txt> nodes_scored_noLocus.txt
fgrep -v locus edges_scored.txt> edges_scored_noLocus.txt

/soft/devel/python-2.7/bin/python ~/PROGRAMS/Interactome/scripts/generate_netscore_files.py -iseed nodes_scored_noLocus.txt -radius 5 -taxid 9606 -rUSER restricted_methods.txt -trans translate_network.txt -node network.nodes -edge network.edges -stype geneid -ttype geneid -score 0.05 -v >& network.log 

/soft/devel/python-2.7/bin/python ~/PROGRAMS/Interactome/scripts/generate_netscore_files.py -iseed nodes_scored_noLocus.txt -radius 5 -taxid 9606 -eAFF -trans translate_long_network.txt -node long_network.nodes -edge long_network.edges -stype geneid -ttype geneid -score 0.01 -v >& long_network.log 

# Transform scores of long_network and network
# Nodes: 0.01 for long_network and 0.02 for network
# edges: 0.2 for long-network and 1.0 for network


/soft/devel/python-2.7/bin/python TranslateNetwork.py -i network.edges -n network.nodes -iformat guild -oformat guild -oe network_geneid_guild.edges -on network_geneid_guild.nodes -trans translate_network.txt

/soft/devel/python-2.7/bin/python TranslateNetwork.py -i long_network.edges -n long_network.nodes -iformat guild -oformat guild -oe long_network_geneid_guild.edges -on long_network_geneid_guild.nodes -trans translate_long_network.txt

#/soft/devel/python-2.7/bin/python TranslateNetwork.py -i network.edges -n network.nodes -iformat guild -oformat netscore -oe network_geneid_netscore.edges -on network_geneid_netscore.nodes -trans translate_network.txt

#/soft/devel/python-2.7/bin/python TranslateNetwork.py -i long_network.edges -n long_network.nodes -iformat guild  -oformat netscore -oe long_network_geneid_netscore.edges -on long_network_geneid_netscore.nodes -trans translate_long_network.txt

#Add and prune the PPI network

#On BIANA codes (unnecessary we will use GeneID on GUILD

#/soft/devel/python-2.7/bin/python PruneAddNetwork.py -e network.edges -n network.nodes -ae long_network.edges -an long_network.nodes -oe network_ppi.edges -on network_ppi.nodes
#/soft/devel/python-2.7/bin/python PruneAddNetwork.py -e network.edges -n network.nodes -ae long_network.edges -an long_network.nodes -oe network_ppi_netscore.edges -on network_ppi_netscore.nodes -oformat netscore

#On geneid codes

/soft/devel/python-2.7/bin/python PruneAddNetwork.py -e network_geneid_guild.edges -n network_geneid_guild.nodes  -ae long_network_geneid_guild.edges -an long_network_geneid_guild.nodes -oe network_geneid_guild_ppi.edges -on network_geneid_guild_ppi.nodes

/soft/devel/python-2.7/bin/python PruneAddNetwork.py -e network_geneid_guild.edges -n network_geneid_guild.nodes  -ae long_network_geneid_guild.edges -an long_network_geneid_guild.nodes -oe network_geneid_netscore_ppi.edges -on network_geneid_netscore_ppi.nodes -oformat netscore


#Modify edge score 1.2 back to 1.0
#Modify Node scores from 2 back to 1 and 0.03 (and 0.02) back to 0.01 (GUILD) or 0.0 (Netscore)
#save files as network_geneid_guild_ppi_bin.nodes  and network_geneid_netscore_ppi_bin.nodes
#check with sort -g -r -k 4 (netscore) and -k 2 (guild)

#Modify ALL Nodes to 0.01 (GUILD) and 0.0 (netscore) 

/soft/devel/python-2.7/bin/python  NullifyScore.py network_geneid_guild_ppi.nodes network_geneid_guild_ppi_null.nodes
/soft/devel/python-2.7/bin/python  NullifyScore.py network_geneid_netscore_ppi.nodes network_geneid_netscore_ppi_null.nodes netscore

#Add and prune the PPI network (without other methods than y2h, otherwise use network_geneid_guild_ppi.edges) and the SPLICING specific network

/soft/devel/python-2.7/bin/python PruneAddNetwork.py -e network_geneid_guild.edges -n network_geneid_guild_ppi_null.nodes -ae  edges_scored_noLocus.txt -an nodes_scored_noLocus.txt -oe ppi_splicing.edges -on ppi_splicing.nodes -v

/soft/devel/python-2.7/bin/python PruneAddNetwork.py -e edges_scored_noLocus.txt -n nodes_scored_noLocus.txt -oe edges_scored_noLocus_netscore.txt -on nodes_scored_noLocus_netscore.txt -oformat netscore 

/soft/devel/python-2.7/bin/python PruneAddNetwork.py -e network_geneid_netscore.edges -n network_geneid_netscore_ppi_null.nodes -ae  edges_scored_noLocus_netscore.txt -an nodes_scored_noLocus_netscore.txt -oe ppi_splicing_netscore.edges -on ppi_splicing_netscore.nodes -v -oformat netscore -iformat netscore

#Run priporitization

mkdir GUILD NETSCORE
cp ppi_splicing.* GUILD
cp ppi_splicing_netscore.* NETSCORE


