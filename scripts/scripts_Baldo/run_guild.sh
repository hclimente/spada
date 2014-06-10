# prepare the data for GUILD

#cut -f 1,2 ppi_splicing.edges > ppi_splicing_guild.edges

# Create 100 random networks for NetZcore

#/soft/devel/python-2.7/bin/python ~/PROJECTS/GUILD/guild/src/create_random_networks_for_netzcore.py ppi_splicing_guild.edges 100

# Run NetScore/NetZcore/NetShort

 ~/PROJECTS/GUILD/guild/guild -s s -n ppi_splicing.nodes -e ppi_splicing_guild.edges -o ppi_splicing_guild_netScore.out -r 3 -i 2

 ~/PROJECTS/GUILD/guild/guild -s z -n ppi_splicing.nodes -e ppi_splicing_guild.edges -o ppi_splicing_guild_netZcore.out -d ppi_splicing_guild.edges. -x 100

 ~/PROJECTS/GUILD/guild/guild -s d -n ppi_splicing.nodes -e ppi_splicing_guild.edges -o ppi_splicing_guild_netShort.out 

# Run netCombo

/soft/devel/python-2.7/bin/python ~/PROGRAMS/Interactome/GUILD/guild/src/combine_scores.py ppi_splicing_guild_netScore.out ppi_splicing_guild_netZcore.out ppi_splicing_guild_netShort.out  ppi_splicing_guild_netCombo.out

# Run Functional Flow

 ~/PROJECTS/GUILD/guild/guild -s f -n ppi_splicing.nodes -e ppi_splicing_guild.edges -o ppi_splicing_guild_fFlow.out -i 5 -t 1.0

# run NetRank

 ~/PROJECTS/GUILD/guild/guild -s r -n ppi_splicing.nodes -e ppi_splicing_guild.edges -o ppi_splicing_guild_netRank.out

# run Full Combo
# first modif fFlow, change "inf" for a large value (i.e. its maximum)
# in this case the maximum is 14.32, so we take 15. name the file as ppi_splicing_guild_fFlow_mod.out

/soft/devel/python-2.7/bin/python ~/PROGRAMS/Interactome/GUILD/guild/src/combine_scores.py ppi_splicing_guild_netScore.out ppi_splicing_guild_netZcore.out ppi_splicing_guild_netShort.out ppi_splicing_guild_fFlow_mod.out ppi_splicing_guild_netRank.out ppi_splicing_guild_FullCombo.out

# reorder the files to get top 1% (around 114 first nodes) or top5% (first 570 nodes)

foreach net ("netScore" "netZcore" "netShort" "netCombo" "fFlow" "netRank" "FullCombo")
sort -g -r -k 2 ppi_splicing_guild_${net}.out  |head -114 > ${net}_top1.txt 
sort -g -r -k 2 ppi_splicing_guild_${net}.out  |head -570 > ${net}_top5.txt 
end

#sort -g -r -k 2 ppi_splicing_guild_netScore.out  |head -114 > netScore_top1.txt 
#sort -g -r -k 2 ppi_splicing_guild_netZcore.out  |head -114 > netZcore_top1.txt 
#sort -g -r -k 2 ppi_splicing_guild_netShort.out  |head -114 > netShort_top1.txt 
#sort -g -r -k 2 ppi_splicing_guild_netCombo.out  |head -114 > netCombo_top1.txt 
#sort -g -r -k 2 ppi_splicing_guild_fFlow.out     |head -114 > netfFlow_top1.txt 
#sort -g -r -k 2 ppi_splicing_guild_netRank.out   |head -114 > netRank_top1.txt 
#sort -g -r -k 2 ppi_splicing_guild_FullCombo.out |head -114 > FullCombo_top1.txt 
#
#sort -g -r -k 2 ppi_splicing_guild_netScore.out  |head -570 > netScore_top5.txt 
#sort -g -r -k 2 ppi_splicing_guild_netZcore.out  |head -570 > netZcore_top5.txt 
#sort -g -r -k 2 ppi_splicing_guild_netShort.out  |head -570 > netShort_top5.txt 
#sort -g -r -k 2 ppi_splicing_guild_netCombo.out  |head -570 > netCombo_top5.txt 
#sort -g -r -k 2 ppi_splicing_guild_fFlow.out     |head -570 > netfFlow_top5.txt 
#sort -g -r -k 2 ppi_splicing_guild_netRank.out   |head -570 > netRank_top5.txt 
#sort -g -r -k 2 ppi_splicing_guild_FullCombo.out |head -570 > FullCombo_top5.txt 


# Select the TOP subnetworks

foreach net ("netScore" "netZcore" "netShort" "netCombo" "fFlow" "netRank" "FullCombo")
foreach top ("top1" "top5")
echo "Test Network $net with top $top percent\n"
echo "Check splicing-special edges edges_scored_noLocus\n"
/soft/devel/python-2.7/bin/python  SelectSubNetwork.py -e ../edges_scored_noLocus.txt    -oe edges_scored_${net}_$top.edge -n ${net}_$top.txt -on edges_scored_${net}_$top.node
/soft/devel/python-2.7/bin/python  TransformNetwork.py -i  edges_scored_${net}_$top.edge -oe edges_scored_${net}_$top.sif  -oformat sif
echo "Check PPI edges \n"
/soft/devel/python-2.7/bin/python  SelectSubNetwork.py -e ../network_geneid_guild.edges  -oe ppi_network_${net}_$top.edge  -n ${net}_$top.txt -on ppi_network_${net}_$top.node
/soft/devel/python-2.7/bin/python  TransformNetwork.py -i  ppi_network_${net}_$top.edge  -oe ppi_network_${net}_$top.sif   -oformat sif
echo "Check all edges (PPI plus splicing-special)\n"
/soft/devel/python-2.7/bin/python  SelectSubNetwork.py -e ../ppi_splicing.edges          -oe all_${net}_$top.edge          -n ${net}_$top.txt -on all_${net}_$top.node
/soft/devel/python-2.7/bin/python  TransformNetwork.py -i  all_${net}_$top.edge          -oe all_${net}_$top.sif           -oformat sif
end
end



