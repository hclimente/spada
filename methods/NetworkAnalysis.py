#!/soft/devel/python-2.7/bin/python

import pandas as pd
import numpy as np
import sys
from libsmartas import *
import networkx as nx


class Table:
	def __init__(self, path="", df=None):
		if path: self._table = pd.DataFrame.from_csv(path, sep="\t")
		elif df: self._table = df
	## GETTERS ##
	def getTable(self): return self._table
	## SETTERS ##
	def setTable(self, table): self._table = table
	def toCsv(self, path): self.getTable().to_csv(path, sep="\t", index=False)

class EdgeTable(Table):
	def __init__(self, path="", df=None):
		Table.__init__(self, path, df)

class NodeTable(Table):
	def __init__(self, path="", df=None):
		Table.__init__(self, path, df)

class Network:
	def __init__(self, nodesTable, edgesTable):
		self._net = nx.Graph()

		# network_input_file		network.edges
		# nodes_input_file 			network.nodes
		# translation_of_nodes_file	translate_network.txt
		# input_format				guild
		# output_format				guild
		# output_file_edges			network_geneid_guild.edges
		# output_file_nodes			network_geneid_guild.nodes

		self.translate()

	def translate(self):
		#check that both nodes and edges exists, as well as translation file.

		new = {}
		with open("Data/Databases/translate_network.txt", "r") as TRANSLATION:
			for line in TRANSLATION:
				elements = line.strip().split("\t")
				node=word[0]
				new.setdefault(node,set())
	    
				for translation in [ x for x in word[1].split("'") if x and x != "," ]:
         			name="_".join([str(x) for x in translation.split()])
					new[node].add(name)
	

	  #   fd=open(options.input_edge,"r")
	  #   out_network=options.output_edge
	  #   if not options.output_edge == sys.stdout: out_network=open(options.output_edge,"w")
	
	  #   nodes=set()
	  #   for line in fd:
	  #     word=line.split()
	  #     p=word[0]
	  #     score=float(word[1])
	  #     q=word[2]
	  #     info="".join([str(word[i]) for i in range(3,len(word))])
	  #     for a in new[p]:
	  #      nodes.add(a)
	  #      for b in new[q]:
	  #         nodes.add(b)
	  #         out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
	  #   fd.close()

	  #   if not options.output_edge == sys.stdout: out_network.close()
	
	  #   out_network=options.output_node
	  #   if not options.output_node == sys.stdout: out_network=open(options.output_node,"w")
	
	  #   if fileExist(options.input_node):
	  #    fd=open(options.input_node,"r")
	  #    for line in fd:
	  #     word=line.split()
	  #     p=word[0]
	  #     score=float(word[1])
	  #     for a in new[p]:
			# out_network.write("{0} {1:10.5f}\n".format(a,score))
	  #    fd.close()
	  #   else:
	  #    for a in nodes:
			# out_network.write("{0} {1:10.5f}\n".format(a,0.01))
	
	  #   out_network.close()

if __name__ == '__main__':
	#path = sys.argv[1]
	path = "/home/hector/SmartAS/Results/TCGA/luad_mE-1.0/"
	samples = 57
	candidates = { "ARL1" : {"nTx": "uc001tib.2", "tTx": "uc001tic.2"} }
	samples = round(samples * 0.1)
	#Create nodes table
	Candidates = Table(path=path + "candidateList.top.tsv")
	Candidates.setTable( Candidates.getTable()[ ["Gene", "Replicates"] ] )
	#Remove less significant switches
	Candidates.setTable( Candidates.getTable()[ Candidates.getTable().Replicates >= samples ] )
	Candidates.setTable( Candidates.getTable().groupby("Gene").sum() )
	Candidates.getTable().Replicates = 0.2 + 0.3 * (Candidates.getTable().Replicates - samples)/(57 - samples)
	#Create edges table
	edges = []
	for gene in candidates.keys():
		InteraX = Table(path=path + "iLoops/InteraXChanges_" + gene + "_" + candidates[gene]["nTx"] + "_" + candidates[gene]["tTx"] + "_full1.tsv")
		#Filter out "locus" genes, as they are not understandable by GUILD
		locus_filter = [ "locus" not in x for x in [ x for x in InteraX.getTable().Partner_gene ] ]
		InteraX.setTable( InteraX.getTable()[ locus_filter ] )
		#Select genes for which there is only a prediction for one isoform
		dRC_filter = np.abs(InteraX.getTable().dRC) <= 50
		for aGene in InteraX.getTable().Partner_gene.unique():
			genenameMask = InteraX.getTable().Partner_gene == aGene
			geneInfo 	 = InteraX.getTable()[ genenameMask ]
			geneScore 	 = 0.0
			if any(geneInfo.Annotation == "Driver"):
				geneScore += 0.5
			if ( genenameMask & dRC_filter ).any():
				geneScore += float( max(np.abs(geneInfo.dRC)) ) / 50
			else:
				geneScore += 0.5
			edges.append( pd.Series([gene, geneScore, aGene], index=["Gene", "Score", "Partner"]) )
	SmartAS_edges = EdgeTable(df=pd.DataFrame(edges))
	SmartAS_nodes = NodeTable(df=Candidates.getTable())
	SmartAS_nodes.toCsv("/home/hector/Desktop/kk_nodes.tsv")
	SmartAS_edges.toCsv("/home/hector/Desktop/kk_edges.tsv")

	# #Falta vaciarlos de LOCUS
	# cmd("fgrep -v locus nodes_scored.txt> nodes_scored_noLocus.txt")
	# cmd("fgrep -v locus edges_scored.txt> edges_scored_noLocus.txt")
	# cmd("/soft/devel/python-2.7/bin/python ~/PROGRAMS/Interactome/scripts/generate_netscore_files.py -iseed nodes_scored_noLocus.txt -radius 5 -taxid 9606 -rUSER restricted_methods.txt -trans translate_network.txt -node network.nodes -edge network.edges -stype geneid -ttype geneid -score 0.05 -v >& network.log")
	# cmd("/soft/devel/python-2.7/bin/python ~/PROGRAMS/Interactome/scripts/generate_netscore_files.py -iseed nodes_scored_noLocus.txt -radius 5 -taxid 9606 -eAFF -trans translate_long_network.txt -node long_network.nodes -edge long_network.edges -stype geneid -ttype geneid -score 0.01 -v >& long_network.log")

	# SmartAS_nodes = NodeTable("long_network.nodes")
	# SmartAS_edges = EdgeTable("long_network.edges")

#gene <- "IGF2BP3"
#tx1 <- "uc003swf.2"
#tx2 <- "uc003swg.2"

kk <- vector("list", length(unique(InteraX$Partner_gene)))

mask1 <- abs(InteraX$dRC) != 9999

for (aGene in unique(InteraX$Partner_gene)){
  mask2 <- InteraX$Partner_gene == aGene
  score <- 0
  if (any(mask1 & mask2)){
    score <- score + max(abs(InteraX$dRC[mask1 & mask2]))/50  
  }
  else{
    score <- score + 0.5
  }
  
  if (any(InteraX$Annotation[mask2] == "Driver")){
    score <- score + 0.5
  }
  kk[[aGene]] <- c(gene, aGene, score)
}

perrote <- do.call("rbind", kk)
write.table(perrote, file=paste0("~/Desktop/Baldo_", gene, ".tsv"), sep="\t", row.names=F, col.names=F, quote=F)