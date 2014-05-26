thisGene="ARL1"

genes = {}

with open("Data/TCGA/UnifiedFasta_iLoops_devel.fa", "r") as KK:
	for line in KK:
		if ">" in line:
			gene = line.strip().split("#")[1].split("|")
			if len(gene) == 1:
				genes[gene[0]] = gene[0]
			elif len(gene) == 2:
				genes[gene[0]] = gene[1]

with open("/home/hector/Desktop/Baldo_" + thisGene + ".tsv", "r") as KK2, open("/home/hector/Desktop/Baldo_" + thisGene + "_2.tsv", "w") as KK3:
	for line in KK2:
		element = line.split("\t")
		KK3.write(genes[element[0]] + "\t" + genes[element[1]] + "\t" + element[2])

