from libs import options
from libs import utils
from methods import method

class Test(method.Method):
	def __init__(self,gn_network,tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		for line in utils.readTable("{}data/{}/annotation.gtf".format(options.Options().wd,options.Options().annotation),header=False):
			if "uc010ssa" in line[8]:
				import pdb
				pdb.set_trace()
			else:
				continue

			seqname 	= line[0]
			source 		= line[1]
			feature 	= line[2]
			start 		= int(line[3])
			end 		= int(line[4])
			score 		= line[5]
			strand 		= line[6]
			frame 		= line[7]
			attribute 	= line[8].split(";")

			tx = [ x.split(" ")[2].strip("\"") for x in attribute if "transcript_id" in x ][0]

			if tx not in self._transcript_network.nodes():
				continue

			if feature=="exon":
				exon = [ start, end ]
				self._transcript_network.update_node(tx,"exonStructure",exon)

			elif feature=="CDS":
				cds.setdefault(tx,[])
				cds[tx].append([start,end])

			self._transcript_network.update_node(tx, "strand", strand)
			self._transcript_network.update_node(tx, "chr", seqname)