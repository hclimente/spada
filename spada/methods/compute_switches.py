from spada.biological_entities.switch import IsoformSwitch
from spada.biological_entities.gene_expression import GeneExpression
from spada.io import io
from spada.methods import method

class ComputeSwitches(method.Method):

	def __init__(self, gn_network, tx_network, ctrlFile, caseFile):

		method.Method.__init__(self, __name__,gn_network,tx_network)

		self._CTRL = open(ctrlFile, "r")
		self._CASE = open(caseFile, "r")

		self._genes.flushSwitches()

	def run(self, minExpression):

		self.findSwitches(minExpression)
		self._genes.calculateCompatibilityTable()
		io.printSwitches(self._genes, self._txs)
		self._genes.saveNetwork("genes.pkl")

	def findSwitches(self, minExpression):

		gene2tx = self.getGene2Tx()
		idsCtrl = self.readSamples(self._CTRL)
		idsCase = self.readSamples(self._CASE)

		# gene -> xpr
		expression = { }

		for (tx,ctrl),(tx2,case) in zip(io.parseExpression(self._CTRL), io.parseExpression(self._CASE)):

			if tx != tx2:
				raise SpadaError("Case and control expresion files mismatch: {} vs. {}.".format(tx, tx2))

			try:
				gene = self._txs._net.node[tx]["gene_id"]
			except KeyError:
				self.logger.debug('Transcript {} not in network'.format(tx))
				continue

			expression.setdefault(gene, GeneExpression(gene2tx[gene], idsCtrl, idsCase))
			expression[gene].addTx(tx, ctrl, case)

			if expression[gene].isComplete:
				switches = expression[gene].detectSwitches(minExpression)

				for (nTx,tTx),samples in switches.items():
					thisSwitch = IsoformSwitch(nTx, tTx, samples)
					nInfo = self._txs.nodes()[nTx]
					tInfo = self._txs.nodes()[tTx]
					thisSwitch.addTxInfo(nInfo, tInfo)

					self._genes.update_node("switches", thisSwitch, full_name = gene)

				expression.pop(gene)

	def readSamples(self, FILE):

		ids = FILE.readline().strip().split('\t')
		ids.pop(0)

		return ids

	def getGene2Tx(self):

		gene2tx = {}

		pairs = [ (t,i["gene_id"]) for t,i in self._txs.nodes(data=True) ]

		for tx,gene in pairs:
			gene2tx.setdefault(gene, set())
			gene2tx[gene].add(tx)

		return(gene2tx)
