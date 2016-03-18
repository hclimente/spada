from libs import options
from libs import utils
from methods import method

class Test(method.Method):
	def __init__(self,gn_network,tx_network,gn_subnetwork=False):
		method.Method.__init__(self, __name__, gn_network, tx_network, gn_subnetwork)

	def run(self):
		from biological_entities import switch

		s= switch.IsoformSwitch( 'uc002pwd.2',  'uc002pwf.2', self._gene_network._net.node["147645"]["isoformSwitches"][0]["patients"])
		nInfo = tx._net.node[s.nTx]
		tInfo = tx._net.node[s.tTx]
		s.addTxs(nInfo,tInfo)
		s.addIsos(nInfo,tInfo,True)

		import pdb
		pdb.set_trace()