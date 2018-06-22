from spada.methods import method

class GetSwitches(method.Method):
	def __init__(self, annotation = 'annotation.pklz'):
		method.Method.__init__(self, __name__, annotation)
		self._genes.flushSwitches()

	def run(self, switchesFile):

		self._genes.readSwitches(switchesFile, self._txs)

if __name__ == '__main__':
	pass
