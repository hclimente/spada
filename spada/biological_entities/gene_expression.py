from spada.io.error import SpadaError

from itertools import combinations
import numpy as np

class GeneExpression:

	def __init__(self, txs, ctrl_ids, case_ids):
		self._allTxs = txs
		self._storedTxs = []

		self._top_ctrl = (None, 0)
		self._top_case = (None, 0)

		self._expressionCtrl = np.array([])
		self._expressionCase = np.array([])
		self._matchedExpressionCtrl = np.array([])

		self._dPSI = np.array([])
		self._wtdPSI = np.array([])
		self._dExp = np.array([])
		self._wtdExp = np.array([])

		self._idsCase = case_ids
		self._idsCtrl = ctrl_ids

		self._complete = False

	def addTx(self, tx, expressionCtrl, expressionCase):

		if self._expressionCtrl.size == 0:
			self._expressionCtrl = np.empty(shape = (0, expressionCtrl.shape[1]))
		self._expressionCtrl = np.vstack((self._expressionCtrl, expressionCtrl))

		medianCtrl = np.median(expressionCtrl)
		if medianCtrl > self._top_ctrl[1]:
			self._top_ctrl = (tx, medianCtrl)

		if self._expressionCase.size == 0:
			self._expressionCase = np.empty(shape = (0, expressionCase.shape[1]))
		self._expressionCase = np.vstack((self._expressionCase, expressionCase))

		medianCase = np.median(expressionCase)
		if medianCase > self._top_case[1]:
			self._top_case = (tx, medianCase)

		self._storedTxs.append(tx)

	@property
	def isComplete(self):

		if not self._complete and not self._allTxs.difference(self._storedTxs):
			self._matchedExpressionCtrl = self.matchExpressions(self._expressionCtrl)
			psiCtrl = self.computePSI(self._matchedExpressionCtrl)
			psiCase = self.computePSI(self._expressionCase)
			self._dPSI = psiCase - psiCtrl
			self._wtdPSI = self.computeExpectedDelta(expression = self._expressionCtrl)

			geneMatchedExpressionCtrl = self._matchedExpressionCtrl.sum(axis = 0)
			geneExpressionCtrl = self._expressionCtrl.sum(axis = 0)
			geneExpressionCase = self._expressionCase.sum(axis = 0)

			self._dExp = geneExpressionCase - geneMatchedExpressionCtrl
			self._wtdExp = self.computeExpectedDelta(psi = geneExpressionCtrl)

			self._complete = True

		return self._complete

	def computePSI(self, expression, nan_rm = False):

		if expression.shape == (0,):
			raise SpadaError("Expression empty.")
		psi = expression / expression.sum(axis = 0)

		if nan_rm:
			nancols = np.where(np.isnan(psi))[1]
			psi = np.delete(psi, nancols, axis = 1)

		return psi

	def computeExpectedDelta(self, expression = None, psi = np.array([]), nan_rm = True):

		if psi.shape == (0,):
			psi = self.computePSI(expression, nan_rm)

		expectedDelta = [ abs(a - b) for a, b in combinations(psi.T, 2) ]
		expectedDelta = np.array(expectedDelta)
		expectedDelta = expectedDelta.T

		return expectedDelta

	def matchExpressions(self, expressionCtrl):

		mask_case = [ i in self._idsCtrl for i in self._idsCase ]
		mask_ctrl = [ self._idsCtrl.index(i) for i in self._idsCase if i in self._idsCtrl ]

		median = np.median(expressionCtrl, axis=1)
		matched = np.tile(median, (len(self._idsCase),1)).T
		matched[:,mask_case] = expressionCtrl[:,mask_ctrl]

		return matched

	def cutoff(self, wt, percent):

		if wt.shape == (0,):
			cutoff = np.inf
		else:
			cutoff = np.percentile(wt, percent, axis = 1)
			cutoff = cutoff.reshape(cutoff.shape[0], 1)

		return cutoff

	def detectSwitches(self, minExpression = 0.1):

		switches = {}

		if self.isComplete:

			bigChange = abs(self._dPSI) > self.cutoff(self._wtdPSI, 95)
			notDE = abs(self._dExp) < self.cutoff(np.expand_dims(self._wtdExp, 0), 95)
			ctrlExpression = (self._matchedExpressionCtrl >= minExpression) & (self._dPSI < 0)
			caseExpression = (self._expressionCase >= minExpression) & (self._dPSI > 0)
			unchanged = np.logical_not(bigChange & notDE & (ctrlExpression | caseExpression))

			# discarded txs by nan
			dpsi = np.copy(self._dPSI)
			dpsi[unchanged] = np.nan

			# find rows with highest and lowest dpsi
			rank = dpsi.argsort(axis = 0)
			txsCtrl = rank[0,:]
			# nans will go to the bottom. Take the last valid index
			# i.e. #rows - #nans - 1
			lastValid = dpsi.shape[0] - np.isnan(dpsi).sum(axis = 0) - 1
			lastValid[lastValid < 0] = 0
			txsCase = rank[lastValid, np.arange(rank.shape[1])]

			for i in range(len(self._idsCase)):
				txCtrl = self._storedTxs[txsCtrl[i]]
				txCase = self._storedTxs[txsCase[i]]

				if txCtrl == txCase:
					continue

				sample = self._idsCase[i]

				switches.setdefault((txCtrl,txCase), set())
				switches[(txCtrl,txCase)].add(sample)

		return switches
