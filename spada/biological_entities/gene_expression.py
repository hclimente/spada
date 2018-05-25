from spada.io.io import SpadaError

from itertools import combinations
import numpy as np

class GeneExpression:

	def __init__(self, txs, ctrl_ids, case_ids):
		self._allTxs = txs
		self._storedTxs = []

		self._expressionCtrl = np.array([])
		self._expressionCase = np.array([])

		self._psiCtrl = np.array([])
		self._psiCase = np.array([])

		self._dPSI = np.array([])
		self._wtdPSI = np.array([])

		self._idsCase = case_ids
		self._idsCtrl = ctrl_ids

		self._complete = False

	def addTx(self, tx, expressionCtrl, expressionCase):

		if self._expressionCtrl.size == 0:
			self._expressionCtrl = np.empty(shape = (0, expressionCtrl.shape[1]))
		self._expressionCtrl = np.vstack((self._expressionCtrl, expressionCtrl))

		if self._expressionCase.size == 0:
			self._expressionCase = np.empty(shape = (0, expressionCase.shape[1]))
		self._expressionCase = np.vstack((self._expressionCase, expressionCase))

		self._storedTxs.append(tx)

	@property
	def isComplete(self):

		if not self._complete and not self._allTxs.difference(self._storedTxs):
			self._psiCtrl, self._psiCase, self._dPSI = self.computeDeltaPSI(self._expressionCtrl, self._expressionCase)
			self._wtdPSI = self.computeWTDeltaPSI(self._expressionCtrl, self._psiCtrl)

			self._complete = True

		return self._complete

	def computePSI(self, expression):
		if not expression.shape:
			raise SpadaError("Expression empty.")
		psi = expression / expression.sum(axis = 0)
		return psi

	def computeDeltaPSI(self, expressionCtrl, expressionCase):

		mask_case = [ i in self._idsCtrl for i in self._idsCase ]
		mask_ctrl = [ self._idsCtrl.index(i) for i in self._idsCase if i in self._idsCtrl ]

		psiCtrl = self.computePSI(expressionCtrl)
		psiCase = self.computePSI(expressionCase)

		median_psi_ctrl = np.median(psiCtrl, axis=1)
		median_psi_ctrl = median_psi_ctrl.reshape(len(median_psi_ctrl), 1)

		dpsi = psiCase - median_psi_ctrl
		dpsi[:,mask_case] = psiCase[:,mask_case] - psiCtrl[:,mask_ctrl]

		return (psiCtrl, psiCase, dpsi)

	def computeWTDeltaPSI(self, expression = None, psi = np.array([])):

		if not psi.shape:
			psi = self.computePSI(expression)

		wt_dpsi = [ abs(a - b) for a, b in combinations(psi.T, 2) ]
		wt_dpsi = np.array(wt_dpsi)
		wt_dpsi = wt_dpsi.T

		return wt_dpsi

	def cutoff(self, wtdPSI, percent):

		cutoff = np.percentile(wtdPSI, percent, axis = 1)
		cutoff = cutoff.reshape(cutoff.shape[0],1)

		return cutoff

	def detectSwitches(self, minExpression = 0.1):

		switches = {}

		if self.isComplete:
			bigChange = abs(self._dPSI) > self.cutoff(self._wtdPSI, 95)
			ctrlExpression = (self._expressionCtrl >= minExpression) & (self._dPSI < 0)
			caseExpression = (self._expressionCase >= minExpression) & (self._dPSI > 0)
			unchanged = np.logical_not(bigChange & (ctrlExpression | caseExpression))

			# discarded txs by nan
			dpsi = np.copy(self._dPSI)
			dpsi[unchanged] = np.nan

			# find rows with highest and lowest dpsi
			rank = dpsi.argsort(axis = 0)
			txsCtrl = rank[0,:]
			# nans will go to the bottom. Take the last valid index
			# i.e. #rows - #nans - 1
			lastValid = dpsi.shape[0] - np.isnan(dpsi).sum(axis = 0) - 1
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
