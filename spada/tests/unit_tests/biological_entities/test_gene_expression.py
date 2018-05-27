from spada.biological_entities.gene_expression import GeneExpression

import numpy as np
import pytest

txs = {'tx1', 'tx2', 'tx3'}
ctrl_ids = ['X','Y','Z','A']
case_ids = ['A','Z','B','C']

def test_init():

	g = GeneExpression(txs, ctrl_ids, case_ids)

	assert g._allTxs == txs
	assert g._storedTxs == []
	assert g._expressionCtrl.size == 0
	assert g._expressionCase.size == 0

def test_addTx():

	expr1 = np.array([[1.2, 3.4, 5.6]])
	expr2 = np.array([[2.1, 4.3, 6.5]])

	g = GeneExpression(txs, ctrl_ids, case_ids)
	g.addTx('tx1', expr1, expr2)
	assert np.all(g._expressionCtrl == expr1)
	assert g._expressionCtrl.shape == (1,3)
	assert np.all(g._expressionCase == expr2)
	assert g._expressionCase.shape == (1,3)
	assert g._storedTxs == ['tx1']

	g.addTx('tx2', expr2, expr1)
	assert np.all(g._expressionCtrl[1,:] == expr2)
	assert g._expressionCtrl.shape == (2,3)
	assert np.all(g._expressionCase[1,:] == expr1)
	assert g._expressionCase.shape == (2,3)
	assert g._storedTxs == ['tx1', 'tx2']

def test_isComplete():

	expr = np.array([[2, 3, 4, 5]])
	g = GeneExpression(txs, ctrl_ids, case_ids)

	assert not g.isComplete
	g.addTx('tx1', expr, expr)
	assert not g.isComplete
	g.addTx('tx2', expr, expr)
	assert not g.isComplete
	g.addTx('tx3', expr, expr)
	assert g.isComplete

def test_computePSI():

	expr = np.array([[2, 3, 4, 5]])

	g = GeneExpression(txs, ctrl_ids, case_ids)

	g.addTx('tx1', expr, expr)
	psi = g.computePSI(g._expressionCtrl)
	assert np.all(psi == np.array([[1, 1, 1, 1]]))

	g.addTx('tx2', expr, expr)
	psi = g.computePSI(g._expressionCtrl)
	assert np.all(psi == np.array([[.5, .5, .5, .5],[.5, .5, .5, .5]]))

	g.addTx('tx3', np.array([[6, 4, 2, 0]]), expr)
	psi = g.computePSI(g._expressionCtrl)
	assert np.all(psi == np.array([[.2, .3, .4, .5],
								   [.2, .3, .4, .5],
								   [.6, .4, .2, 0]]))

def test_matchExpressions():

	expr1 = np.array([[2, 3, 4, 5]])
	expr2 = np.array([[6, 4, 2, 0]])

	g1 = GeneExpression(txs, ctrl_ids, ctrl_ids)
	g1.addTx('tx1', expr1, expr2)
	g1.addTx('tx2', expr1, expr2)
	g1.addTx('tx3', expr2, expr1)
	psiCtrl = g1.computePSI(g1.matchExpressions(g1._expressionCtrl))
	psiCase = g1.computePSI(g1._expressionCase)
	dpsi = psiCase - psiCtrl

	assert dpsi.shape == (3, 4)
	assert np.all(dpsi[:,0] == [3/7  - .2, 3/7  - .2, 1/7  - .6])
	assert np.all(dpsi[:,1] == [4/11 - .3, 4/11 - .3, 3/11 - .4])
	assert np.all(dpsi[:,2] == pytest.approx([-.15, -.15, .3], .0001))
	assert np.all(dpsi[:,3] == pytest.approx([-.5, -.5, 1], .0001))

	g2 = GeneExpression(txs, ctrl_ids, case_ids)
	g2.addTx('tx1', expr1, expr2)
	g2.addTx('tx2', expr1, expr2)
	g2.addTx('tx3', expr2, expr1)
	psiCtrl = g2.computePSI(g2.matchExpressions(g2._expressionCtrl))
	psiCase = g2.computePSI(g2._expressionCase)
	dpsi = psiCase - psiCtrl

	assert dpsi.shape == (3, 4)
	assert np.all(dpsi[:,0] == [3/7  - .5,  3/7  - .5,  1/7  - .0])
	assert np.all(dpsi[:,1] == [4/11 - .4,  4/11 - .4,  3/11 - .2])
	assert np.all(dpsi[:,2] == pytest.approx([1/4 - .35, 1/4 - .35, 1/2 - .3], .0001))
	assert np.all(dpsi[:,3] == pytest.approx([0	  - .35, 0   - .35, 1   - .3], .0001))

def test_computeWTDeltaPSI():

	g = GeneExpression(txs, ctrl_ids, case_ids)

	expr1 = np.array([[2, 3, 4, 5]])
	expr2 = np.array([[6, 4, 2, 0]])

	g.addTx('tx1', expr1, expr2)
	g.addTx('tx2', expr1, expr2)
	g.addTx('tx3', expr2, expr1)

	psi_ctrl = np.array([[.5, .3, .4],
						 [.3, .5, .15]])
	wt_dpsi = g.computeWTDeltaPSI(psi = psi_ctrl)
	assert wt_dpsi.shape == (2, 3)
	assert np.all(wt_dpsi == pytest.approx(np.array([[.2, .1, .1],
													 [.2, .15, .35]]), .0001))

def test_cutoff():

	g = GeneExpression(txs, ctrl_ids, case_ids)

	x = np.array([[0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1],
				  [0,  0,  0,  0,  0,  1,  1,  1,  1,  1, 1]])
	assert np.all(g.cutoff(x, 45) == [[.45], [.5]])
	assert np.all(g.cutoff(x, 80) == [[.8], [1]])
	assert np.all(g.cutoff(x, 16.34) == pytest.approx(np.array([[.1634], [0]]), .0001))

def test_detectSwitches():

	expression = np.array([[0,  .5,  1, 2],
						   [2, .05, .4, 5],
						   [8,   3, .1, 0]])
	dpsi = np.array([[  1,- .5, .5,  0],
					 [ .5,-  1,  1,  1],
				 	 [-.3,   1,-.25, 0]])
	wtdpsi = np.array([[  0, 0, .15, .25, .05, .11, .09],
					   [  0, 0, .15, .25, .05, .11, .09],
				 	   [  0, 0, .15, .25, .05, .11, .09]])

	g = GeneExpression(txs, ctrl_ids, case_ids)
	g._storedTxs = ['tx1', 'tx2', 'tx3']
	g._complete = True
	g._matchedExpressionCtrl = expression
	g._expressionCase = expression
	g._wtdPSI = wtdpsi
	g._dPSI = dpsi

	switches = g.detectSwitches()

	assert switches[('tx3','tx2')] == set(['A','B'])
	assert switches[('tx1','tx3')] == set('Z')
