from spada.biological_entities import protein
from spada.io import io
from spada.methods import get_switches, method

import pickle
import pytest
import os

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/../../data/"

def test_printSwitches():

	g = get_switches.GetSwitches(dataPath + 'annotation.pklz')
	g.run(dataPath + 'switches')

	io.printSwitches(g._genes, g._txs)

	switches = [ x for x in io.readTable("switches_spada.tsv") ]

	assert len(switches) == 8

	os.remove("switches_spada.tsv")

def test_getGene2Tx():

	m = method.Method('test_getGene2Tx', dataPath + 'annotation.pklz')
	gene2tx = io.getGene2Tx(m._txs)

	assert len(m._genes.nodes()) == len(gene2tx)
	assert len([ t for t,i in m._txs.transcripts() if i['canonical'] ]) == len([ t for g,T in gene2tx.items() for t in T ])

def test_readSamples():

	EXPR = open(dataPath + 'expression', "r")
	ids = io.readSamples(EXPR)

	assert ids == ['1','2','3']

def test_printSwitchesToGtf():
	
	g = get_switches.GetSwitches(dataPath + 'annotation.pklz')
	g.run(dataPath + 'switches')

	switches = set([ s for g,i,s in g._genes.switches(g._txs)])

	io.printSwitchesToGtf(g._genes, g._txs, switches)

	for line in io.readGTF('switches_spada.gtf'):

		if line['feature'] == 'transcript':
			assert line['transcript_id'] in g._txs.nodes()
			assert line['strand'] == g._txs[line['transcript_id']]['strand']
		elif line['feature'] == 'exon':
			assert [int(line['start']), int(line['end'])] in g._txs[line['transcript_id']]['exons']
			assert line['strand'] == g._txs[line['transcript_id']]['strand']
		else:
			assert line['strand'] == g._txs[line['transcript_id']]['strand']
			p = protein.Protein(line['transcript_id'], g._txs[line['transcript_id']])
			regions = p.getFeature(line['feature'], line['{}_id'.format(line['feature'])])
			gpos = map(lambda x: [ y.genomicPosition for y in x ], regions)
			bounds = [ (min(x), max(x)) for x in gpos ]

			if line['strand'] == '+':
				assert (int(line['start']), int(line['end'])) in bounds
			else:
				assert (int(line['end']), int(line['start'])) in bounds

	os.remove('switches_spada.gtf')