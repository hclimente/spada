from spada.bio import protein, transcript
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

def test_printSwitchesToGff():
	
	g = get_switches.GetSwitches(dataPath + 'annotation.pklz')
	g.run(dataPath + 'switches')

	io.printSwitchesToGff(g._genes, g._txs)

	for line in io.readGTF('switches_spada.gff', tag_sep = '='):

		if line['feature'] == 'mRNA':
			assert line['ID'] in g._txs.nodes()
			assert line['strand'] == g._txs[line['ID']]['strand']
		elif line['feature'] == 'exon':
			assert [line['start'], line['end']] in g._txs[line['Parent']]['exons']
			assert line['strand'] == g._txs[line['Parent']]['strand']
		elif line['feature'] == 'CDS':
			assert (line['start'], line['end']) == tuple(g._txs[line['Parent']]['CDS'])
			assert line['strand'] == g._txs[line['Parent']]['strand']
		elif line['feature'] == 'start_codon':
			if line['strand'] == '+':
				assert line['start'] ==  g._txs[line['Parent']]['start_codon']
			elif line['strand'] == '-':
				assert line['end'] ==  g._txs[line['Parent']]['start_codon']
			assert line['strand'] == g._txs[line['Parent']]['strand']
		elif line['feature'] in ('five_prime_UTR', 'three_prime_UTR'):
			pass
		elif line['feature'] == 'translated_nucleotide_match':
			parent = line['CDS_matches'].replace('cds_', '')
			assert line['strand'] == g._txs[parent]['strand']
			p = protein.Protein(parent, g._txs[parent])
			regions = p.getFeature(line['Note'], line['Alias'])
			gpos = map(lambda x: [ y.genomicPosition for y in x ], regions)
			bounds = [ (min(x), max(x)) for x in gpos ]

			assert (line['start'], line['end']) in bounds

	os.remove('switches_spada.gff')