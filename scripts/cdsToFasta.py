#!/soft/devel/python-2.7/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

equal = 0
notOrf = 0
badCoordinates = 0

with open("knownGene.txt","r") as CDSFile, open("proteins_cds.fa","w") as PROTEINS, open("transcripts_cds.fa","w") as TRANSCRIPTS:
	currentChr = ""
	chromosome = ""
	for line in CDSFile:
		elements=line.strip().split("\t")

		name		=elements[0] 
		chrom		=elements[1] 
		strand		=elements[2] 
		txStart		=int(elements[3])
		txEnd		=int(elements[4])
		cdsStart	=int(elements[5])
		cdsEnd		=int(elements[6])
		exonCount	=int(elements[7])
		exonStarts	=map(int, filter(None, elements[8].split(",") ) )
		exonEnds	=map(int, filter(None, elements[9].split(",") ) )
		proteinID	=elements[10]
		alignID		=elements[11]
		
		if chrom != currentChr:
			currentChr = chrom
			print currentChr

			CHROMOSOME = open("/projects/rg/TCGA/download/UCSC/hg19/fasta/" + currentChr + ".fa", "rU")
			chromosome = SeqIO.to_dict(SeqIO.parse(CHROMOSOME, "fasta"))[currentChr]
			CHROMOSOME.close()

		if cdsStart == cdsEnd:
			print(name + " couldn't be processed: cdsStart == cdsEnd")
			equal += 1
			continue
		
		orf = ""
		cds = chromosome[cdsStart:cdsEnd]
		norExonStarts = map(lambda x: x - cdsStart, exonStarts)
		norExonEnds = map(lambda x: x - cdsStart, exonEnds)

		if strand == "-":
			cds = cds.reverse_complement()
			norExonStarts = reversed( map(lambda x: cdsEnd - x, exonEnds) )
			norExonEnds = reversed( map(lambda x: cdsEnd - x, exonStarts) )

		for exonSt, exonEnd in zip(norExonStarts, norExonEnds):
			if exonEnd < 0 or exonSt > len(cds):
				print(name + " couldn't be processed: exonEnd < 0 or exonSt > len(cds)")
				badCoordinates += 1
				continue
			
			if not orf:
				exonSt = 0
			if exonEnd > len(cds):
				exonEnd = len(cds) - 1
			
			orf += str(cds[exonSt:exonEnd].seq)
	
		if not orf:
			print(name + " couldn't be processed: not orf")
			notOrf += 1
			continue
		record = SeqRecord(Seq(orf, generic_dna), id=name)

		PROTEINS.write(">" + record.id + "\n")
		PROTEINS.write(str(record.seq.translate(stop_symbol="")) + "\n")

		TRANSCRIPTS.write(">" + record.id + "\n")
		TRANSCRIPTS.write(orf + "\n")

print("Equal: " + str(equal))
print("not orf: " + str(notOrf))
print("bad coordinates: " + str(badCoordinates))