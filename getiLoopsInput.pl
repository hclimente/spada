#!/soft/devel/perl-5.18.1/bin/perl

print "Hola2";

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Bio::SeqIO;
use File::Copy;
use File::Path 'remove_tree';

print "Hola";
 
my $expressedTranscripts = $ARGV[0];
my $candidateTranscripts = $ARGV[1];
my $getExpressedGenes = $ARGV[2];
open EXPRESSED, $expressedTranscripts or die $!;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
								 -host => 'ensembldb.ensembl.org',
								 -user => 'anonymous'
							    );
$registry->set_reconnect_when_lost();

print $getExpressedGenes;

if($getExpressedGenes){
	foreach my $line (<EXPRESSED>)  {

		my $stable_id = (split(/\./, $line))[0];
		my $transcript_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
		my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);
	
		next if !defined $transcript;
		next if !defined $transcript->translate();
		
		my $out_seq = Bio::SeqIO->new(
		                              -file => '>>Results/iLoops/ExpressedTranscripts.fasta',
		                              -format => 'fasta'
		                             );
		
		my $seq = Bio::Seq->new(
							 -seq => $transcript->translate()->seq(),
                        	 -id  => $stable_id,
							 -accession_number => $stable_id,
							);

		$out_seq->write_seq($seq);
	}
} else {
	copy("old/iLoops/input/ExpressedTranscripts.fasta","Results/iLoops/input/ExpressedTranscripts.fasta") or die "Copy failed: $!";
}

close EXPRESSED or die $!;
open CANDIDATES, $candidateTranscripts or die $!;
open GFF_TRACK, '>Results/candidates.gff' or die $!;

print GFF_TRACK "##gff-version 3\n";

foreach my $line (<CANDIDATES>)  {
	my @candidates = split(/\t/, $line);
	my $delete = 0;
	foreach my $rawCandidate (@candidates){
		
		my $fileNumber = 1;
		my $numberOfCandidates = 0;

		my $candidate 		    = (split(/\./, $rawCandidate))[0];
		my $transcript_adaptor  = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
		my $transcript 		    = $transcript_adaptor->fetch_by_stable_id($candidate);

		if(defined $transcript){
			my @exons_transcript    = @{ $transcript->get_all_Exons() };
			my $nb_exons_transcript = scalar @exons_transcript;
	
			mkdir "Results/iLoops/input/".$candidate;
	
			my $tChromosome	= "chr".$transcript->seq_region_name();
			my $tStart      = $transcript->start();
			my $tEnd        = $transcript->end();
			my $tStrand     = ( $transcript->strand()==1 ? "+" : "-" );
			my $tFrame     	= ".";#$transcript->frame();
			
			#print GFF_TRACK "$tChromosome\t.\tmRNA\t$tStart\t$tEnd\t.\t$tStrand\t$tFrame\tID=$candidate\n";
			foreach my $exon (@exons_transcript){
			    my $exon_id = $exon->stable_id() or die("exon id is required");
				my $start   = $exon->start();
				my $end     = $exon->end();
				my $frame   = $exon->frame();
				print GFF_TRACK "$tChromosome\t.\texon\t$start\t$end\t.\t$tStrand\t$frame\tParent=$candidate\n";
			}
		
			open EXPRESSED, $expressedTranscripts or die $!;
			open PAIRS, '>>Results/iLoops/input/'.$candidate.'/'.$candidate.'_'.$fileNumber.'.net' or die $!;
	
			foreach my $rawExpressed (<EXPRESSED>)  {
				my $expressed = (split(/\./, $rawExpressed))[0];
				print PAIRS $candidate."\t".$expressed."\n";
				$numberOfCandidates++;
				if($numberOfCandidates>=10000){
					$fileNumber++;
					$numberOfCandidates = 0;
					close PAIRS or die $!;
					open PAIRS, '>>Results/iLoops/input/'.$candidate.'/'.$candidate.'_'.$fileNumber.'.net' or die $!;
				}
			}
	
			close EXPRESSED or die $!;
			close PAIRS or die $!;

		} else {
			print "\t*".$candidate." not defined.\n";
			$delete = 1;
		}
	}

	if($delete){
		foreach my $rawCandidate (@candidates){
			my $candidate = (split(/\./, $rawCandidate))[0];
			remove_tree('Results/iLoops/input/'.$candidate) or die "Couldn't delete directory ".$candidate;
		}
	}
}

close CANDIDATES or die $!;
close GFF_TRACK or die $!;