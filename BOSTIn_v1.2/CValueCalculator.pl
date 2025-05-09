#!/usr/bin/perl
use strict;
use warnings;
use FAST::Bio::SeqIO;

#Read the Fasta file
my $alignFile = FAST::Bio::SeqIO->new(-file => $ARGV[0], -format => 'Fasta', -alphabet => 'dna');
my $FreqFile = "$ARGV[1].SiteSaturation.TotalCScore.txt";
my $TaxFile = "$ARGV[1].SiteSaturation.TaxaCScore.txt";
#Open output files for writing.
open (FREQ, '>', $FreqFile) || die ("Can not open $FreqFile\n");

open (TAXA, '>', $TaxFile) || die ("Can not open $TaxFile\n");
print TAXA "FileName\tTiTv\tTaxonPValue\tCScore\n";
#Set Purines and Pyrimidines to measure TiTv
my @AllSeqs = ();
my @Pyramidine = ("T","C");
my @Purine= ("A","G");
my @ValidSeqs = ("A","T","C","G","-");
my %PhySeqs;

my %purine = map { $_ => 1 } @Purine;
my %pyramidine = map { $_ => 1 } @Pyramidine;
my %valids = map {$_ => 1} @ValidSeqs;

#Read Fasta File, placing the taxa name and the sequence in separate bins, and then splitting the sequence into single character strings inside the array @sequence, attached together as part of the hash PhySeqs
while ( my $seq = $alignFile->next_seq() ) {
	my $ID = $seq->display_id;
	my $subseq = $seq->subseq(1, $seq->length());
	my @sequence = split(//,$subseq);
	$PhySeqs{$ID}=[@sequence];
}

#How many taxa are in the file?
my $size = keys %PhySeqs;
print "there are $size keys in the hash\n";
my $a = 0;
my @AllPDists;
my @AllTiFreqs;
my @AllTiTv;

#For Each Taxa, do a pairwise comparison against each other taxa - compk. Count Transitions, Transversions, no Change, Gaps and Ambiguous characters.
foreach my $k (keys %PhySeqs){
	my @vals = @{ $PhySeqs{$k}};
	my $b = 0;
	my @TaxPDists;
	my @TaxTiFreqs;
	my @TaxTiTv;
	foreach my $compk (keys %PhySeqs){
		next if ($k eq $compk);
		my @compvals = @{ $PhySeqs{$compk}};
		my $i = 0;
		my $pairwise_ti = 0;
		my $pairwise_tv = 0;
		my $SameCounter = 0;
		my $GapCounter = 0;
		my $AmbigCounter = 0;
		
		foreach (@vals){
	    	my $query = $_;
	    	my $comp_query = $compvals[$i];
	    	if ($valids{$query} && $valids{$comp_query}){
	    	if ($query eq "-" || $comp_query eq "-") {
	   	     	$GapCounter++;
	   		}
	    	elsif ($query eq $comp_query) {
	       	 	$SameCounter++;
	    	}
	    	elsif ($purine{$query} && $purine{$comp_query}) {
                $pairwise_ti++;
            }
            elsif ($pyramidine{$query} && $pyramidine{$comp_query}) {
                $pairwise_ti++;
            }
            elsif ($purine{$query} && $pyramidine{$comp_query}) {
                $pairwise_tv++;
            }
            elsif ($pyramidine{$query} && $purine{$comp_query}) {
                $pairwise_tv++;

            }
	    	else {
	        	$AmbigCounter++;
	    	}
	    	}
	    $i++;
		}
		$b++;
		my $TiTv;
		my $P;

		if ($pairwise_ti == 0 || $pairwise_tv == 0) {
			$TiTv = 0;
			$P = 0;

		}
		else{
			$TiTv = $pairwise_ti/$pairwise_tv;
			$P = ($pairwise_ti+$pairwise_tv)/($pairwise_ti+$pairwise_tv+$SameCounter);			
			}

		push (@AllTiTv, $TiTv);
		push (@TaxTiTv, $TiTv);
		push (@AllPDists, $P);
		push (@TaxPDists, $P);
	}
	$a++;
	my $TaxStdTiTv = get_stddev(\@TaxTiTv);
	my $TaxPValue = get_stddev(\@TaxPDists);
	my $TaxCScore = $TaxStdTiTv/$TaxPValue;
	print "\n Finished assessing $k. Taxon $a of $size \n";
	print TAXA "$k\t$TaxStdTiTv\t$TaxPValue\t$TaxCScore\n";
}
my $StdTiTv = get_stddev(\@AllTiTv);
my $PValue = get_stddev(\@AllPDists);
my $CScore = $StdTiTv/$PValue;

print FREQ "$ARGV[0] C-Score Summary\nStdDevTiTv\t$StdTiTv\nPValue\t$PValue\nCScore\t$CScore\n";

sub avg {
  my ($avg_array_ref) = @_;
  my $avg_count = @$avg_array_ref;
  my $avg_sum = 0;
  for my $avg_num (@$avg_array_ref) {
      $avg_sum += $avg_num;
  }
  return $avg_sum / $avg_count;
}

sub get_stddev {
  return sqrt(get_disp(@_));
}

sub get_disp {
  my ($dis_array_ref) = @_;
  my $mean = avg($dis_array_ref);
  my $dis_count = @$dis_array_ref;
  my $dis_sum = 0;
  
  for my $dis_num (@$dis_array_ref) {
      $dis_sum += (($dis_num - $mean) ** 2);
  }
  return $dis_sum / $dis_count;
}
