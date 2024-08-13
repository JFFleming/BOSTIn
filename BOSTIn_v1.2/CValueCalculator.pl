#!/usr/bin/perl
use strict;
use warnings;
use FAST::Bio::SeqIO;

my $alignFile = FAST::Bio::SeqIO->new(-file => $ARGV[0], -format => 'Fasta', -alphabet => $ARGV[2]);
my $FreqFile = "$ARGV[1].SiteSaturation.TotalCValue.txt";
my $TaxFile = "$ARGV[1].SiteSaturation.TaxaCValue.txt";
open (FREQ, '>', $FreqFile) || die ("Can not open $FreqFile\n");
print FREQ "Taxa\tStdTiTv\tPValue\tCValue\n";

#print FREQ "Taxa1\tTaxa2\tTi\tTv\tSame\tAmbig\tGap\tTiTv\tPDist\tCValue\n";
open (TAXA, '>', $TaxFile) || die ("Can not open $TaxFile\n");
print TAXA "FileName\tTransitionFreq\tTiTv\tTaxonPValue\tC-Value\n";
my @AllSeqs = ();
my @Pyramidine = ("T", "C");
my @Purine= ("A","G");
my %PhySeqs;

while ( my $seq = $alignFile->next_seq() ) {
	my $ID = $seq->display_id;
	my $subseq = $seq->subseq(1, $seq->length());
	my @sequence = split(//,$subseq);
#	print $sequence[0];
#	print "\n";
	$PhySeqs{$ID}=[@sequence];
}

#foreach my $key (keys %PhySeqs){
#	my @unlock = @{$PhySeqs{$key}};
#	print "$unlock[1]\n";
#}

my $size = keys %PhySeqs;
print "there are $size keys in the hash\n";
my $a = 0;
my @AllPDists;
my @AllTiFreqs;
my @AllTiTv;

foreach my $k (keys %PhySeqs){
	my @vals = @{ $PhySeqs{$k}};
#		print "@vals \n";
	my $b = 0;
	my @TaxPDists;
	my @TaxTiFreqs;
	my @TaxTiTv;
	foreach my $compk (keys %PhySeqs){
		next if ($k eq $compk);
		my @compvals = @{ $PhySeqs{$compk}};
		my $i = 0;
		my $pairwiseTI = 0;
		my $pairwiseTV = 0;
		my $SameCounter = 0;
		my $GapCounter = 0;
		my $AmbigCounter = 0;
		foreach (@vals){
			my $query = $_;
			my $compquery = $compvals[$i];
#				print "$i $compquery \t $query\n";
			if ($query eq "-" || $compvals[$i] eq "-"){
				$GapCounter++;
			}
			elsif ($compquery eq $query){
				$SameCounter++;
			}
			elsif ($compquery ~~ @Purine && $query ~~ @Purine){
				$pairwiseTI++;
#				print "$compquery $query\t";
#				$SmallCounter++;
			}
			elsif ($query ~~ @Pyramidine && $compquery ~~ @Pyramidine){
                $pairwiseTI++;
            }
           	elsif ($query ~~ @Pyramidine && $compquery ~~ @Purine){
                $pairwiseTV++;
            }
            else {
            	$AmbigCounter++;
            }
			$i++;
		}
		$b++;
		my $TiTv;
		my $P;
#			print "LOOP $k SUBLOOP $compk OVER - Pairwise Transitions are $pairwiseTI Pairwise Transversions are $pairwiseTV \n";
		if ($pairwiseTI == 0 || $pairwiseTV == 0) {
			$TiTv = 0;
			$P = 0;

		}
		else{
			$TiTv = $pairwiseTI/$pairwiseTV;
			$P = ($pairwiseTI+$pairwiseTV)/($pairwiseTI+$pairwiseTV+$SameCounter);			
			}
#		print "$k \t $compk $pairwiseTI \t $pairwiseTV \t $SameCounter \t $AmbigCounter \t $GapCounter\t $TiTv\t$P\n";

		push (@AllTiTv, $TiTv);
		push (@TaxTiTv, $TiTv);
		push (@AllPDists, $P);
		push (@TaxPDists, $P);
	}
	$a++;
	my $TaxStdTiTv = get_stddev(\@TaxTiTv);
	my $TaxPValue = get_stddev(\@TaxPDists);
	my $TaxCValue = $TaxStdTiTv/$TaxPValue;
	print "\nLOOP $k OVER \n";
	print TAXA "$k\t$TaxStdTiTv\t$TaxPValue\t$TaxCValue\n";
#	print "$k => @vals \n";
}
#my $checksize = scalar(@AllTiFreqs);
#print $checksize;
#my $AverageTiTv = avg(\@AllTiTv);
my $StdTiTv = get_stddev(\@AllTiTv);
#my $AverageTiFreq = avg(\@AllTiFreqs);
#my $StdTiFreq =  get_stddev(\@AllTiFreqs);
my $PValue = get_stddev(\@AllPDists);
my $CValue = $StdTiTv/$PValue;

print FREQ "$ARGV[0]\t$StdTiTv\t$PValue\t$CValue\n";

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
