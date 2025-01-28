#!/usr/bin/perl
use strict;
use warnings;
use FAST::Bio::SeqIO;

#Read the Fasta File, and open the output files for writing
my $alignFile = FAST::Bio::SeqIO->new(-file => $ARGV[0], -format => 'Fasta', -alphabet => protein);
my $FreqFile = "$ARGV[1].SiteSaturation.TotalFrequencies.txt";
my $TaxFile = "$ARGV[1].SiteSaturation.TaxaFrequencies.txt";
open (FREQ, '>', $FreqFile) || die ("Can not open $FreqFile\n");
print FREQ "FileName\tExchangeFreq\tExchangeFreqStDev\tDE-Score\n";
open (TAXA, '>', $TaxFile) || die ("Can not open $TaxFile\n");
print TAXA "FileName\tExchangeFreq\tExchangeFreqStDev\tDE-Score\n";

#Set the Dayhoff Groups
my @AllSeqs = ();
my @Small= ("A","G","P","S","T");
my @AcidAmide= ("D","E","N","Q");
my @Basic= ("H","K","R");
my @Hydrophobic= ("I","L","V","M");
my @Aromatic= ("F","W","Y");
my @Sulfur= ("C");
my %PhySeqs;

#Read the sequences in, the split into sequence and taxa ID. Then split the sequence into individual characters and store these together in a hash.
while ( my $seq = $alignFile->next_seq() ) {
	my $ID = $seq->display_id;
	my $subseq = $seq->subseq(1, $seq->length());
	my @sequence = split(//,$subseq);
#	print $sequence[0];
#	print "\n";
	$PhySeqs{$ID}=[@sequence];
}

#debugging
#foreach my $key (keys %PhySeqs){
#	my @unlock = @{$PhySeqs{$key}};
#	print "$unlock[1]\n";
#}

#Print the number of taxa in the dataset
my $size = keys %PhySeqs;
print "there are $size keys in the hash\n";
my $a = 0;
my @AllTvs;
my @AllTis;
my @AllTiFreqs;

#Foreach taxa in the dataset, make pairwise comparisons against each character in the sequence one by one, against each other taxa - compk

foreach my $k (keys %PhySeqs){
	my @vals = @{ $PhySeqs{$k}};
#		print "@vals \n";
	my $b = 0;
	my @TaxTvs;
	my @TaxTis;
	my @TaxTiFreqs;
	foreach my $compk (keys %PhySeqs){
		next if ($k eq $compk);
		my @compvals = @{ $PhySeqs{$compk}};
		my $i = 0;
		my $pairwiseTI = 0;
		my $pairwiseTV = 0;
		my $SameCounter = 0;
		my $SmallCounter = 0;
		my $GapCounter = 0;
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
			elsif ($compquery ~~ @Small && $query ~~ @Small){
				$pairwiseTI++;
#				print "$compquery $query\t";
#				$SmallCounter++;
			}
			elsif ($query ~~ @AcidAmide && $compquery ~~ @AcidAmide){
                                $pairwiseTI++;
                        }
                        elsif ($query ~~ @Basic && $compquery ~~ @Basic){
                                $pairwiseTI++;
                        }
                        elsif ($query ~~ @Hydrophobic && $compquery ~~ @Hydrophobic){
                                $pairwiseTI++;
                        }
                        elsif ($query ~~  @Aromatic && $compquery ~~ @Aromatic){
                                $pairwiseTI++;
                        }
                        elsif ($query ~~  @Sulfur && $compquery ~~ @Sulfur){
                                $pairwiseTI++;
                        }
                        else{
                                $pairwiseTV++;
                        }
			$i++;
		}
		$b++;
#			print "LOOP $k SUBLOOP $compk OVER - Pairwise Transitions are $pairwiseTI Pairwise Transversions are $pairwiseTV \n";
		if ($pairwiseTI == 0 && $pairwiseTV == 0) {
#				print "GAP \n";
			$GapCounter++;
		}
		else{
#			print "TI: $pairwiseTI TV: $pairwiseTV \n";
		my $TiFreq = $pairwiseTI/($pairwiseTI+$pairwiseTV);
#		print "TI $pairwiseTI TV $pairwiseTV \t $SmallCounter \t $SameCounter \n";
#Having looked at each taxa for a pairwise comparison, shove the values into a taxa specific array for later and an all values array for later.
		push (@AllTis, $pairwiseTI);
		push (@AllTvs, $pairwiseTV);
		push (@AllTiFreqs, $TiFreq);
		push (@TaxTis, $pairwiseTI);
		push (@TaxTvs, $pairwiseTV);
		push (@TaxTiFreqs, $TiFreq);
		}
	}
	$a++;
 #Calculate Taxa specific Between and Within Dayhoff Group exchange frequencies, standard deviations and the taxa specific DE-Score.
	my $TaxAverageTransI = avg(\@TaxTis);
	my $TaxAverageTransV = avg(\@TaxTvs);
	my $TaxAverageTiFreq = avg(\@TaxTiFreqs);
	my $TaxStdTransI = get_stddev(\@TaxTis);
	my $TaxStdTransV = get_stddev(\@TaxTvs);
	my $TaxStdTiFreq =  get_stddev(\@TaxTiFreqs);
	my $TaxDist = $TaxAverageTiFreq - 0.177;
	my $TaxDE = $TaxDist/(0.255*$size**-0.15);
#	print "\nLOOP $k OVER \n";
	print TAXA "$k\t$TaxAverageTiFreq\t$TaxStdTiFreq\t$TaxDE\n";
#	print "$k => @vals \n";
}
#my $checksize = scalar(@AllTiFreqs);
#print $checksize;
#Calculate same values as above but for the whole dataset.
my $AverageTransI = avg(\@AllTis);
my $AverageTransV = avg(\@AllTvs);
my $AverageTiFreq = avg(\@AllTiFreqs);
my $StdTransI = get_stddev(\@AllTis);
my $StdTransV = get_stddev(\@AllTvs);
my $StdTiFreq =  get_stddev(\@AllTiFreqs);
my $FreqDist = $AverageTiFreq - 0.177;
my $totalDE = $FreqDist/(0.255*$size**-0.15);

print FREQ "$ARGV[0]\t$AverageTiFreq\t$StdTiFreq\t$totalDE\n";

#Subroutine to calculate the average
sub avg {
  my ($avg_array_ref) = @_;
  my $avg_count = @$avg_array_ref;
  my $avg_sum = 0;
  for my $avg_num (@$avg_array_ref) {
      $avg_sum += $avg_num;
  }
  return $avg_sum / $avg_count;
}

#Subroutine to calculate the standard deviation
sub get_stddev {
  return sqrt(get_disp(@_));
}

#Subroutine to calculate the disparity for the standard deviation
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
