#!/usr/bin/perl
use strict;
use warnings;

#Open Input Files
my $LBScoreSummary = "$ARGV[0].LBSummary.txt";
my $LBTaxa = "$ARGV[0].LBi-scores";
my $LBYellow = "$ARGV[0].LBScore.Taxa.yellowflags.txt";
my $LBRed = "$ARGV[0].LBScore.Taxa.redflags.txt";

#Open Whole dataset summary file and extract overall LB Score
open (my $summFile, '<', $LBScoreSummary);
my %stats;
while(<$summFile>){
	my $line = $_;
	chomp($line);
	my @splitter = split("\t", $line);
	my $statValue = $splitter[1];
	$stats{$splitter[0]} = $statValue;
}
close($summFile);

#Open taxa specific summary file and extract taxa specific scores
open (my $taxFile, '<', $LBTaxa);

my %taxaHash;
my $dummy = <$taxFile>;
while(<$taxFile>){
	my $taxline = $_;
	chomp($taxline);
	my @taxsplitter = split("\t",$taxline);
	my $taxValue = $taxsplitter[2];
	$taxaHash{$taxsplitter[0]} = $taxValue;
}
close($taxFile);

#Calculate Red Flag and Yellow Flag taxa
my @yellowFlag;
my @redFlag;
#my $bigsd2 = $stats{Mean}+(2*$stats{standardDeviation});
#my $smallsd2 = $stats{Mean}-(2*$stats{standardDeviation});
my $UQ = $stats{UpperQuartile};
my $UQStdDev = $stats {stdDevUpperQuartile};
my $UQBound1 = $UQ+$UQStdDev;
my $UQBound2 = $UQ+(2*$UQStdDev);
#print $bigsd2, "\n";
#print $smallsd2, "\n";

foreach my $key (keys %taxaHash){
#	print $taxaHash{$key}, "\n";
	if ($taxaHash{$key} > $UQBound2){
		push (@redFlag, $key);
		}
	if ($taxaHash{$key} > $UQBound1){
		push (@yellowFlag, $key);
	}
}
open (yellowFile, '>', $LBYellow);
open (redFile, '>', $LBRed);

foreach (@yellowFlag){
print yellowFile "$_ \t $taxaHash{$_} \n";
}
foreach (@redFlag){
print redFile "$_ \t $taxaHash{$_} \n";
}

print "LB Score
	To identify long branched taxa, BOSTIn uses the LB-score. 
	The LB-Score measures the percentage deviation of each taxon from the average patristic distance, and so is independent of the actual topology of the tree itself, making it quite useful to identify long branches. 
	BostIn rapidly generates a Neighbour-Joining tree to calculate the LB-Score. This produces an LB-Score that is normally significantly similar, even under large amounts of Long Branch Attraction, but it won't be as accurate as an LB-Score generated under the best possible model. 
	For the purposes of defining the sextile of taxa most likely to cause a long branch attraction artifact, however, it ought to suffice. To read more about this, see the BOSTIn manuscript when it appears in pre-print (I'll add a reference here later!)
	The taxa specific LB-Scores in your dataset range from $stats{Minimum} to $stats{Maximum} , with a mean of $stats{Mean} and a standard deviation of $stats{standardDeviation} . 
	
	Typically, LB-Scores identify suspect long-branched taxa by assessing which taxa are outside of one, and then two standard deviations of the mean of the Upper Quartile. This is because it is the extremes of Branch Length heterogeneity that can cause the greatest problems.
	Your Upper Quartile starts at $stats{UpperQuartile}, with your Lower Quartile at $stats{LowerQuartile}.
	We've identified these as yellow flags and red flags respectively, as with the other measurements in BOSTIn.
	Your Red Flag Taxa are:
	@redFlag
	Your Yellow Flag Taxa are:
	@yellowFlag
";	
