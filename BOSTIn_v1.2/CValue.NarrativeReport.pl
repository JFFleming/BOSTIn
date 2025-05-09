#!/usr/bin/perl
use strict;
use warnings;
my $prefix = $ARGV[0];
my $summInput = "$prefix.SiteSaturation.TotalCScore.txt";
open(my $summFile, '<', $summInput);
my $fulldummy = <$summFile>;

my %stats;

while(<$summFile>){
	my $line = $_;
	chomp($line);
	my @splitter = split("\t", $line);
	my $statValue = $splitter[1];
	$stats{$splitter[0]} = $statValue;
}
close($summFile);

print "C Score
	As you are using nucleotide data, we selected the C Score - the Convergence Factor. This is based on the ratio of the standard deviation of the Transition/Transversion ratio of the dataset and the standard deviation of the uncorrected p-distance. It also has the advantage that it can be used to detect taxa that are particularly saturated. We call this the taxon C Score.
	The C Score of this entire dataset is $stats{CScore}, where the standard deviation of the transition transversion ratios is $stats{StdDevTiTv} and the standard deviation of the uncorrected p-distances is $stats{PValue}.";

if ($stats{CScore} > 20){
 print "As the C Score is greater than 20, this means that saturation is unlikely to be a problem with this dataset, as the Observed and Expected Ti/Tv ratio is likely greater than 1 (Struck et al 2008).\n";
}
elsif ($stats{CScore} > 10){
	print "As the C Score is between 10 and 20, it means that this dataset is not entirely saturated, but that saturation might affect topological reconstruction. As it is still quite low, you might want to observe the distribution of the taxa C Scores, to see if any taxa are contributing particularly to the overall score.\n";
}
else{
	print "As the C Score is below 10, this means that saturation might be a serious problem in this dataset. Check the taxa C Scores to see which taxa are contributing to this.\n";
}

my $t_summary = "$prefix.SiteSaturation.TaxaCScore.txt";
my @t_stats = summaryStats($t_summary);

my $tC_mean = $t_stats[0];
my $tC_median = $t_stats[1];
my $tC_sd = $t_stats[2]; 
my $tC_20 = $t_stats[4];
my $tC_10 = $t_stats[5];
my $tC_lowsig = $tC_mean-(2*$tC_sd);
my $tC_highsig = $tC_mean+(2*$tC_sd);

my $red_tfile = "$prefix.SiteSaturation.TaxaCScore.redflags.txt";
my @tC_red = redLorryYellow($red_tfile);

my $yellow_tfile = "$prefix.SiteSaturation.TaxaCScore.yellowflags.txt";
my @tC_yellow = redLorryYellow($yellow_tfile);

print "Taxa C Score
	A histogram of taxa C Scores can be seen here:
	
	$prefix.SiteSaturation.TaxaFrequencies.histogram.pdf
	
	The C Score is unique among the BOSTIn metrics. in that it has a more absolute way of assessing saturation. Sequences with a C Score below 20 are potentially saturated to an uninformative degree, whereas taxa with a C Score below 10 show significant amounts of saturation.
	The Percentage of taxa in your dataset with a C Score below 20 is $tC_20.
	The Percentage of taxa in your dataset with a C Score below 10 is $tC_10.\n";
	
print "From here, we can identify those potentially problematic sequences, assigning Yellow Flags to those with C Scores between 10 and 20, and Red Flags to those less than 10.
	The Red Flag taxa are:
	@tC_red
	The Yellow Flag taxa are:
	@tC_yellow
";

sub summaryStats {
	open (my $input, '<', $_[0]);
	my @summary;
	my $dummy = <$input>;
	while(<$input>){
		my $line = $_;
		my @splitSum = split(/\t/, $line);
		my $stat = $splitSum[1];
		chomp($stat);
		my $rounded = sprintf('%.5f', $stat);
		push (@summary, $rounded);
	}
	return @summary;
	close($input);
}

sub redLorryYellow {
	open(my $input, '<', $_[0]);
	my @lorries;
	while(<$input>){
		my $line = $_;
		my @splitLorry = split(/\t/, $line);
		my $tax_id = $splitLorry[0];
		chomp($tax_id);
		push (@lorries, $tax_id);
	}
	return @lorries;
	close($input)
}
