#!/usr/bin/perl
use strict;
use warnings;
my $prefix = $ARGV[0];
my $summInput = "$prefix.SiteSaturation.TotalCValue.txt";
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

print "C Factor
	As you are using nucleotide data, we selected the C Factor - the Convergence Factor. This is based on the ratio of the standard deviation of the Transition/Transversion ratio of the dataset and the standard deviation of the uncorrected p-distance. It also has the advantage that it can be used to detect taxa that are particularly saturated. We call this the taxon C Factor.
	The C Factor of this entire dataset is $stats{CValue}, where the standard deviation of the transition transversion ratios is $stats{StdDevTiTv} and the standard deviation of the uncorrected p-distances is $stats{PValue}.";

if ($stats{CValue} > 20){
 print "As the C Factor is greater than 20, this means that saturation is unlikely to be a problem with this dataset, as the Observed and Expected Ti/Tv ratio is likely greater than 1 (Struck et al 2008).\n";
}
elsif ($stats{CValue} > 10){
	print "As the C Factor is between 10 and 20, it means that this dataset is not entirely saturated, but that saturation might affect topological reconstruction. As it is still quite low, you might want to observe the distribution of the taxa C Factors, to see if any taxa are contributing particularly to the overall score.\n";
}
else{
	print "As the C Factor is below 10, this means that saturation might be a serious problem in this dataset. Check the taxa C Factors to see which taxa are contributing to this.\n";
}

my $t_summary = "$prefix.SiteSaturation.TaxaCValue.summary.txt";
my @t_stats = summaryStats($t_summary);

my $tC_mean = $t_stats[0];
my $tC_median = $t_stats[1];
my $tC_sd = $t_stats[2]; 
my $tC_20 = $t_stats[4];
my $tC_10 = $t_stats[5];
my $tC_lowsig = $tC_mean-(2*$tC_sd);
my $tC_highsig = $tC_mean+(2*$tC_sd);

#print "$tC_mean $tC_median $tC_sd";

my $red_tfile = "$prefix.SiteSaturation.TaxaCValue.redflags.txt";
my @tC_red = redLorryYellow($red_tfile);

my $yellow_tfile = "$prefix.SiteSaturation.TaxaCValue.yellowflags.txt";
my @tC_yellow = redLorryYellow($yellow_tfile);

print "Taxa C Factor
	A histogram of taxa C Factors can be seen here:
	
	$prefix.SiteSaturation.TaxaFrequencies.histogram.pdf
	
	The C Factor is unique among the BOSTIn metrics. in that it has a more absolute way of assessing saturation. Sequences with a C Factor below 20 are potentially saturated to an uninformative degree, whereas taxa with a C Factor below 10 show significant amounts of saturation.
	The Percentage of taxa in your dataset with a C Factor below 20 is $tC_20.
	The Percentage of taxa in your dataset with a C Factor below 10 is $tC_10.\n";
	
print "From here, we can identify those potentially problematic sequences, assigning Yellow Flags to those with C Factors between 10 and 20, and Red Flags to those less than 10.
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
#		print $stat;
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
