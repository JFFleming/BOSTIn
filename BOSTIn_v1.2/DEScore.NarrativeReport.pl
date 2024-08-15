#!/usr/bin/perl
use strict;
use warnings;

#Take the output Prefix as the command line argument to parse downstream inputs and make outputs.
my $prefix = $ARGV[0];
my $summInput = "$prefix.SiteSaturation.TotalFrequencies.txt";
open(my $summFile, '<', $summInput);
my $fulldummy = <$summFile>;

my $totalDE;
my $TIFreq;
#Read the Summary file produced by the DE Score calculator and extract the relevant output values.
while(<$summFile>){
	my $summline = $_;
	chomp($summline);
	my @splitter = split(/\t/, $summline);
	$totalDE = $splitter[4];
	$TIFreq = $splitter[1];
}
close($summFile);

print "DE-Score
	As you are using amino acid data, we selected the DE-Score - the Dayhoff Exchange Score. This is based on the ratio of pairwise within Dayhoff Category changes and between Dayhoff Category changes throughout the alignment. It also has the advantage that it can be used to detect taxa that are particularly saturated. We call this the tDE-Score.
	The DE-Score of this entire dataset is $totalDE, with a Dayhoff Category Exchange Ratio of $TIFreq. ";

if ($totalDE > 0.5){
 print "As the DE-Score is very high, this means that saturation is unlikely to be a problem with this dataset.\n";
}
elsif ($totalDE > 0.3){
	print "As the DE-Score is above 0.3 but below 0.5, it means that this dataset is not entirely saturated, but as it is still quite low, so you might want to observe the distribution of the tDE-Scores, to see if any taxa are contributing particularly to the score.\n";
}
else{
	print "As the DE-Score is below 0.3, this means that saturation might be a serious problem in this dataset. Check the tDE-Scores to see which taxa are contributing to this.\n";
}
#Open the R Summary file and extract the relevant taxon specific files from the R Summary output.
my $t_summary = "$prefix.SiteSaturation.TaxaFrequencies.summary.txt";
my @t_stats = summaryStats($t_summary);

my $tDE_mean = $t_stats[0];
my $tDE_median = $t_stats[1];
my $tDE_sd = $t_stats[2]; 
my $tDE_lowsig = $tDE_mean-(2*$tDE_sd);
my $tDE_highsig = $tDE_mean+(2*$tDE_sd);

#Use the Red Lorry Yellow Lorry subroutine to extract the Red Flag and Yellow Flag taxa from the txt files produced by the R script.
my $red_tfile = "$prefix.SiteSaturation.TaxaFrequencies.SiteSat.Taxa.redflags.txt";
my @tDE_red = redLorryYellow($red_tfile);

my $yellow_tfile = "$prefix.SiteSaturation.TaxaFrequencies.SiteSat.Taxa.yellowflags.txt";
my @tDE_yellow = redLorryYellow($yellow_tfile);

print "tDE-Score
	A histogram of tDE-Scores can be seen here:
	
	$prefix.SiteSaturation.TaxaFrequencies.histogram.pdf
	
	The difference between the median, $tDE_median ,  and the mean, $tDE_mean , can tell us about the distribution of tDE-Scores, and help us identify taxa that are particularly prone to causing topological artifacts due to site saturation.";


if ($tDE_median > $tDE_highsig){	
	print "The mean is significantly smaller than the median, which suggests that some sequences in the dataset are very significantly saturated. Explore the histogram to learn more, and consider removing the red flag taxa.";
}
elsif ($tDE_median > $tDE_mean){	
	print "The mean is smaller than the median, which suggests that some sequences in the dataset are significantly saturated, but the mean and median are still quite close. This might not be cause for concern.";
}
elsif ($tDE_median < $tDE_lowsig){
	print "The median is significantly smaller than the mean, which suggests that particularly saturated taxa might not cause issues in your dataset.";
}
elsif ($tDE_median < $tDE_mean){
	print "The median is smaller than the mean, which suggests that saturated taxa might not be an issue, but check the histogram and red flag taxa (if any) to see if there are one or two outliers in your dataset.";
}
else{
	print "Your median and mean are equivalent, which suggests that saturation is well distributed across your dataset.";
}
print "From the median, we can identify potentially problematic sequences by observing how many standard deviations and median absolute deviations an individual tDE-Score is from the median. Red flag taxa are more than 2 standard deviations from the median, while Yellow Flag taxa are more than 2 median absolute deviations from the median.
	The Red Flag taxa are:
	@tDE_red
	The Yellow Flag taxa are:
	@tDE_yellow
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
