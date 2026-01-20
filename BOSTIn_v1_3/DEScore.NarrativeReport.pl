#!/usr/bin/perl
use strict;
use warnings;
my $prefix = $ARGV[0];
my $summInput = "$prefix.SiteSaturation.TotalFrequencies.txt";
open(my $summFile, '<', $summInput);
my $fulldummy = <$summFile>;

my $totalDE;
my $TIFreq;
my $crittDE;

while(<$summFile>){
	my $summline = $_;
	chomp($summline);
	my @splitter = split(/\t/, $summline);
	if ($summline =~ /^DE-Score/){
		$totalDE = $splitter[1];
		}
	if ($summline =~ /^Dayhoff Category Exchange Frequency/){
	$TIFreq = $splitter[1];
	}
	if ($summline =~ /^Critical DE-Score/){
	$crittDE = $splitter[1];
	}
}
close($summFile);
$crittDE=sprintf("%.5f", $crittDE);

print "DE-Score
	As you are using amino acid data, we selected the DE-Score - the Dayhoff Exchange Score. This is based on the ratio of pairwise within Dayhoff Category changes and between Dayhoff Category changes throughout the alignment. It also has the advantage that it can be used to detect taxa that are particularly saturated. We call this the tDE-Score.
	The DE-Score of this entire dataset is $totalDE, with a Dayhoff Category Exchange Ratio of $TIFreq.
	The DE-Score is the average of the taxon-specific DE-Scores.
        When evaluating the total DE-Score, consider two things. First, a DE-Score below 0 indicates that the data is completely saturated, and likely uninformative - for best results, seek a DE-Score of 0.354 or higher. This means that the alignment is 2 minimum-information steps away from total saturation (2x0.177 away from DE-Score 0). When comparing between alignments, the one with a higher DE-Score has less entropic site saturation.";

if ($totalDE > 0.353){
 print "As the overall DE-Score is higher than 0.354, this means that saturation is unlikely to be a problem with this dataset overall, but pockets of entropic saturation may still be an issue. Look at the tDE-Scores to check.\n";
}
elsif ($totalDE < 0.354){
	print "As the DE-Score is less than 0.354, saturation might be a serious problem in this dataset. Check the tDE-Scores to see which taxa are contributing to this.\n";
}

my $t_summary = "$prefix.SiteSaturation.TaxaFrequencies.summary.txt";
my @t_stats = summaryStats($t_summary);

my $tDE_mean = $t_stats[0];
my $tDE_median = $t_stats[1];
my $tDE_sd = $t_stats[2];
my $tDE_lowsig = $tDE_mean-(2*$tDE_sd);
my $tDE_highsig = $tDE_mean+(2*$tDE_sd);

my $red_tfile = "$prefix.SiteSaturation.TaxaFrequencies.SiteSat.Taxa.redflags.txt";
my @tDE_red = redLorryYellow($red_tfile);

my $yellow_tfile = "$prefix.SiteSaturation.TaxaFrequencies.SiteSat.Taxa.yellowflags.txt";
my @tDE_yellow = redLorryYellow($yellow_tfile);

print "tDE-Score
	A histogram of tDE-Scores can be seen here:
	
	$prefix.SiteSaturation.TaxaFrequencies.histogram.pdf
        For further analysis, the histograms show the distribution of taxon-specific DE-Scores across the taxa in your dataset, which can help you work out which taxa might end up causing problems in your final phylogenetic analysis.
        We can use the median and the mean of the tDE-Score to better understand what the distribution of site saturation in your dataset looks like.
        tDE-Score median: $tDE_median
	tDE-Score mean: $tDE_mean
";

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
print "Beyond this, we can identify potentially problematic sequences in two ways. First, we can observe which specific taxa have DE-Scores below 0. We mark these as highly saturated Red Flag taxa. However, as we are comparing within the same dataset to seek out problematic taxa, we can also mark taxa with Dayhoff Category Exchange Frequencies below 0.354 as slightly saturated Yellow Flags. We can do this by generating a Critical DE-Score, which is the DE-Score if the Dayhoff Category Exchange Frequencies is 0.354.
	The Critical DE-Score for this dataset is: $crittDE 
	The Red Flag taxa are:\n";
print join("\n", @tDE_red);
print "\nThe Yellow Flag taxa are:\n";
print join("\n", @tDE_yellow);
print "\n";

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
