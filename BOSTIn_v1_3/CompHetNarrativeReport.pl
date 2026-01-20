#!/usr/bin/perl
use strict;
use warnings;

#This Opening section opens the nRCFV file, searches for the nRCFV and prints it to the report file.
my $nRCFV;
my $RCFV_File;
my $prefix = $ARGV[0];
my $datatype = $ARGV[1];
open ($RCFV_File, '<', "$prefix.RCFV.txt");

while(<$RCFV_File>){
	my $line = $_;
#	print $line;
	if ($line =~ /^nRCFV/){
		my @splitLine = split (/\t/, $line);
		$nRCFV = $splitLine[1];
		chomp($nRCFV);
		}
}
close($RCFV_File);
print ("nRCFV \n The total nRCFV of your dataset was $nRCFV. The nRCFV is the average of the taxon-specific nRCFV values and the character-specific nRCFV values, though, so it doesn't really mean a great deal without the broader context of those values, unless you are comparing between different alignments. If you are, the one with a lower nRCFV has less compositional heterogeneity.\n\n");

my $t_summary = "$prefix.ntRCFV.summary.txt";
my @t_stats = summaryStats($t_summary);
#print @t_stats;

#This section uses the red Lorry Yellow subroutine to grab each taxa with a conspicuously high tsRCFV value.
my $red_tfile = "$prefix.ntRCFV.CompHet.Taxa.redflags.txt";
my @tRCFV_red = redLorryYellow($red_tfile);

my $yellow_tfile = "$prefix.ntRCFV.CompHet.Taxa.yellowflags.txt";
my @tRCFV_yellow = redLorryYellow($yellow_tfile);

print("ntsRCFV
\tA histogram of your taxon-specific nRCFV (ntRCFV) values can be found here:
	
$prefix.ntRCFV.histogram.pdf
	
\tThis shows the distribution of tsnRCFV values across the taxa in your dataset, which, when compared to the average, can help you work out which taxa might end up causing problems in your final phylogenetic analysis. High tsnRCFV values indicate high levels of compositional heterogeneity.
\tWe can use the difference between the median and the mean of the tsnRCFV values to better understand what the distribution of compositional heterogeneity in your dataset looks like.
\tThe mean of your taxon-specific nRCFV is $t_stats[0], while the median is $t_stats[1].");
if ($t_stats[0] > $t_stats[1]){
	print("Your mean is larger than the median, which indicates that your data contains some sequences that are significantly more compositionally heteogeneous than others. This could be a warning sign. \n");
}
elsif ($t_stats[0] < $t_stats[1]){	
	print("Your mean is smaller than the median, which indicates that your data is distributed towards sequences that are compositionally homogenous. This is a positive sign!\n");
}
else{
	print("Your mean and median are equivalent, which suggests that there aren't any great outliers in either direction. This could be a good sign, or could indicate a spread of very heterogeneous and very homogenous sequences. Check your histogram to see whether the distribution looks normal or bifurcating!\n");
}
print("\tFrom the median, we can identify potentially problematic taxa using two metrics, the standard deviation and the median absolute deviation. In a normal distribution, 95% of data should exist within 2 standard deviations of the mean, but single instances of very large values can inflate the mean by quite a bit. If we instead measure 2 standard deviations from the median, we can more comfortably identify any potential outlier. The median absolute deviation, meanwhile, is a bit like the median's version of a standard deviation. It is the median of how far away each value in the dataset is from the median. It is a stricter measurement, and so all taxa that are 2 standard deviations from the median are marked with a Red Flag so that you can seriously consider excluding them from your analysis, whereas all taxa that are 2 median absolute deviations from the median are marked with a Yellow Flag to warn you that they might be a problem.
The Red Flag taxa are:
@tRCFV_red
The Yellow Flag taxa are:
@tRCFV_yellow\n\n");

#This section uses the red Lorry Yellow subroutine to grab each character with a conspicuously high csRCFV value.

my @cRCFV_red;
my $red_cfile = "$prefix.ncsRCFV.CompHet.Character.redflags.txt";
@cRCFV_red = redLorryYellow($red_cfile);

my $yellow_cfile = "$prefix.ncsRCFV.CompHet.Character.yellowflags.txt";
my @cRCFV_yellow = redLorryYellow($yellow_cfile);

my $c_summary = "$prefix.ntRCFV.summary.txt";
my @c_stats = summaryStats($c_summary);


print("ncsRCFV
	A histogram of your character-specific nRCFV (ncsRCFV) values can be found here:
	
	$prefix.ncsRCFV.histogram.pdf
	
	This shows the distribution of csnRCFV values across the characters in your dataset, which, when compared to the average, can help you work out which characters might end up causing problems in your final phylogenetic analysis. High csnRCFV values indicate high levels of compositional heterogeneity.
	We can use the difference between the median and the mean of the csnRCFV values to better understand what the distribution of compositional heterogeneity in your dataset looks like."); 

if ($c_stats[0] > $c_stats[1]){
	print("Your mean is larger than the median, which indicates that your data contains some sites that are significantly more compositionally heteogeneous than others. This could be a warning sign.");
}
elsif ($c_stats[0] < $c_stats[1]){
	print("Your mean is smaller than the median, which indicates that your data is distributed towards sites that are compositionally homogenous. This is a positive sign!");
}
else{
	print("Your mean and median are equivalent, which suggests that there aren't any great outliers in either direction. This could be a good sign, or could indicate a spread of very heterogeneous and very homogenous characters. Check your histogram to see whether the distribution looks normal or bifurcating!\n");
}
print("Like with the tsnRCFV, we can assess the standard deviation and median absolute deviation from the median to establish whether there are any problematic characters. Unfortunately, how to deal with difficult characters is a bit more complex! 
	The Red Flag characters are:
	@cRCFV_red
	The Yellow Flag characters are:
	@cRCFV_yellow \n");
if (scalar(@cRCFV_red) > 4 && $datatype=~"protein"){
	print("\tAs 5 or more characters show high compositional heterogeneity, it might be a good idea to consider an approach such as Dayhoff 6-state recoding, and comparing the trees produced by that process against ones produced by your data without any modification. Be aware, while simplifying the amino acid alphabet down from 20 to 6 characters reduces compositional heterogeneity, it might also mask useful phylogenetic information! Check Hernandez et al (2019)\n");
}
elsif (scalar(@cRCFV_red) > 1 && $datatype=~"dna"){
	print("\tAs 2 or more characters show high compositional heterogeneity, it might be a good idea to consider an approach such as 2-state recoding to purines a pyrimidines, and comparing the trees produced by that process against ones produced by your data without any modification. It might be the result of a significant AT or CG bias in your data, which is sometimes biological, and sometimes artefactual. Be aware, while simplifying the nucleotide alphabet down from 4 to 2 characters reduces compositional heterogeneity, it might also mask useful phylogenetic information! Check Hernandez et al (2019)\n");
}
elsif (scalar(@cRCFV_red) > 0){
 	print("\tHere, only a few characters show high levels of compositional heterogeneity. It might be worth considering a site-sensitive model, such as CAT+GTR for your analysis. You might not need to employ a recoding strategy, as it is possible that you will lose information by masking the variation within your recoded groups.\n");
}
else{
	print("\tThere weren't any particularly problematic characters identified. It doesn't seem like site to site compositional heterogeneity is going to be a problem for your analysis!\n");
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
