#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $compHet;
my $branchLength;
my $siteSat;
GetOptions('ch' => \$compHet, 'blh' => \$branchLength, 's' => \$siteSat) or die ("Error - no command arguments. \n");

print "Other things found on the command line:\n" if $ARGV[0];
foreach (@ARGV)
{
  print "$_\n";
}

if ($compHet){
	do_compHet();
}

if ($branchLength){
	do_branchLength();
}

if ($siteSat) {
	do_siteSat();
}


sub do_compHet {
	if ($ARGV[0]=~ /dna/){
		print "RCFV Reader DNA!\n";
        system("RCFV_Reader.pl dna $ARGV[1] $ARGV[2]");
		system("nRCFVAnalysis_Script.R $ARGV[2].ntRCFV.txt $ARGV[2].ncsRCFV.txt");
		system("CompHetNarrativeReport.pl $ARGV[2] dna > $ARGV[2].CompositionalHeterogeneity.NarrativeReport.txt");
	}
	elsif ($ARGV[0]=~ /protein/){
		print "RCFV Reader Protein!\n";
		system("RCFV_Reader.pl protein $ARGV[1] $ARGV[2]");
		system("nRCFVAnalysis_Script.R $ARGV[2].ntRCFV.txt $ARGV[2].ncsRCFV.txt");
		system("CompHetNarrativeReport.pl $ARGV[2] protein > $ARGV[2].CompositionalHeterogeneity.NarrativeReport.txt");
	}
	else {
		print "To run a Compositional Heterogeneity analysis in BOSTIn, you need to specify whether you are using Amino Acid or Nucleotide data. For amino acid data, use:
		
		BOSTIn.pl --ch protein <AlignmentFile> <Prefix for Output files>
		
		For nucleotide data:
		
		BOSTIn.pl --ch dna <AlignmentFile> <Prefix for Output files>\n";
	}
}

sub do_branchLength {
	if ($ARGV[0]=~ /dna/){
		print "LB-Score DNA!\n";
		system("CombinedNJ_LBScore.R $ARGV[1] $ARGV[2].tre $ARGV[2].LBi-scores $ARGV[2].LBSummary.txt DNA");
		system("LBNarrativeReport.pl $ARGV[2] > $ARGV[2].BranchHeterogeneity.NarrativeReport.txt");
	}
	elsif ($ARGV[0]=~ /protein/){
		print "LB-Score AA!\n";
		system("CombinedNJ_LBScore.R $ARGV[1] $ARGV[2].tre $ARGV[2].LBi-scores $ARGV[2].LBSummary.txt AA");
		system("LBNarrativeReport.pl $ARGV[2] > $ARGV[2].BranchHeterogeneity.NarrativeReport.txt");
	}
	else {
		print "To run a Branch Length analysis in BOSTIn, you need to specify whether you are using Amino Acid or Nucleotide data. For amino acid data, use:
		BOSTIn.pl --blh protein <AlignmentFile> <Prefix for Output files>
		
		For nucleotide data:
		
		BOSTIn.pl --blh dna <AlignmentFile> <Prefix for Output files>\n";
	}
}

sub do_siteSat {
	if ($ARGV[0]=~ /dna/){
		print "Site Saturation DNA!\n";
		system ("CValueCalculator.pl $ARGV[1] $ARGV[2]");
		system ("CFactorAnalysis_Script.r $ARGV[2].SiteSaturation.TaxaCScore.txt");
		system ("CValue.NarrativeReport.pl $ARGV[2] > $ARGV[2].SiteSaturation.NarrativeReport.txt");	
		}
	elsif ($ARGV[0]=~ /protein/){
		print "Site Saturation AA!\n";
		system ("DEScoreCalculator.pl $ARGV[1] $ARGV[2]");
		system ("DEScoreAnalysis_Script.r $ARGV[2].SiteSaturation.TaxaFrequencies.txt");
		system ("DEScore.NarrativeReport.pl $ARGV[2] > $ARGV[2].SiteSaturation.NarrativeReport.txt");
	}
	else {
		print "To run a site saturation analysis in BOSTIn, you need to specify whether you are using Amino Acid or Nucleotide data. For amino acid data, use:
		
		BOSTIn.pl --s protein <AlignmentFile> <Prefix for Output files>
		
		For nucleotide data:
		
		BOSTIn.pl --s dna <AlignmentFile> <Prefix for Output files>\n";
	}
}

sub do_branchLengthSiteSat{
	print "combination Site Sat Branch Length analysis for DNA to save NJ Tree construction time has not yet been implemented.";
}
