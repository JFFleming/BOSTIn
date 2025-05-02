#!/usr/bin/perl
use strict;
use warnings;
use FAST::Bio::SeqIO;
my $align_check = $ARGV[0] or die "No Alignment File Provided \n";
my $output_check = $ARGV[1] or die "No Output Prefix provided \n";
my $align_file = FAST::Bio::SeqIO->new(-file => $ARGV[0], -format => 'Fasta', -alphabet => 'protein');
my $freq_file = "$ARGV[1].SiteSaturation.TotalFrequencies.txt";
my $tax_file = "$ARGV[1].SiteSaturation.TaxaFrequencies.txt";

open (FREQ, '>', $freq_file) || die ("Can not open $freq_file\n");
print FREQ "FileName\tExchangeFreq\tExchangeFreqStDev\tDE-Score\n";
open (TAXA, '>', $tax_file) || die ("Can not open $tax_file\n");
print TAXA "FileName\tExchangeFreq\tExchangeFreqStDev\tDE-Score\n";
my @all_seqs = ();
my @small= ("A","G","P","S","T");
my @acid_amide= ("D","E","N","Q");
my @basic= ("H","K","R");
my @hydrophobic= ("I","L","V","M");
my @aromatic= ("F","W","Y");
my @sulfur= ("C");

# Create hashes for fast lookup of each category
my %small = map { $_ => 1 } @small;
my %acid_amide = map { $_ => 1 } @acid_amide;
my %basic = map { $_ => 1 } @basic;
my %hydrophobic = map { $_ => 1 } @hydrophobic;
my %aromatic = map { $_ => 1 } @aromatic;
my %sulfur = map { $_ => 1 } @sulfur;

my %phy_seqs;

print "
DE-Score Calculator, James F. Fleming and Torsten H. Struck, 2024.
Welcome to the DE-Score Calculator. DE-Score Calculator accepts amino acid FASTA files as input, and then outputs 2 files.
To run DE-Score Calculator on your data, use the following command:
perl DEScoreCalculator.pl <filename> <prefix for output files>

The 2 output files are:
- <prefix for output files>.SiteSaturation.TaxaFrequencies.txt: This file gives the Within Category/Between Category exchange frequency, the standard deviation of that frequency and the taxon-specific DE-Score for each taxa.
- <prefix for output files>.SiteSaturation.TotalFrequencies.txt: This file gives the Within Category/Between Category exchange frequency, the standard deviation of that frequency and the DE-Score for the whole dataset.

If you have any questions, queries or comments, please don't hesitate to get in touch at:
jfleming\@jamstec.go.jp
";
if($ARGV[0]=""||$ARGV[1]=""){
    print "Are you sure you formatted your input correctly? Check the description above to be sure.";
}

### This reads the Fasta file in ###
while ( my $seq = $align_file->next_seq() ) {
    my $id = $seq->display_id;
    my $subseq = $seq->subseq(1, $seq->length());
    my @sequence = split(//,$subseq);
    $phy_seqs{$id}=[@sequence];
}

my $size = keys %phy_seqs;
print "There are $size taxa in this input dataset\n";
my $a = 0;
my @all_tvs;
my @all_tis;
my @all_ti_freqs;

### This loop counts pairwise Dayhoff Exchanges ###
foreach my $k (keys %phy_seqs){
    my @vals = @{ $phy_seqs{$k}};
    my $b = 0;
    my @tax_tvs;
    my @tax_tis;
    my @tax_ti_freqs;
    
    foreach my $compk (keys %phy_seqs){
        next if ($k eq $compk);
        my @comp_vals = @{ $phy_seqs{$compk}};
        my $i = 0;
        my $pairwise_ti = 0;
        my $pairwise_tv = 0;
        my $same_counter = 0;
        my $small_counter = 0;
        my $gap_counter = 0;
        
        foreach (@vals){
            my $query = $_;
            my $comp_query = $comp_vals[$i];
            
            if ($query eq "-" || $comp_query eq "-"){
                $gap_counter++;
            }
            elsif ($comp_query eq $query){
                $same_counter++;
            }
            # Check for amino acids in the same category
            elsif ($small{$query} && $small{$comp_query}) {
                $pairwise_ti++;
            }
            elsif ($acid_amide{$query} && $acid_amide{$comp_query}) {
                $pairwise_ti++;
            }
            elsif ($basic{$query} && $basic{$comp_query}) {
                $pairwise_ti++;
            }
            elsif ($hydrophobic{$query} && $hydrophobic{$comp_query}) {
                $pairwise_ti++;
            }
            elsif ($aromatic{$query} && $aromatic{$comp_query}) {
                $pairwise_ti++;
            }
            elsif ($sulfur{$query} && $sulfur{$comp_query}) {
                $pairwise_ti++;
            }
            else{
                $pairwise_tv++;
            }
            $i++;
        }
        $b++;
        if ($pairwise_ti == 0 && $pairwise_tv == 0) {
            $gap_counter++;
        }
### This part pushes each pairwise TI/TV to the taxon specific and whole dataset arrays. ###
        else{
	        my $ti_freq = $pairwise_ti/($pairwise_ti+$pairwise_tv);
	        push (@all_tis, $pairwise_ti);
	        push (@all_tvs, $pairwise_tv);
	        push (@all_ti_freqs, $ti_freq);
	        push (@tax_tis, $pairwise_ti);
	        push (@tax_tvs, $pairwise_tv);
            push (@tax_ti_freqs, $ti_freq);
        }
    }
### Here we do the taxon DE-Score calculation ###
    $a++;
    my $tax_average_transi = avg(\@tax_tis);
    my $tax_average_transv = avg(\@tax_tvs);
    my $tax_average_ti_freq = avg(\@tax_ti_freqs);
    my $tax_std_transi = get_stddev(\@tax_tis);
    my $tax_std_transv = get_stddev(\@tax_tvs);
    my $tax_std_ti_freq =  get_stddev(\@tax_ti_freqs);
    my $tax_dist = $tax_average_ti_freq - 0.177;
    my $tax_DE = $tax_dist/(0.255*$size**-0.15);
    print "\n Finished assessing $k. Taxon $a of $size";
    print TAXA "$k\t$tax_average_ti_freq\t$tax_std_ti_freq\t$tax_DE\n";
}
print "\n";
### Here we do the whole dataset DE-Score calculation ###
my $average_transi = avg(\@all_tis);
my $average_transv = avg(\@all_tvs);
my $average_ti_freq = avg(\@all_ti_freqs);
my $std_transi = get_stddev(\@all_tis);
my $std_transv = get_stddev(\@all_tvs);
my $std_ti_freq =  get_stddev(\@all_ti_freqs);
my $freq_dist = $average_ti_freq - 0.177;
my $total_DE = $freq_dist/(0.255*$size**-0.15);

print FREQ "$ARGV[0]\t$average_ti_freq\t$std_ti_freq\t$total_DE\n";

### These are the mathematics subroutines for averages, standard deviations and disparity ###
sub avg {
  my ($avg_array_ref) = @_;
  my $avg_count = @$avg_array_ref;
  my $avg_sum = 0;
  for my $avg_num (@$avg_array_ref) {
      $avg_sum += $avg_num;
  }
  if ($avg_sum == 0){
      return $avg_sum;
  }
  else{
      return $avg_sum / $avg_count;
  }
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
  if ($dis_sum == 0){
      return $dis_sum;
  }
  else{
      return $dis_sum / $dis_count;
  }
}
