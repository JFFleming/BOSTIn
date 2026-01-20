#!/usr/bin/perl
use strict;
use warnings;

my $version   = "1.6";

print "
DE-Score Calculator Version $version, James F. Fleming and Torsten H. Struck, 2025.
Welcome to the DE-Score Calculator. DE-Score Calculator accepts amino acid FASTA files as input, and then outputs 2 files.
To run DE-Score Calculator on your data, use the following command:
perl DEScoreCalculator.pl <filename> <prefix for output files>

The 2 output files are:
- <prefix for output files>.SiteSaturation.TaxaFrequencies.txt: This file gives the Within Category/Between Category exchange frequency, the standard deviation of that frequency and the taxon-specific DE-Score for each taxa.
- <prefix for output files>.SiteSaturation.TotalFrequencies.txt: This file gives the Within Category/Between Category exchange frequency, the standard deviation of that frequency and the DE-Score for the whole dataset.

If you have any questions, queries or comments, please don't hesitate to get in touch at:
j.fleming\@ub.edu
";

# Check which input files are provided, and trigger the automatic input if necessary.
my ($align_check, $output_check);
if (@ARGV < 2 ) {
print "Only one or two arguments were provided, sorry. Please note the above command to run DE-Score Calculator. You may have forgotten to specify the alignment file or the output prefix.\n";
exit(1);
} 
else {
    $align_check  = $ARGV[0];
    $output_check = $ARGV[1];
}
my %phy_seqs = parse_fasta($align_check);

my $freq_file = "$output_check.SiteSaturation.TotalFrequencies.txt";
my $tax_file  = "$output_check.SiteSaturation.TaxaFrequencies.txt";

# Create the critical DE-Score
my $size = keys %phy_seqs;
my $normalisation_constant = 0.255 * $size**-0.15;
my $crit_DaCER             = 0.265;
my $saturation_DaCER       = 0.177;
my $crit = ($crit_DaCER - $saturation_DaCER) / $normalisation_constant;
my $round_crit = sprintf("%.5f", $crit);

open(FREQ, '>', $freq_file) or die "Cannot open $freq_file\n";
print FREQ "Summary of the DE-Score Analysis\n";

open(TAXA, '>', $tax_file) or die "Cannot open $tax_file\n";
print TAXA "FileName\tExchangeFreq\tExchangeFreqStDev\tDE-Score\n";

#Set up the 6 Dayhoff Categories
my @small        = qw(A G P S T);
my @acid_amide   = qw(D E N Q);
my @basic        = qw(H K R);
my @hydrophobic  = qw(I L V M);
my @aromatic     = qw(F W Y);
my @sulfur       = qw(C);
my @valid_seqs   = qw(A G P S T D E N Q H K R I L V M F W Y C);

# Create hashes for lookup
my %valids = map { $_ => 1 } @valid_seqs;

my %category;
$category{$_} = 'small'        for @small;
$category{$_} = 'acid_amide'   for @acid_amide;
$category{$_} = 'basic'        for @basic;
$category{$_} = 'hydrophobic'  for @hydrophobic;
$category{$_} = 'aromatic'     for @aromatic;
$category{$_} = 'sulfur'       for @sulfur;

print "There are $size taxa in this input dataset\n";
print "The Critical DE-Score for this dataset is thereby $round_crit \n";
print TAXA "CRITICAL\t0.266\tN/A\t$round_crit\n";

my $a = 0;
my @all_tvs;
my @all_tis;
my @all_ti_freqs;

my @taxa_ids = sort keys %phy_seqs;
my %taxon_tis;
my %taxon_tvs;
my %taxon_ti_freqs;

# Loop avoiding redundant comparisons (i < j only)
for (my $i = 0; $i < @taxa_ids; $i++) {
    my $id1   = $taxa_ids[$i];
    my @vals1 = @{ $phy_seqs{$id1} };
    my $tracker = $i+1;
	print "Comparing $taxa_ids[$i] Taxon $tracker of $size \n";
    for (my $j = $i + 1; $j < @taxa_ids; $j++) {
        my $id2   = $taxa_ids[$j];
        my @vals2 = @{ $phy_seqs{$id2} };

        my ($pairwise_ti, $pairwise_tv, $pos) = (0, 0, 0);

        for my $aa1 (@vals1) {
            my $aa2 = $vals2[$pos];

            if ($valids{$aa1} && $valids{$aa2}) {
                if ($aa1 eq "-" || $aa2 eq "-") {
                    # skip gaps
                } elsif ($aa1 eq $aa2) {
                    # identical, skip
                } elsif ($category{$aa1} && $category{$aa2} && $category{$aa1} eq $category{$aa2}) {
                    $pairwise_ti++;
                } else {
                    $pairwise_tv++;
                }
            }
            $pos++;
        }

        next if ($pairwise_ti == 0 && $pairwise_tv == 0);

        my $ti_freq = ($pairwise_ti == 0) ? 0 : $pairwise_ti / ($pairwise_ti + $pairwise_tv);

        push @all_tis, $pairwise_ti;
        push @all_tvs, $pairwise_tv;
        push @all_ti_freqs, $ti_freq;

        push @{ $taxon_tis{$id1} }, $pairwise_ti;
        push @{ $taxon_tvs{$id1} }, $pairwise_tv;
        push @{ $taxon_ti_freqs{$id1} }, $ti_freq;

        push @{ $taxon_tis{$id2} }, $pairwise_ti;
        push @{ $taxon_tvs{$id2} }, $pairwise_tv;
        push @{ $taxon_ti_freqs{$id2} }, $ti_freq;
    }
}

foreach my $taxon_id (@taxa_ids) {
    $a++;
    my $tis_ref   = $taxon_tis{$taxon_id}     || [];
    my $tvs_ref   = $taxon_tvs{$taxon_id}     || [];
    my $freq_ref  = $taxon_ti_freqs{$taxon_id}|| [];

    my $avg_freq  = avg($freq_ref);
    my $std_freq  = get_stddev($freq_ref);
    my $tax_dist  = $avg_freq - $saturation_DaCER;
    my $tax_DE    = $tax_dist / $normalisation_constant;
	my $round_avg_freq = sprintf("%.5f", $avg_freq);
	my $round_std_freq = sprintf("%.5f", $std_freq);
	my $round_tax_DE = sprintf("%.5f", $tax_DE);
#    print "\n Finished assessing $taxon_id. Taxon $a of $size";
    print TAXA "$taxon_id\t$round_avg_freq\t$round_std_freq\t$round_tax_DE\n";
}

print "\n";

# Calculate the whole dataset DE-Score
my $average_ti_freq = avg(\@all_ti_freqs);
my $std_ti_freq     = get_stddev(\@all_ti_freqs);
my $freq_dist       = $average_ti_freq - $saturation_DaCER;
my $total_DE        = $freq_dist / $normalisation_constant;
my $round_average_ti_freq = sprintf("%.5f", $average_ti_freq);
my $round_std_ti_freq     = sprintf("%.5f", $std_ti_freq);
my $round_total_DE = sprintf("%.5f", $total_DE);

print FREQ "Dayhoff Category Exchange Frequency:\t$round_average_ti_freq\nExchange Frequency Standard Deviation:\t$round_std_ti_freq\nCritical DE-Score:\t$round_crit\nDE-Score:\t$round_total_DE\n";

### Math functions ###
sub avg {
    my ($arr_ref) = @_;
    my $count = @$arr_ref;
    return 0 if $count == 0;
    my $sum = 0;
    $sum += $_ for @$arr_ref;
    return $sum / $count;
}

sub get_stddev {
    return sqrt(get_disp(@_));
}

sub get_disp {
    my ($arr_ref) = @_;
    my $mean  = avg($arr_ref);
    my $count = @$arr_ref;
    return 0 if $count == 0;
    my $sum = 0;
    $sum += ($_ - $mean) ** 2 for @$arr_ref;
    return $sum / $count;
}

### FASTA Parser ###
sub parse_fasta {
    my ($file) = @_;
    my %sequences;
    open(my $fh, "<", $file) or die "Cannot open FASTA file $file\n";

    my ($id, $seq) = ("", "");

    while (my $line = <$fh>) {
        chomp($line);
        if ($line =~ /^>(\S+)/) {
            if ($id ne "") {
                $sequences{$id} = [split //, $seq];
            }
            $id = $1;
            $seq = "";
        } else {
            $seq .= $line;
        }
    }

    if ($id ne "") {
        $sequences{$id} = [split //, $seq];
    }

    close($fh);
    return %sequences;
}
