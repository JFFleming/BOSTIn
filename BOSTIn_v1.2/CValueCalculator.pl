#!/usr/bin/perl
use strict;
use warnings;

#Read the Fasta file
my $align_check = $ARGV[0];
my %phy_seqs = parse_fasta($align_check);
my $freq_file = "$ARGV[1].SiteSaturation.TotalCScore.txt";
my $tax_file = "$ARGV[1].SiteSaturation.TaxaCScore.txt";
#Open output files for writing.
open (FREQ, '>', $freq_file) || die ("Can not open $freq_file\n");

open (TAXA, '>', $tax_file) || die ("Can not open $tax_file\n");
print TAXA "FileName\tTiTv\tTaxonPValue\tCScore\n";
#Set purines and Pyrimidines to measure TiTv

my @pyramidine = qw (T C);
my @purine= qw (A G);
my @valid_seqs = qw (A T C G);

my %valids = map { $_ => 1 } @valid_seqs;

my %category;
$category{$_} = 'pyramidine'        for @pyramidine;
$category{$_} = 'purine'   for @purine;

my $a = 0;
my @all_tvs;
my @all_tis;
my @all_ti_tv;
my @all_same;
my @all_p_dist;

my @taxa_ids = sort keys %phy_seqs;
my $size = keys %phy_seqs;
my %taxon_tis;
my %taxon_tvs;
my %taxon_ti_tv;
my %taxon_same;
my %taxon_p_dist;

# Loop avoiding redundant comparisons (i < j only)
for (my $i = 0; $i < @taxa_ids; $i++) {
    my $id1   = $taxa_ids[$i];
    my @vals1 = @{ $phy_seqs{$id1} };
    my $tracker = $i+1;
	print "Comparing $taxa_ids[$i] Taxon $tracker of $size \n";
    for (my $j = $i + 1; $j < @taxa_ids; $j++) {
        my $id2   = $taxa_ids[$j];
        my @vals2 = @{ $phy_seqs{$id2} };

        my ($pairwise_ti, $pairwise_tv, $same, $pos) = (0, 0, 0, 0);

        for my $nucl1 (@vals1) {
            my $nucl2 = $vals2[$pos];

            if ($valids{$nucl1} && $valids{$nucl2}) {
                if ($nucl1 eq "-" || $nucl2 eq "-") {
                    # skip gaps
                } elsif ($nucl1 eq $nucl2) {
                    # identical, skip
                    $same++;
                } elsif ($category{$nucl1} && $category{$nucl2} && $category{$nucl1} eq $category{$nucl2}) {
                    $pairwise_ti++;
                } else {
                    $pairwise_tv++;
                }
            }
            $pos++;
        }

        next if ($pairwise_ti == 0 && $pairwise_tv == 0);

        my $ti_tv = ($pairwise_ti == 0) ? 0 : $pairwise_ti / $pairwise_tv;
        print $ti_tv, "\n";
		my $p_dist = ($pairwise_ti + $pairwise_tv == 0) ? 0 : ($pairwise_ti + $pairwise_tv)/($pairwise_ti+$pairwise_tv+$same);
		print $p_dist, "\n";

        push @all_tis, $pairwise_ti;
        push @all_tvs, $pairwise_tv;
        push @all_ti_tv, $ti_tv;
        push @all_same, $same;
        push @all_p_dist, $p_dist;

        push @{ $taxon_tis{$id1} }, $pairwise_ti;
        push @{ $taxon_tvs{$id1} }, $pairwise_tv;
        push @{ $taxon_ti_tv{$id1} }, $ti_tv;
		push @{	$taxon_same{$id1} }, $same;
		push @{ $taxon_p_dist{$id1} }, $p_dist; 	
		
        push @{ $taxon_tis{$id2} }, $pairwise_ti;
        push @{ $taxon_tvs{$id2} }, $pairwise_tv;
        push @{ $taxon_ti_tv{$id2} }, $ti_tv;
        push @{	$taxon_same{$id2} }, $same;
        push @{ $taxon_p_dist{$id2} }, $p_dist;
    }
}

foreach my $taxon_id (@taxa_ids) {
    $a++;
    my $ti_ref   = $taxon_tis{$taxon_id}     || [];
    my $tv_ref   = $taxon_tvs{$taxon_id}     || [];
    my $ti_tv_ref  = $taxon_ti_tv{$taxon_id}|| [];
    my $same_ref =  $taxon_same{$taxon_id}|| [];
	my $p_dist_ref = $taxon_p_dist{$taxon_id} || [];

    my $avg_freq  = avg($ti_tv_ref);
    my $std_ti_tv  = get_stddev($ti_tv_ref);
    print $std_ti_tv, "\n";
    my $tax_p_dist = get_stddev($p_dist_ref);
    print $tax_p_dist, "\n";
    my $tax_c_score = ($std_ti_tv == 0) ? 0 :$std_ti_tv/$tax_p_dist;
	
#    print "\n Finished assessing $taxon_id. Taxon $a of $size";
	print TAXA "$taxon_id\t$std_ti_tv\t$tax_p_dist\t$tax_c_score\n";
}

my $average_ti_tv_freq = avg(\@all_ti_tv);
my $std_ti_tv_freq     = get_stddev(\@all_ti_tv);
my $p_value = get_stddev(\@all_p_dist);
my $c_score = ($std_ti_tv_freq == 0) ? 0 : $std_ti_tv_freq/$p_value;

print FREQ "$ARGV[0] C-Score Summary\nStdDevTiTv\t$std_ti_tv_freq\nPValue\t$p_value\nCScore\t$c_score\n";

sub avg {
  my ($avg_array_ref) = @_;
  my $avg_count = @$avg_array_ref;
  my $avg_sum = 0;
  for my $avg_num (@$avg_array_ref) {
      $avg_sum += $avg_num;
  }
  return $avg_sum / $avg_count;
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
  return $dis_sum / $dis_count;
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
