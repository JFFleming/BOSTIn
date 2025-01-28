#!/usr/bin/perl
use 5.014;
#use warnings;
use FAST::Bio::SeqIO;

print "
RCFV Reader, James F. Fleming & Torsten H Struck, 2022.
Welcome to RCFV Reader. RCFV Reader accepts FASTA files as input, and then outputs 6 files.
To run RCFV Reader on nucleotide data, input the command as follows:

perl RCFVReader.pl dna <filename> <prefix for output files>

for amino acid data, input the command like this:

perl RCFVReader.pl protein <filename> <prefix for output files>

the 6 output files are:
- RCFV.txt - this file contains the total RCFV and nRCFV of the entire dataset.
- csRCFV.txt - this file contains a list of the character-specific RCFVs of your input dataset
- ncsRCFV.txt - as per csRCFV.txt, but with ncsRCFV values
- tRCFV.txt - this file contains a list of the taxon-specific RCFVs of your input dataset
- ntRCFV.txt - as per tRCFV.txt, but with ntRCFV values
- Frequencies.txt - this file contains the per taxon character frequencies, and the mean per taxon character frequencies across the dataset (as used by RCFV, NOT the total dataset character frequencies)

I hope you enjoy this fast new way to use normalised Relative Frequency Composition Values! If you have any questions or queries, or notice any bugs, contact me at:
j.f.fleming\@nhm.uio.no
";
#This if statement checks that the alphabet variable is valid
if ($ARGV[0]!=~/dna/||/protein/){
	print "WARNING: Are you sure you remembered to specify protein or dna for amino acid or nucleotide data?\n";}
#Set up the input and output files
my $fasta  = FAST::Bio::SeqIO->new(-file => $ARGV[1], -format => 'Fasta', -alphabet => $ARGV[0]);
#print "my $fasta  = FAST::Bio::SeqIO->new(-file => $ARGV[1], -format => 'Fasta', -alphabet => $ARGV[0]);";
my $seqnum=0;
my $RCFV_File = "$ARGV[2].RCFV.txt";
#my $csRCFV_File = "$ARGV[2].csRCFV.txt";
#my $tRCFV_File = "$ARGV[2].tRCFV.txt";
my $ncsRCFV_File = "$ARGV[2].ncsRCFV.txt";
my $ntRCFV_File = "$ARGV[2].ntRCFV.txt";
my $outfile = "$ARGV[2].Frequencies.txt";
open(OUT, '>', $outfile) or die $!;
open(RCFV_OUT, '>', $RCFV_File) or die $!;
#open(CSRCFV_OUT, '>', $csRCFV_File) or die $!;
#open(TRCFV_OUT, '>', $tRCFV_File) or die $!;
open(NCSRCFV_OUT,  '>', $ncsRCFV_File) or die $!;
open(NTRCFV_OUT, '>', $ntRCFV_File) or die $!;

#Nucleotide RCFV Calculator
if ($ARGV[0]=~ /dna/){
#print "boop";
my %A_Lens;
my %C_Lens;
my %G_Lens;
my %T_Lens;
my %A_Freqs;
my %C_Freqs;
my %G_Freqs;
my %T_Freqs;
my @ID_List;
my @Length;

#print OUT "NAME\tFreq(A)\tFreq(C)\tFreq(G)\tFreq(T)\n";
#Calculate frequencies and lengths
while ( my $seq = $fasta->next_seq() ) {
    my $stats;
    $stats->{len} = length($seq->seq);
    push (@Length, length($seq->seq));
    $stats->{$_}++ for split //, $seq->seq;
#    say ++$seqnum, " @$stats{qw(len A C G T a c g t)}";
    my $seqTotal = @$stats{qw(A)} + @$stats{qw(C)} + @$stats{qw(T)} + @$stats{qw(G)} +@$stats{qw(a)} +@$stats{qw(c)} +@$stats{qw(t)} +@$stats{qw(g)};
    my $seqA_Freq = (@$stats{qw(A)} + @$stats{qw(a)})/$seqTotal;
    my $seqC_Freq = (@$stats{qw(C)} + @$stats{qw(c)})/$seqTotal;
    my $seqG_Freq = (@$stats{qw(G)} + @$stats{qw(g)})/$seqTotal;
    my $seqT_Freq = (@$stats{qw(T)} + @$stats{qw(t)})/$seqTotal;
#	print OUT $seq->id, "\t", $seqA_Freq, "\t", $seqC_Freq, "\t", $seqG_Freq, "\t", $seqT_Freq, "\n";
	push (@ID_List, ($seq->id));
    $A_Lens{$seq->id} = @$stats{qw(A)}+@$stats{qw(a)};
	$C_Lens{$seq->id} = @$stats{qw(C)}+@$stats{qw(c)};
	$G_Lens{$seq->id} = @$stats{qw(G)}+@$stats{qw(g)};
	$T_Lens{$seq->id} = @$stats{qw(T)}+@$stats{qw(t)};
    $A_Freqs{$seq->id} = $seqA_Freq;
    $C_Freqs{$seq->id} = $seqC_Freq;
    $G_Freqs{$seq->id} = $seqG_Freq;
    $T_Freqs{$seq->id} = $seqT_Freq;
}

my %check;
@check{@Length} = (1) x @Length;
if (keys %check == 1){
	print "Fasta File Aligned\n";
	} 
	else{
	print "Are you sure this file is aligned? Sequences seem to have differing lengths.\n" and die $!;}


my $A_total = eval join '+', values %A_Freqs;
my $C_total = eval join '+', values %C_Freqs;
my $G_total = eval join '+', values %G_Freqs;
my $T_total = eval join '+', values %T_Freqs;
my $Total = $A_total + $C_total + $G_total + $T_total;

#print $Total, "\t", $A_total, "\t", $C_total, "\t", $G_total, "\t", $T_total, "\n";
#Calculate mean frequencies
my $MeanA_Freq = ($A_total/$Total);
my $MeanC_Freq = ($C_total/$Total);
my $MeanG_Freq = ($G_total/$Total);
my $MeanT_Freq = ($T_total/$Total);

print OUT "Mean_Freq_Across_Taxa\t", $MeanA_Freq, "\t", $MeanC_Freq, "\t", $MeanG_Freq, "\t", $MeanT_Freq, "\n";

my %A_RCFV_List;
my %C_RCFV_List;
my %G_RCFV_List;
my %T_RCFV_List;
#Calculate csRCFV
foreach my $Akey (keys %A_Freqs){
	my $Mu =  abs($A_Freqs{$Akey}-$MeanA_Freq);
	$A_RCFV_List{$Akey} = $Mu;
#	print $Mu, "\n";
	}
my $A_Size = keys %A_RCFV_List;
my $A_RCFV_Total = eval join '+', values %A_RCFV_List;
my $A_RCFV = $A_RCFV_Total/$A_Size;

foreach my $Ckey (keys %C_Freqs){
	my $Mu =  abs($C_Freqs{$Ckey}-$MeanC_Freq);
	$C_RCFV_List{$Ckey} = $Mu;
#	print $Mu, "\n";
	}
my $C_RCFV_Total = eval join '+', values %C_RCFV_List;
my $C_RCFV = $C_RCFV_Total/$A_Size;

foreach my $Gkey (keys %G_Freqs){
	my $Mu =  abs($G_Freqs{$Gkey}-$MeanG_Freq);
	$G_RCFV_List{$Gkey} = $Mu;
#	print $Mu, "\n";
	}
my $G_RCFV_Total = eval join '+', values %G_RCFV_List;
my $G_RCFV = $G_RCFV_Total/$A_Size;

foreach my $Tkey (keys %T_Freqs){
	my $Mu =  abs($T_Freqs{$Tkey}-$MeanT_Freq);
	$T_RCFV_List{$Tkey} = $Mu;
#	print $Mu, "\n";
	}
my $T_RCFV_Total = eval join '+', values %T_RCFV_List;
my $T_RCFV = $T_RCFV_Total/$A_Size;

#print TRCFV_OUT "\ntRCFV values:\n";
print NTRCFV_OUT "\nntRCFV values:\n";

#Calculate tRCFV
my $i= 0;
my $tRCFV;
foreach(@ID_List){
	my $unlock=$_;
	$tRCFV = ($A_RCFV_List{$unlock}+$C_RCFV_List{$unlock}+$G_RCFV_List{$unlock}+$T_RCFV_List{$unlock})/$A_Size;
	#print TRCFV_OUT $_,"\t",$tRCFV,"\n";
	print NTRCFV_OUT $_,"\t",(0.25*$tRCFV*$A_Size*sqrt($Length[0]))/100,"\n";
	$i++;
	}

my $RCFV = $A_RCFV + $C_RCFV + $G_RCFV + $T_RCFV;
my $nRCFV = $RCFV/(($Length[0]**-0.5)*400);

#print CSRCFV_OUT "character RCFV values:\ncsRCFV(A)\t", $A_RCFV, "\ncsRCFV(C)\t", $C_RCFV, "\ncsRCFV(G)\t", $G_RCFV, "\ncsRCFV(T)\t", $T_RCFV, "\n";

print NCSRCFV_OUT "character RCFV values:\ncsRCFV(A)\t", $A_RCFV/(($Length[0]**-0.5)*100), "\ncsRCFV(C)\t", $C_RCFV/(($Length[0]**-0.5)*100), "\ncsRCFV(G)\t", $G_RCFV/(($Length[0]**-0.5)*100), "\ncsRCFV(T)\t", 
$T_RCFV/(($Length[0]**-0.5)*100), "\n";

print RCFV_OUT "total RCFV\nRCFV\t", $RCFV, "\nnRCFV\t", $nRCFV, "\n";
}

#AA data
elsif ($ARGV[0]=~/protein/){
my %A_Len;
my %V_Len;
my %L_Len;
my %I_Len;
my %P_Len;
my %M_Len;
my %F_Len;
my %W_Len;
my %G_Len;
my %S_Len;
my %T_Len;
my %C_Len;
my %N_Len;
my %Q_Len;
my %Y_Len;
my %D_Len;
my %E_Len;
my %K_Len;
my %R_Len;
my %H_Len;

my %A_Freqs;
my %V_Freqs;
my %L_Freqs;
my %I_Freqs;
my %P_Freqs;
my %M_Freqs;
my %F_Freqs;
my %W_Freqs;
my %G_Freqs;
my %S_Freqs;
my %T_Freqs;
my %C_Freqs;
my %N_Freqs;
my %Q_Freqs;
my %Y_Freqs;
my %D_Freqs;
my %E_Freqs;
my %K_Freqs;
my %R_Freqs;
my %H_Freqs;

my @ID_List;
my @Length;

print OUT "NAME\tFreq(A)\tFreq(V)\tFreq(L)\tFreq(I)\tFreq(P)\tFreq(M)\tFreq(F)\tFreq(W)\tFreq(G)\tFreq(S)\tFreq(T)\tFreq(C)\tFreq(N)\tFreq(Q)\tFreq(Y)\tFreq(D)\tFreq(E)\tFreq(K)\tFreq(R)\tFreq(H)\n";
while ( my $seq = $fasta->next_seq() ) {
    my $stats;
    $stats->{len} = length($seq->seq);
    push (@Length, length($seq->seq));
    $stats->{$_}++ for split //, $seq->seq;
#   say ++$seqnum, " @$stats{qw(len A V L I P M F W G S T C N Q Y D E K R H)}";
#Calculate frequencies and lengths
    my $seqTotal = @$stats{qw(A)} + @$stats{qw(V)} + @$stats{qw(L)} + @$stats{qw(I)} + @$stats{qw(P)} + @$stats{qw(M)} + @$stats{qw(F)} + @$stats{qw(W)} + @$stats{qw(G)} + @$stats{qw(S)} + @$stats{qw(T)} + 
@$stats{qw(C)} + @$stats{qw(N)} + @$stats{qw(Q)} + @$stats{qw(Y)} + @$stats{qw(D)} + @$stats{qw(E)} + @$stats{qw(K)} + @$stats{qw(R)} + @$stats{qw(H)};
    my $seqA_Freq = @$stats{qw(A)}/$seqTotal;
    my $seqV_Freq = @$stats{qw(V)}/$seqTotal;
    my $seqL_Freq = @$stats{qw(L)}/$seqTotal;
    my $seqI_Freq = @$stats{qw(I)}/$seqTotal;
    my $seqP_Freq = @$stats{qw(P)}/$seqTotal;
    my $seqM_Freq = @$stats{qw(M)}/$seqTotal;
    my $seqF_Freq = @$stats{qw(F)}/$seqTotal;
    my $seqW_Freq = @$stats{qw(W)}/$seqTotal;
    my $seqG_Freq = @$stats{qw(G)}/$seqTotal;
    my $seqS_Freq = @$stats{qw(S)}/$seqTotal;
    my $seqT_Freq = @$stats{qw(T)}/$seqTotal;
    my $seqC_Freq = @$stats{qw(C)}/$seqTotal;
    my $seqN_Freq = @$stats{qw(N)}/$seqTotal;
    my $seqQ_Freq = @$stats{qw(Q)}/$seqTotal;
    my $seqY_Freq = @$stats{qw(Y)}/$seqTotal;
    my $seqD_Freq = @$stats{qw(D)}/$seqTotal;
    my $seqE_Freq = @$stats{qw(E)}/$seqTotal;
    my $seqK_Freq = @$stats{qw(K)}/$seqTotal;
    my $seqR_Freq = @$stats{qw(R)}/$seqTotal;
    my $seqH_Freq = @$stats{qw(H)}/$seqTotal;
    
	print OUT $seq->id, "\t", $seqA_Freq ,"\t",    $seqV_Freq ,"\t",    $seqL_Freq ,"\t",    $seqI_Freq ,"\t",    $seqP_Freq ,"\t",    $seqM_Freq ,"\t",    $seqF_Freq ,"\t",    $seqW_Freq ,"\t",    $seqG_Freq 
,"\t",    $seqS_Freq ,"\t",    $seqT_Freq ,"\t",    $seqC_Freq ,"\t",    $seqN_Freq ,"\t",    $seqQ_Freq ,"\t",    $seqY_Freq ,"\t",    $seqD_Freq ,"\t",    $seqE_Freq ,"\t",    $seqK_Freq ,"\t",    $seqR_Freq 
,"\t",    $seqH_Freq , "\n";
	push (@ID_List, ($seq->id));
    $A_Len{$seq->id} = @$stats{qw(A)};
    $V_Len{$seq->id} = @$stats{qw(V)};
    $L_Len{$seq->id} = @$stats{qw(L)};
    $I_Len{$seq->id} = @$stats{qw(I)};
    $P_Len{$seq->id} = @$stats{qw(P)};
    $M_Len{$seq->id} = @$stats{qw(M)};
    $F_Len{$seq->id} = @$stats{qw(F)};
    $W_Len{$seq->id} = @$stats{qw(W)};
    $G_Len{$seq->id} = @$stats{qw(G)};
    $S_Len{$seq->id} = @$stats{qw(S)};
    $T_Len{$seq->id} = @$stats{qw(T)};
    $C_Len{$seq->id} = @$stats{qw(C)};
    $N_Len{$seq->id} = @$stats{qw(N)};
    $Q_Len{$seq->id} = @$stats{qw(Q)};
    $Y_Len{$seq->id} = @$stats{qw(Y)};
    $D_Len{$seq->id} = @$stats{qw(D)};
    $E_Len{$seq->id} = @$stats{qw(E)};
    $K_Len{$seq->id} = @$stats{qw(K)};
    $R_Len{$seq->id} = @$stats{qw(R)};
    $H_Len{$seq->id} = @$stats{qw(H)};

    $A_Freqs{$seq->id} = $seqA_Freq;
    $V_Freqs{$seq->id} = $seqV_Freq;
    $L_Freqs{$seq->id} = $seqL_Freq;
    $I_Freqs{$seq->id} = $seqI_Freq;
    $P_Freqs{$seq->id} = $seqP_Freq;
    $M_Freqs{$seq->id} = $seqM_Freq;
    $F_Freqs{$seq->id} = $seqF_Freq;
    $W_Freqs{$seq->id} = $seqW_Freq;
    $G_Freqs{$seq->id} = $seqG_Freq;
    $S_Freqs{$seq->id} = $seqS_Freq;
    $T_Freqs{$seq->id} = $seqT_Freq;
    $C_Freqs{$seq->id} = $seqC_Freq;
    $N_Freqs{$seq->id} = $seqN_Freq;
    $Q_Freqs{$seq->id} = $seqQ_Freq;
    $Y_Freqs{$seq->id} = $seqY_Freq;
    $D_Freqs{$seq->id} = $seqD_Freq;
    $E_Freqs{$seq->id} = $seqE_Freq;
    $K_Freqs{$seq->id} = $seqK_Freq;
    $R_Freqs{$seq->id} = $seqR_Freq;
    $H_Freqs{$seq->id} = $seqH_Freq;

}

my %check;
@check{@Length} = (1) x @Length;
if (keys %check == 1){
	print "Fasta File Aligned\n";
	} 
	else{
	print "Are you sure this file is aligned? Sequences seem to have differing lengths.\n" and die $!;}

my $A_total = eval join '+', values %A_Freqs;
my $V_total = eval join '+', values %V_Freqs;
my $L_total = eval join '+', values %L_Freqs;
my $I_total = eval join '+', values %I_Freqs;
my $P_total = eval join '+', values %P_Freqs;
my $M_total = eval join '+', values %M_Freqs;
my $F_total = eval join '+', values %F_Freqs;
my $W_total = eval join '+', values %W_Freqs;
my $G_total = eval join '+', values %G_Freqs;
my $S_total = eval join '+', values %S_Freqs;
my $T_total = eval join '+', values %T_Freqs;
my $C_total = eval join '+', values %C_Freqs;
my $N_total = eval join '+', values %N_Freqs;
my $Q_total = eval join '+', values %Q_Freqs;
my $Y_total = eval join '+', values %Y_Freqs;
my $D_total = eval join '+', values %D_Freqs;
my $E_total = eval join '+', values %E_Freqs;
my $K_total = eval join '+', values %K_Freqs;
my $R_total = eval join '+', values %R_Freqs;
my $H_total = eval join '+', values %H_Freqs;

my $Total = $A_total + $V_total + $L_total + $I_total + $P_total + $M_total + $F_total + $W_total + $G_total + $S_total + $T_total + $C_total + $N_total + $Q_total + $Y_total + $D_total + $E_total + $K_total + 
$R_total + $H_total;


#print $Total, "\t", $A_total, "\t", $V_total ,"\t", $L_total ,"\t", $I_total ,"\t", $P_total ,"\t", $M_total ,"\t", $F_total ,"\t", $W_total ,"\t", $G_total ,"\t", $S_total ,"\t", $T_total ,"\t", $C_total ,"\t", 
$N_total ,"\t", $Q_total ,"\t", $Y_total ,"\t", $D_total ,"\t", $E_total ,"\t", $K_total ,"\t", $R_total ,"\t", $H_total, "\n";
#Calculate mean frequencies
my $MeanA_Freq = ($A_total/$Total);
my $MeanV_Freq = ($V_total/$Total);
my $MeanL_Freq = ($L_total/$Total);
my $MeanI_Freq = ($I_total/$Total);
my $MeanP_Freq = ($P_total/$Total);
my $MeanM_Freq = ($M_total/$Total);
my $MeanF_Freq = ($F_total/$Total);
my $MeanW_Freq = ($W_total/$Total);
my $MeanG_Freq = ($G_total/$Total);
my $MeanS_Freq = ($S_total/$Total);
my $MeanT_Freq = ($T_total/$Total);
my $MeanC_Freq = ($C_total/$Total);
my $MeanN_Freq = ($N_total/$Total);
my $MeanQ_Freq = ($Q_total/$Total);
my $MeanY_Freq = ($Y_total/$Total);
my $MeanD_Freq = ($D_total/$Total);
my $MeanE_Freq = ($E_total/$Total);
my $MeanK_Freq = ($K_total/$Total);
my $MeanR_Freq = ($R_total/$Total);
my $MeanH_Freq = ($H_total/$Total);


print OUT "Mean_Freq_Across_Taxa\t", $MeanA_Freq, "\t", $MeanV_Freq, "\t", $MeanL_Freq, "\t", $MeanI_Freq, "\t",
$MeanP_Freq, "\t", $MeanM_Freq, "\t", $MeanF_Freq, "\t", $MeanW_Freq, "\t",
$MeanG_Freq, "\t", $MeanS_Freq, "\t", $MeanT_Freq, "\t", $MeanC_Freq, "\t",
$MeanN_Freq, "\t", $MeanQ_Freq, "\t", $MeanY_Freq, "\t", $MeanD_Freq, "\t",
$MeanE_Freq, "\t", $MeanK_Freq, "\t", $MeanR_Freq, "\t", $MeanH_Freq, "\t", "\n";


my %A_RCFV_List;
my %V_RCFV_List;
my %L_RCFV_List;
my %I_RCFV_List;
my %P_RCFV_List;
my %M_RCFV_List;
my %F_RCFV_List;
my %W_RCFV_List;
my %G_RCFV_List;
my %S_RCFV_List;
my %T_RCFV_List;
my %C_RCFV_List;
my %N_RCFV_List;
my %Q_RCFV_List;
my %Y_RCFV_List;
my %D_RCFV_List;
my %E_RCFV_List;
my %K_RCFV_List;
my %R_RCFV_List;
my %H_RCFV_List;

#Calculate csRCFV
foreach my $Akey (keys %A_Freqs){
	my $AMu =  abs($A_Freqs{$Akey}-$MeanA_Freq);
	$A_RCFV_List{$Akey} = $AMu;
#	print $Mu, "\n";
	}
my $A_Size = keys %A_RCFV_List;
#print $A_Size;
my $A_RCFV_Total = eval join '+', values %A_RCFV_List;
my $A_RCFV = $A_RCFV_Total/$A_Size;

foreach my $Vkey (keys %V_Freqs){
	my $VMu =  abs($V_Freqs{$Vkey}-$MeanV_Freq);
	$V_RCFV_List{$Vkey} = $VMu;
#	print $Mu, "\n";
	}
my $V_Size = keys %V_RCFV_List;
#print $V_Size;
my $V_RCFV_Total = eval join '+', values %V_RCFV_List;
my $V_RCFV = $V_RCFV_Total/$V_Size;

foreach my $Lkey (keys %L_Freqs){
	my $LMu =  abs($L_Freqs{$Lkey}-$MeanL_Freq);
	$L_RCFV_List{$Lkey} = $LMu;
#	print $Mu, "\n";
	}
my $L_RCFV_Total = eval join '+', values %L_RCFV_List;
my $L_RCFV = $L_RCFV_Total/$A_Size;

foreach my $Ikey (keys %I_Freqs){
	my $IMu =  abs($I_Freqs{$Ikey}-$MeanI_Freq);
	$I_RCFV_List{$Ikey} = $IMu;
#	print $Mu, "\n";
	}
my $I_RCFV_Total = eval join '+', values %I_RCFV_List;
my $I_RCFV = $I_RCFV_Total/$A_Size;

foreach my $Pkey (keys %P_Freqs){
	my $PMu =  abs($P_Freqs{$Pkey}-$MeanP_Freq);
	$P_RCFV_List{$Pkey} = $PMu;
#	print $Mu, "\n";
	}
my $P_RCFV_Total = eval join '+', values %P_RCFV_List;
my $P_RCFV = $P_RCFV_Total/$A_Size;

foreach my $Mkey (keys %M_Freqs){
	my $MMu =  abs($M_Freqs{$Mkey}-$MeanM_Freq);
	$M_RCFV_List{$Mkey} = $MMu;
#	print $Mu, "\n";
	}
my $M_RCFV_Total = eval join '+', values %M_RCFV_List;
my $M_RCFV = $M_RCFV_Total/$A_Size;

foreach my $Fkey (keys %F_Freqs){
	my $FMu =  abs($F_Freqs{$Fkey}-$MeanF_Freq);
	$F_RCFV_List{$Fkey} = $FMu;
#	print $Mu, "\n";
	}
my $F_RCFV_Total = eval join '+', values %F_RCFV_List;
my $F_RCFV = $F_RCFV_Total/$A_Size;

foreach my $Wkey (keys %W_Freqs){
	my $WMu =  abs($W_Freqs{$Wkey}-$MeanW_Freq);
	$W_RCFV_List{$Wkey} = $WMu;
#	print $Mu, "\n";
	}
my $W_RCFV_Total = eval join '+', values %W_RCFV_List;
my $W_RCFV = $W_RCFV_Total/$A_Size;

foreach my $Gkey (keys %G_Freqs){
	my $GMu =  abs($G_Freqs{$Gkey}-$MeanG_Freq);
	$G_RCFV_List{$Gkey} = $GMu;
#	print $Mu, "\n";
	}
my $G_RCFV_Total = eval join '+', values %G_RCFV_List;
my $G_RCFV = $G_RCFV_Total/$A_Size;

foreach my $Skey (keys %S_Freqs){
	my $SMu =  abs($S_Freqs{$Skey}-$MeanS_Freq);
	$S_RCFV_List{$Skey} = $SMu;
#	print $Mu, "\n";
	}
my $S_RCFV_Total = eval join '+', values %S_RCFV_List;
my $S_RCFV = $S_RCFV_Total/$A_Size;

foreach my $Tkey (keys %T_Freqs){
	my $TMu =  abs($T_Freqs{$Tkey}-$MeanT_Freq);
	$T_RCFV_List{$Tkey} = $TMu;
#	print $Mu, "\n";
	}
my $T_RCFV_Total = eval join '+', values %T_RCFV_List;
my $T_RCFV = $T_RCFV_Total/$A_Size;

foreach my $Ckey (keys %C_Freqs){
	my $CMu =  abs($C_Freqs{$Ckey}-$MeanC_Freq);
	$C_RCFV_List{$Ckey} = $CMu;
#	print $Mu, "\n";
	}
my $C_RCFV_Total = eval join '+', values %C_RCFV_List;
my $C_RCFV = $C_RCFV_Total/$A_Size;

foreach my $Nkey (keys %N_Freqs){
	my $NMu =  abs($N_Freqs{$Nkey}-$MeanN_Freq);
	$N_RCFV_List{$Nkey} = $NMu;
#	print $Mu, "\n";
	}
my $N_RCFV_Total = eval join '+', values %N_RCFV_List;
my $N_RCFV = $N_RCFV_Total/$A_Size;

foreach my $Qkey (keys %Q_Freqs){
	my $QMu =  abs($Q_Freqs{$Qkey}-$MeanQ_Freq);
	$Q_RCFV_List{$Qkey} = $QMu;
#	print $Mu, "\n";
	}
my $Q_RCFV_Total = eval join '+', values %Q_RCFV_List;
my $Q_RCFV = $Q_RCFV_Total/$A_Size;

foreach my $Ykey (keys %Y_Freqs){
	my $YMu =  abs($Y_Freqs{$Ykey}-$MeanY_Freq);
	$Y_RCFV_List{$Ykey} = $YMu;
#	print $Mu, "\n";
	}
my $Y_RCFV_Total = eval join '+', values %Y_RCFV_List;
my $Y_RCFV = $Y_RCFV_Total/$A_Size;

foreach my $Dkey (keys %D_Freqs){
	my $DMu =  abs($D_Freqs{$Dkey}-$MeanD_Freq);
	$D_RCFV_List{$Dkey} = $DMu;
#	print $Mu, "\n";
	}
my $D_RCFV_Total = eval join '+', values %D_RCFV_List;
my $D_RCFV = $D_RCFV_Total/$A_Size;

foreach my $Ekey (keys %E_Freqs){
	my $EMu =  abs($E_Freqs{$Ekey}-$MeanE_Freq);
	$E_RCFV_List{$Ekey} = $EMu;
#	print $Mu, "\n";
	}
my $E_RCFV_Total = eval join '+', values %E_RCFV_List;
my $E_RCFV = $E_RCFV_Total/$A_Size;

foreach my $Kkey (keys %K_Freqs){
	my $KMu =  abs($K_Freqs{$Kkey}-$MeanK_Freq);
	$K_RCFV_List{$Kkey} = $KMu;
#	print $Mu, "\n";
	}
my $K_RCFV_Total = eval join '+', values %K_RCFV_List;
my $K_RCFV = $K_RCFV_Total/$A_Size;

foreach my $Rkey (keys %R_Freqs){
	my $RMu =  abs($R_Freqs{$Rkey}-$MeanR_Freq);
	$R_RCFV_List{$Rkey} = $RMu;
#	print $Mu, "\n";
	}
my $R_RCFV_Total = eval join '+', values %R_RCFV_List;
my $R_RCFV = $R_RCFV_Total/$A_Size;

foreach my $Hkey (keys %H_Freqs){
	my $HMu =  abs($H_Freqs{$Hkey}-$MeanH_Freq);
	$H_RCFV_List{$Hkey} = $HMu;
#	print $Mu, "\n";
	}
my $H_RCFV_Total = eval join '+', values %H_RCFV_List;
my $H_RCFV = $H_RCFV_Total/$A_Size;

#print TRCFV_OUT "tRCFV values:\n";
print NTRCFV_OUT "ntRCFV values:\n";

#Calculate tRCFV
my $i= 0;
my $tRCFV;
foreach(@ID_List){
	my $unlock=$_;
	$tRCFV = ($A_RCFV_List{$unlock} + $V_RCFV_List{$unlock} + $L_RCFV_List{$unlock} + $I_RCFV_List{$unlock} + $P_RCFV_List{$unlock} + $M_RCFV_List{$unlock} + $F_RCFV_List{$unlock} + $W_RCFV_List{$unlock} + 
$G_RCFV_List{$unlock} + $S_RCFV_List{$unlock} + $T_RCFV_List{$unlock} + $C_RCFV_List{$unlock} + $N_RCFV_List{$unlock} + $Q_RCFV_List{$unlock} + $Y_RCFV_List{$unlock} + $D_RCFV_List{$unlock} + 
$E_RCFV_List{$unlock} + $K_RCFV_List{$unlock} + $R_RCFV_List{$unlock} + $H_RCFV_List{$unlock})/$A_Size;
	#print TRCFV_OUT $_,"\t",$tRCFV,"\n";
	print NTRCFV_OUT $_,"\t",(0.05*$tRCFV*$A_Size*sqrt($Length[0]))/100,"\n";
	$i++;
	}

my $RCFV = $A_RCFV + $V_RCFV + $L_RCFV + $I_RCFV + $P_RCFV + $M_RCFV + $F_RCFV + $W_RCFV + $G_RCFV + $S_RCFV + $T_RCFV + $C_RCFV + $N_RCFV + $Q_RCFV + $Y_RCFV + $D_RCFV + $E_RCFV + $K_RCFV + $R_RCFV + $H_RCFV;
#print $Length[0];
my $nRCFV = $RCFV/(($Length[0]**-0.5)*2000);
#my $baseline = $Length[0];
#print $baseline**-0.5, "\n";

#print CSRCFV_OUT "character RCFV values:\ncsRCFV(A)\t", $A_RCFV, "\ncsRCFV(V)\t", $V_RCFV, "\ncsRCFV(L)\t", $L_RCFV, "\ncsRCFV(I)\t", $I_RCFV, "\ncsRCFV(P)\t", $P_RCFV, "\ncsRCFV(M)\t", $M_RCFV, 
#"\ncsRCFV(F)\t", $F_RCFV, "\ncsRCFV(W)\t", $W_RCFV, "\ncsRCFV(G)\t", $G_RCFV, "\ncsRCFV(S)\t", $S_RCFV, "\ncsRCFV(T)\t", $T_RCFV, "\ncsRCFV(C)\t", $C_RCFV, "\ncsRCFV(N)\t", $N_RCFV, "\ncsRCFV(Q)\t", $Q_RCFV, 
#"\ncsRCFV(Y)\t", $Y_RCFV, "\ncsRCFV(D)\t", $D_RCFV, "\ncsRCFV(E)\t", $E_RCFV, "\ncsRCFV(K)\t", $K_RCFV, "\ncsRCFV(R)\t", $R_RCFV, "\ncsRCFV(H)\t", $H_RCFV, "\n";
print NCSRCFV_OUT "character nRCFV values:\nncsRCFV(A)\t", $A_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(V)\t", $V_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(L)\t", $L_RCFV/(($Length[0]**-0.5)*100), 
"\nncsRCFV(I)\t", $I_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(P)\t", $P_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(M)\t", $M_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(F)\t", $F_RCFV/(($Length[0]**-0.5)*100), 
"\nncsRCFV(W)\t", $W_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(G)\t", $G_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(S)\t", $S_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(T)\t", $T_RCFV/(($Length[0]**-0.5)*100), 
"\nncsRCFV(C)\t", $C_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(N)\t", $N_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(Q)\t", $Q_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(Y)\t", $Y_RCFV/(($Length[0]**-0.5)*100), 
"\nncsRCFV(D)\t", $D_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(E)\t", $E_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(K)\t", $K_RCFV/(($Length[0]**-0.5)*100), "\nncsRCFV(R)\t", $R_RCFV/(($Length[0]**-0.5)*100), 
"\nncsRCFV(H)\t", $H_RCFV/(($Length[0]**-0.5)*100), "\n";


print RCFV_OUT "total RCFV\nRCFV\t", $RCFV, "\nnRCFV\t", $nRCFV, "\n";
}

else {
print "Are you sure you remembered to specify protein or dna for amino acid or nucleotide data?" and die $!;
}

&end;

sub end{
		
		TIMER:
		# set timer
		my ( $user, $system, $cuser, $csystem ) = times ;
		
		print "\n\n\t-------------------------------------------\n\tRCFV CALCULATION COMPLETE Ta'ra\n\t-------------------------------------------\n\t";
		
		print <<TIME;
		
		***  time used: $user sec  ***
		
TIME
		

		
		
		
}

