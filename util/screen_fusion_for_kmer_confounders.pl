#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw (min max);
use POSIX;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_retriever;


my $usage = "usage: $0 discafuse_preds.prelim.dat reference_cdna.fasta KMER_SIZE WINDOW_SIZE\n\n";

my $fusion_preds_file = $ARGV[0] or die $usage;
my $ref_cdna_fasta = $ARGV[1] or die $usage;
my $KMER_SIZE = $ARGV[2] or die $usage;
my $SEARCH_DIST_AROUND_BREAKPOINT = $ARGV[3] or die $usage;


main: {

    my $fasta_retriever = new Fasta_retriever($ref_cdna_fasta);
    
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    open (my $ofh, ">$fusion_preds_file.kmer_confound_info") or die "Error, cannot write to file $fusion_preds_file.kmer_confound_info";

    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        
        my $hit_info_left = $x[3];
        my $hit_info_right = $x[8];

        # all the fusions are already oriented left to right at the gene level, so take breakpoint coords according to the (+)--(+) orientation.
        my ($left_acc, $left_lend, $left_rend, $left_orient) = split(/:/, $hit_info_left);
        my $left_break_coord = $left_rend;
        
        my ($right_acc, $right_lend, $right_rend, $right_orient) = split(/:/, $hit_info_right);
        my $right_break_coord = $right_lend;
        
        my $left_trans_seq = $fasta_retriever->get_seq($left_acc);
        my $right_trans_seq = $fasta_retriever->get_seq($right_acc);

        
        my %shared_kmers = &get_shared_kmers($left_trans_seq, $right_trans_seq, $left_break_coord, $right_break_coord, $SEARCH_DIST_AROUND_BREAKPOINT);
        
        if (%shared_kmers) {
            ## filter and annotate
            my $kmer_text = "";
            foreach my $kmer (keys %shared_kmers) {
                my ($A_pos_list_aref, $B_pos_list_aref) = @{$shared_kmers{$kmer}};
                $kmer_text .= "\t$kmer(" . join(",", @$A_pos_list_aref) . "|" . join(",", @$B_pos_list_aref) . ")";

            }
            print $ofh "#$line" . $kmer_text . "\n";
        }
        else {
            print $ofh $line . "\n";
            print $line . "\n"; # unfiltered goes to stdout.
        }
        
    }
    close $fh;
    close $ofh;


    exit(0);
    
}

####
sub get_shared_kmers {
    my ($seqA, $seqB, $brkpt_A, $brktp_B, $windowsize) = @_;
    
    my %kmersA = &get_kmers($seqA, $brkpt_A, $windowsize);

    my %kmersB = &get_kmers($seqB, $brktp_B, $windowsize);
    
    my %joint_kmers;
    foreach my $kmer (keys %kmersA) {
        if (exists $kmersB{$kmer}) {
            my $A_pos_list = $kmersA{$kmer};
            my $B_pos_list = $kmersB{$kmer};
            
            $joint_kmers{$kmer} = [$A_pos_list, $B_pos_list];
        }
    }

    return(%joint_kmers);
}


####
sub get_kmers {
    my ($seq, $centerpt, $windowsize) = @_;

    my %kmers_to_pos_list;

    
    my $beg_pos = max(0, $centerpt - ceil($windowsize/2.0));
    my $end_pos = min(length($seq), $centerpt + ceil($windowsize/2.0));
    
    # print "WINDOW ($centerpt): $beg_pos to $end_pos\n";
    
    for (my $i = $beg_pos; $i <= $end_pos - $KMER_SIZE; $i++) {
        
        my $kmer = uc substr($seq, $i, $KMER_SIZE);
        
        push (@{$kmers_to_pos_list{$kmer}}, $i);
    }

    return(%kmers_to_pos_list);
}
