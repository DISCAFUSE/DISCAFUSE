#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/../PerlLib";
require "overlapping_nucs.ph";

my $MAX_ALLOWED_PCT_OVERLAP = 10;
my $MIN_PER_ID = 90;
my $MIN_ANCHOR_LEN = 30;
my $PCT_SCORE_CONSIDER_EQUIV = 98;

my $OVERLAP_PENALTY = 2;
my $GAP_PENALTY = 2;

my $DEBUG = 1;

my $usage = "\n\n\tusage: $0 blat.gff3 cdna_targets.fasta\n\n";

my $blat_gff3 = $ARGV[0] or die $usage;
my $cdna_targets_fasta = $ARGV[1] or die $usage;


main: {


    my $headers_file = "$cdna_targets_fasta.headers";
    unless (-s $headers_file) {
        die "Error, cannot locate $headers_file";
    }
    
    my %trans_id_to_gene_id = &get_trans_to_gene_mapping($headers_file);
    
    my %trans_to_gene_alignments = &get_trans_to_gene_alignments($blat_gff3, \%trans_id_to_gene_id);
    
    &report_fusion_candidates(\%trans_to_gene_alignments);
    
        
    exit(0);


}


####
sub report_fusion_candidates {
    my ($trans_to_gene_alignments_href) = @_;

    my @final_fusion_preds;


    foreach my $trans_id (keys %$trans_to_gene_alignments_href) {
        
        my @all_aligns = @{$trans_to_gene_alignments_href->{$trans_id}};
        
        my @candidate_fusions_to_filter;
        

        foreach my $gene_orient ('+', '-') {
            
            my @aligns = grep { $_->{gene_orient} eq $gene_orient } @all_aligns;

            unless (@aligns) { next; }
            
            @aligns = sort {$a->{trans_lend}<=>$b->{trans_lend}} @aligns;
                        
            my @fusion_preds = &aligns_to_fusion_preds(@aligns);
            
            # sorting by score
            @fusion_preds = reverse sort {$a->[2]<=>$b->[2]} @fusion_preds;
            
            my %seen;
            
                        
            foreach my $fusion_pred (@fusion_preds) {
                
                my ($left_gene, $right_gene, $score) = @$fusion_pred;
                
                ## report candidate
                my $fusion_name = join("--", $left_gene->{gene_name}, $right_gene->{gene_name});
                
                if ($seen{$fusion_name}) { next; }
                $seen{$fusion_name} = 1;
                
                push (@candidate_fusions_to_filter, { 
                
                    trans_id => $trans_id,
                    left_gene => $left_gene,
                    right_gene => $right_gene,
                    fusion_name => $fusion_name,
                    score => $score,
                    
                      } );
                
            }
        } # end of foreach orient
        
        unless (@candidate_fusions_to_filter) { next; }


        my @filtered_fusion_preds = &filter_fusion_candidates_by_score(@candidate_fusions_to_filter);

        push (@final_fusion_preds, @filtered_fusion_preds);
        
        
    } # end of foreach trans
    

    
    # sort lexically by fusion name - easier to compare lists
    
    @final_fusion_preds = sort {$a->{fusion_name} cmp $b->{fusion_name}} @final_fusion_preds;
    
    
    foreach my $fusion_pred (@final_fusion_preds) {
        
        my ($fusion_name, $trans_id, $left_gene, $right_gene, $score) = ($fusion_pred->{fusion_name},
                                                                         $fusion_pred->{trans_id},
                                                                         $fusion_pred->{left_gene},
                                                                         $fusion_pred->{right_gene},
                                                                         $fusion_pred->{score});
        
        print join("\t", $fusion_name, $trans_id,
                   
                   ## left gene info
                   $left_gene->{gene_name}, 
                   join(":", $left_gene->{hit_id}, $left_gene->{hit_lend}, $left_gene->{hit_rend}, $left_gene->{gene_orient}),
                   $left_gene->{trans_lend}, $left_gene->{trans_rend}, $left_gene->{per_id} . "\%ID",
                   
                   ## right gene info
                   $right_gene->{gene_name},
                   join(":", $right_gene->{hit_id}, $right_gene->{hit_lend}, $right_gene->{hit_rend}, $right_gene->{gene_orient}),
                   $right_gene->{trans_lend}, $right_gene->{trans_rend}, $right_gene->{gene_orient}, $right_gene->{per_id} . "\%ID",
                   
                   $score) . "\n";
        
    }
    
    

    return;
}


####
sub get_percent_overlap {
    my ($coords_A_aref, $coords_B_aref) = @_;

    if (&coordsets_overlap($coords_A_aref, $coords_B_aref)) {
        
        my $len_A = $coords_A_aref->[1] - $coords_A_aref->[0] + 1;
        my $len_B = $coords_B_aref->[1] - $coords_B_aref->[0] + 1;
        
        my $shorter_len = ($len_A < $len_B) ? $len_A : $len_B;

        my $overlap_len = &nucs_in_common($coords_A_aref->[0], $coords_A_aref->[1], $coords_B_aref->[0], $coords_B_aref->[1]);

        my $pct_overlap = $overlap_len / $shorter_len * 100;
        
        return($pct_overlap);

    }
    else {
        return(0);
    }
}
                




####
sub get_trans_to_gene_alignments {
    my ($blat_gff3, $trans_id_to_gene_id_href) = @_;

    my %trans_to_gene_hits;

    my $debug_ofh;
    if ($DEBUG) {
        open ($debug_ofh, ">$blat_gff3.gene_mappings") or die $!;
    }


    open (my $fh, $blat_gff3) or die "Error, cannot open file $blat_gff3";
    while (<$fh>) {
        chomp;
        my $line = $_;
        
        my @x = split(/\t/);
        my $target_trans_id = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $per_id = $x[5];

        if ($per_id < $MIN_PER_ID) { next; }

        my $orient = $x[6];
        my $info = $x[8];

        $info =~ /Target=(\S+) (\d+) (\d+)/ or die "Error, cannot parse hit info from $info";
        
        my $trans_id = $1;
        my $trans_end5 = $2;
        my $trans_end3 = $3;
        
        #$trans_id .= "::$orient";
        
        my $gene_name = $trans_id_to_gene_id_href->{$target_trans_id} or die "Error, no gene name for $target_trans_id";
        
        push (@{$trans_to_gene_hits{$trans_id}}, { gene_name => $gene_name,
                                                   trans_lend => $trans_end5,
                                                   trans_rend => $trans_end3,
                                                   gene_orient => $orient,
                                                   per_id => $per_id,
                                                   
                                                   hit_lend => $lend,
                                                   hit_rend => $rend,
                                                   hit_id => $target_trans_id,
              }
            );


        print $debug_ofh join("\t", $gene_name, $trans_id, $trans_end5, $trans_end3, $orient, $per_id, $line) . "\n" if $DEBUG;

    }
    print STDERR "\n-done parsing gff3\n\n" if $DEBUG;
    close $debug_ofh if $DEBUG;
    
    return(%trans_to_gene_hits);
}


####
sub get_trans_to_gene_mapping {
    my ($headers_file) = @_;
    
    my %trans_id_to_gene_id;

    open (my $fh, $headers_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($trans_id, $gene_id, $gene_name) = split(/\s+/);
        
        if ($gene_name) {
            $trans_id_to_gene_id{$trans_id} = $gene_name;
        }
        else {
            $trans_id_to_gene_id{$trans_id} = $gene_id;
        }

    }
    close $fh;


    return(%trans_id_to_gene_id);
}


####
sub aligns_to_fusion_preds {
    my @aligns = @_;

    my @fusion_preds;

    for (my $i = 0; $i < $#aligns; $i++) {
        
        my $align_i = $aligns[$i];
        
        my $gene_name_i = $align_i->{gene_name};
        my $lend_i = $align_i->{trans_lend};
        my $rend_i = $align_i->{trans_rend};
        my $per_id_i = $align_i->{per_id};
        my $seg_len_i = $rend_i - $lend_i + 1;
        
        if ($seg_len_i < $MIN_ANCHOR_LEN) { next; }
        
        
        for (my $j = $i + 1; $j <= $#aligns; $j++) {
            
            my $align_j = $aligns[$j];
            
            my $gene_name_j = $align_j->{gene_name};
            my $lend_j = $align_j->{trans_lend};
            my $rend_j = $align_j->{trans_rend};
            my $per_id_j = $align_j->{per_id};
            my $seg_len_j = $rend_j - $lend_j + 1;
            
            
            if ($gene_name_i eq $gene_name_j) { next; } # no selfies!

            if ($seg_len_j < $MIN_ANCHOR_LEN) { next; }
            
            my $pct_overlap = &get_percent_overlap([$lend_i, $rend_i], [$lend_j, $rend_j]);
            
            my $delta = abs($lend_j - $rend_i);
            my ($overlap_len, $gap_len) = ($lend_j <= $rend_i) ? ($delta, 0) : (0, $delta);
            
            
            if ($pct_overlap < $MAX_ALLOWED_PCT_OVERLAP) {
                
                my ($left_gene, $right_gene) = ($align_i->{gene_orient} eq '+') ? ($align_i, $align_j) : ($align_j, $align_i);
                
                my $score = ($seg_len_i * $per_id_i) + ($seg_len_j * $per_id_j) 
                    - ($overlap_len * $per_id_i * $OVERLAP_PENALTY) - ($overlap_len * $per_id_j * $OVERLAP_PENALTY)
                    - ($gap_len * 100 * $GAP_PENALTY); 
                
                #print Dumper($left_gene) . Dumper($right_gene) . " Score: $score, delta: $delta, overlap: $overlap_len, gap: $gap_len\n";
                
                push (@fusion_preds, [$left_gene, $right_gene, $score]);
                
            }
        }
    }
    
    return(@fusion_preds);
}

####
sub filter_fusion_candidates_by_score {
    my (@fusion_structs) = @_;

    @fusion_structs = reverse sort {$a->{score}<=>$b->{score}} @fusion_structs;
    
    my $top_scoring_pred = shift @fusion_structs;
    my $top_score = $top_scoring_pred->{score};
    
    unless ($top_score > 0) { return (); }

    my @best_preds = ($top_scoring_pred);

    while (@fusion_structs) {
        my $fusion_pred = shift @fusion_structs;
        if ($fusion_pred->{score} / $top_score * 100 >= $PCT_SCORE_CONSIDER_EQUIV) {
            push (@best_preds, $fusion_pred);
        }
        else {
            last;
        }
    }

    #return(@fusion_structs);
    
    return(@best_preds);
}
