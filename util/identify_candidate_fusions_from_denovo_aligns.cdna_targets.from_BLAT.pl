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
    



    if (0) {
        
        ## Tier and merge strategy

    
        my %tiered_n_merged_gene_alignments = &tier_and_merge_best_alignments(\%trans_to_gene_alignments);
        
        foreach my $trans_id (keys %tiered_n_merged_gene_alignments) {
            
            my @hits = @{$tiered_n_merged_gene_alignments{$trans_id}};
            @hits = sort {$a->{trans_lend}<=>$b->{trans_lend}} @hits;
            
            if (scalar @hits < 2) {
                # shouldn't happen as filtering is done before this.
                die "Error, have fewer than 2 candidate fusion parts";
            }
            
            ## Report fusion pairs:
            my $geneA = shift @hits;
            while (@hits) {
                my $geneB = shift @hits;
                
                if ($geneA->{gene_orient} eq $geneB->{gene_orient}) {
                    
                    my ($left_gene, $right_gene) = ($geneA->{gene_orient} eq '+') ? ($geneA, $geneB) : ($geneB, $geneA);
                    
                    my $fusion_name = join("--", $left_gene->{gene_name}, $right_gene->{gene_name});
                    
                    print join("\t", $fusion_name, $trans_id,
                               $left_gene->{gene_name}, $left_gene->{trans_lend}, $left_gene->{trans_rend}, $left_gene->{gene_orient}, $left_gene->{per_id} . "\%ID",
                               $right_gene->{gene_name}, $right_gene->{trans_lend}, $right_gene->{trans_rend}, $right_gene->{gene_orient}, $right_gene->{per_id} . "\%ID") . "\n";
                    
                }
                
                $geneA = $geneB;
            }
        }
    }
    
    exit(0);


}


####
sub report_fusion_candidates {
    my ($trans_to_gene_alignments_href) = @_;

    my @final_fusion_preds;


    foreach my $trans_id (keys %$trans_to_gene_alignments_href) {
        
        my @aligns = @{$trans_to_gene_alignments_href->{$trans_id}};
        
        @aligns = sort {$a->{trans_lend}<=>$b->{trans_lend}} @aligns;

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


                if ($seg_len_j < $MIN_ANCHOR_LEN) { next; }

                my $pct_overlap = &get_percent_overlap([$lend_i, $rend_i], [$lend_j, $rend_j]);

                my $delta = abs($lend_j -1 - $rend_i);
                my ($overlap_len, $gap_len) = ($lend_j >= $rend_i) ? ($delta, 0) : (0, $delta);
                

                if ($pct_overlap < $MAX_ALLOWED_PCT_OVERLAP) {

                    my ($left_gene, $right_gene) = ($align_i->{gene_orient} eq '+') ? ($align_i, $align_j) : ($align_j, $align_i);
                    
                    my $score = ($seg_len_i * $per_id_i) + ($seg_len_j * $per_id_j) 
                        - ($overlap_len * $per_id_i/2) - ($overlap_len * $per_id_j/2)
                        - ($gap_len * 100); 

                    push (@fusion_preds, [$left_gene, $right_gene, $score]);
                    
                }
            }
        }
        
        @fusion_preds = reverse sort {$a->[2]<=>$b->[2]} @fusion_preds;

        my %seen;
        
        foreach my $fusion_pred (@fusion_preds) {

            my ($left_gene, $right_gene, $score) = @$fusion_pred;
            
            ## report candidate
            my $fusion_name = join("--", $left_gene->{gene_name}, $right_gene->{gene_name});
            
            if ($seen{$fusion_name}) { next; }
            $seen{$fusion_name} = 1;
            
            push (@final_fusion_preds, { 
                
                trans_id => $trans_id,
                left_gene => $left_gene,
                right_gene => $right_gene,
                fusion_name => $fusion_name,
                score => $score,
                
                  } );
            
        }
        
    }
    

    
    # sort lexically by fusion name - easier to compare lists
    
    @final_fusion_preds = sort {$a->{fusion_name} cmp $b->{fusion_name}} @final_fusion_preds;
    
    
    foreach my $fusion_pred (@final_fusion_preds) {
        
        my ($fusion_name, $trans_id, $left_gene, $right_gene, $score) = ($fusion_pred->{fusion_name},
                                                                         $fusion_pred->{trans_id},
                                                                         $fusion_pred->{left_gene},
                                                                         $fusion_pred->{right_gene},
                                                                         $fusion_pred->{score});
        
        print join("\t", $fusion_name, $trans_id,
                   $left_gene->{gene_name}, $left_gene->{trans_lend}, $left_gene->{trans_rend}, $left_gene->{gene_orient}, $left_gene->{per_id} . "\%ID",
                   $right_gene->{gene_name}, $right_gene->{trans_lend}, $right_gene->{trans_rend}, $right_gene->{gene_orient}, $right_gene->{per_id} . "\%ID",
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
sub tier_and_merge_best_alignments {
    my ($trans_to_gene_alignments_href) = @_;

    my %trans_id_to_tiers;

    foreach my $trans_id (keys %$trans_to_gene_alignments_href) {
        
        my @hits = @{$trans_to_gene_alignments_href->{$trans_id}};
        
        my @merged_gene_hits = @hits; #&merge_hits_by_gene(@hits);

        unless (scalar @merged_gene_hits > 1) {
            next;
        }

        foreach my $hit (@merged_gene_hits) {
            $hit->{score} = ($hit->{trans_rend} - $hit->{trans_lend} + 1) * $hit->{per_id};
        }

        @hits = reverse sort {$a->{score}<=>$b->{score}} @merged_gene_hits;
        
        
        my @tiered_hits;
        foreach my $hit (@hits) {
            print STDERR "Testing for tier addition.  Current tier is: " . Dumper(\@tiered_hits) if $DEBUG;

            if (! &overlap_existing_tier_element(\@tiered_hits, $hit)) {
                push (@tiered_hits, $hit);
            }
        }
        if (scalar @tiered_hits > 1) {
            $trans_id_to_tiers{$trans_id} = [@tiered_hits];
        }
    }

    return(%trans_id_to_tiers);
    
}


####
sub overlap_existing_tier_element {
    my ($tiered_aref, $hit) = @_;

    foreach my $tier_entry (@$tiered_aref) {
        
        my $tier_gene_name = $tier_entry->{gene_name};
        my $tier_lend = $tier_entry->{trans_lend};
        my $tier_rend = $tier_entry->{trans_rend};
        my $tier_seg_len = $tier_rend - $tier_lend + 1;


        my $hit_gene_name = $hit->{gene_name};
        my $hit_lend = $hit->{trans_lend};
        my $hit_rend = $hit->{trans_rend};
        my $hit_len = $hit_rend - $hit_lend + 1;


        print STDERR "TESTING $tier_gene_name vs. $hit_gene_name : ($tier_lend-$tier_rend) vs. ($hit_lend, $hit_rend)\n" if $DEBUG;

        
        if (&coordsets_overlap([$tier_lend, $tier_rend], [$hit_lend, $hit_rend])) {

            my $smaller_len = ($tier_seg_len < $hit_len) ? $tier_seg_len : $hit_len;

            my $pct_overlap = &nucs_in_common($tier_lend, $tier_rend, $hit_lend, $hit_rend) / $smaller_len * 100;
            
            print STDERR "($tier_lend-$tier_rend) vs. ($hit_lend, $hit_rend)\tPCT_Overlap: $pct_overlap\n" if $DEBUG;

            if ($pct_overlap > $MAX_ALLOWED_PCT_OVERLAP) {
                return(1);
            }
        }
        else {
            print STDERR "($tier_lend-$tier_rend) vs. ($hit_lend, $hit_rend)\tNO_overlap\n" if $DEBUG;
        }

    }
    return(0); # no sufficient overlap found
    
}
        
        



####
sub merge_hits_by_gene {
    my @hits = @_;

    my %gene_to_hit_list;

    foreach my $hit (@hits) {
        my $gene_name = $hit->{gene_name};
        
        push (@{$gene_to_hit_list{$gene_name}}, $hit);
    }


    ## merge them
    my @merged_gene_hits;

    foreach my $gene_name (keys %gene_to_hit_list) {
        my @hits = @{$gene_to_hit_list{$gene_name}};
        
        ## todo: check for actual overlaps
        ## first, keeping it super easy and selecting range of matches 

        my @lend_coords;
        my @rend_coords;
        
        my $sum_len = 0;
        my $sum_per_id_n_len = 0;

        foreach my $hit (@hits) {
            push (@lend_coords, $hit->{trans_lend});
            push (@rend_coords, $hit->{trans_rend});
        
            my $seg_len = $hit->{trans_rend} - $hit->{trans_lend} + 1;
            my $per_id = $hit->{per_id};

            $sum_per_id_n_len += $seg_len * $per_id;

            $sum_len += $seg_len;
        }

        @lend_coords = sort {$a<=>$b} @lend_coords;
        @rend_coords = sort {$a<=>$b} @rend_coords;

        my $range_lend = shift @lend_coords;
        my $range_rend = pop @rend_coords;


        my $avg_per_id = sprintf("%.1f", $sum_per_id_n_len / $sum_len);
        my $hit = shift @hits;
        $hit->{trans_lend} = $range_lend;
        $hit->{trans_rend} = $range_rend;
        $hit->{per_id} = $avg_per_id;

        push (@merged_gene_hits, $hit);
    }

    return(@merged_gene_hits);
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
        
        $trans_id .= "::$orient";
        
        my $gene_name = $trans_id_to_gene_id_href->{$target_trans_id} or die "Error, no gene name for $target_trans_id";
        
        push (@{$trans_to_gene_hits{$trans_id}}, { gene_name => $gene_name,
                                                   trans_lend => $trans_end5,
                                                   trans_rend => $trans_end3,
                                                   gene_orient => $orient,
                                                   per_id => $per_id,
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
