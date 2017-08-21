#!/usr/bin/env perl

use strict; use warnings;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

###################################################################
######### Setup ###################################################
###################################################################
my $error_message = "Usage: $0 --in_mp <mpileup> --out_ct <counts_tot.tsv> --out_cp <counts_pos.tsv>\n";
my ($in_mp, $out_ct, $out_cp);

GetOptions (
    "in_mp=s" => \$in_mp,
    "out_ct=s" => \$out_ct,
    "out_cp=s" => \$out_cp
) or die $error_message;

if (!($in_mp && $out_ct && $out_cp)) { die $error_message }

my $baseq_min = 30;
my $cov_min = 1;
my $cov_max = 100;
my $soft_masked = 1;

###################################################################
######### Main ####################################################
###################################################################
# total
my $pos_result = &get_position_count($in_mp, $soft_masked, $cov_min, $cov_max, $baseq_min);
my $tot_result = &get_total_count($pos_result);
my $tot_final = &summarise_total($tot_result);
open(OUT_CT, ">", $out_ct)  or die "Can't open $out_ct for writing";
&print_total($tot_final, \*OUT_CT);
close(OUT_CT);

## position
open(OUT_CP, ">", $out_cp)  or die "Can't open $out_cp for writing";
my $pos_final = &summarise_pos($pos_result);
&print_pos($pos_final, \*OUT_CP);
close(OUT_CP);

###################################################################
######### Functions ###############################################
###################################################################

sub print_pos {
    my ($final, $fh) = @_;
    print $fh "Type\tPosition\tAbsolute\tRelative\n";
    foreach my $type (sort keys %{$final}) {
        my $positions = ${$final}{$type};
        foreach my $pos (keys %$positions) {
            my $abs = $positions->{$pos}->{'absolute'};
            my $rel = $positions->{$pos}->{'relative'};
            print $fh "$type\t" . $pos . "\t" . $abs . "\t" . $rel . "\n";
        }
    }
}

sub print_total {
    my ($final, $fh) = @_;
    print $fh "Type\tTotal_Type\tTotal\n";
    foreach my $type (sort keys %{$final}) {
        print $fh "$type\t" . $final->{$type}->{"total_type"} . "\t" . $final->{$type}->{"total"} . "\n";
    }
}

sub summarise_total {
    my ($result) = @_;
    my %final;

    my $mutations = %{$result}{"type"};
    my $nt = %{$result}{"nt"};

    foreach my $type (keys %$mutations) {
        my $total_type = %{$mutations}{$type};
        $type =~ /(.)\_(.)/;
        my $total = %{$nt}{$1};

        $final{$type}{"total_type"} = $total_type;
        $final{$type}{"total"} = $total;
    }
    return(\%final);
}

sub summarise_pos {
    my ($result) = @_;
    my %final;

    my $positions = %{$result}{"type"};
    foreach my $pos (keys %$positions) {
        my $mutations = $$positions{$pos};
        foreach my $type (keys %$mutations) {
            my $total_type = $$result{"type"}{$pos}{$type};
            $type =~ /(.)\_(.)/;
            my $total = $$result{"nt"}{$pos}{$1};
            my $rel_value = $total_type / $total;
            $final{$type}{$pos}{"relative"} = $rel_value;
            $final{$type}{$pos}{"absolute"} = $total_type;
        }
    }
    return(\%final);
}

sub get_total_count {
    my ($href) = @_;
    my (%nt_total, %type_total);

    my $positions = $href->{'type'};
    foreach my $pos (keys %$positions) {
        my $mutations = $positions->{$pos};
        foreach my $type (keys %$mutations) {
            $type_total{$type} += $mutations->{$type};
        }
    }

    $positions = $href->{'nt'};
    foreach my $pos (keys %$positions) {
        my $mutations = $positions->{$pos};
        foreach my $base (keys %$mutations) {
            $nt_total{$base} += $mutations->{$base};
        }
    }
    return({ 'type' => \%type_total,
             'nt'   => \%nt_total});
}

sub get_position_count {
    my ($mpileup, $soft_masked, $cov_min, $cov_max, $baseq_min) = @_;
    my %result;

    open (MPILEUP, $mpileup) or die "Can't open $mpileup!!\n";
    while (<MPILEUP>) {
        s/\r|\n//g;
        # chr | loc       | ref | cov | dbases | baseq | mapq  | readpos
        # 1   | 215622909 | T   | 5   | .,...  | :<@=@ | F]]S] | 47,18,5,3,2
        my ($chr, $loc, $ref, $cov, $dbases, $baseq, $mapq, $readpos) = split /\t/;
        my $cbases = clean_bases($dbases);

        my @tmp = split //, $cbases; # get each base into an array element
        my @readpos_a = split/,/, $readpos; # same for basepos
        my $length_mutation = @tmp; # can just use length of $cbases probably...

        if ($soft_masked == 0) {
            $ref = uc($ref);
        }

        # checkpoint
        if ($cov == $length_mutation && $cov >= $cov_min && $cov <= $cov_max && $ref =~ /[ACTG]/) {
            my @cbase_a = split //, $cbases;
            my @baseq_a = split //, $baseq;

            # Loop over cbase
            for (my $i = 0; $i < @cbase_a; $i++) {
                my $qscore = ord($baseq_a[$i]) - 33;

                if ($qscore > $baseq_min) {
                    my $ch = $cbase_a[$i];
                    my $rp = $readpos_a[$i];

                    if ($ch =~ /[ATCG]/) { # mismatch on forward strand
                        my $type = $ref . "_" . $ch;
                        $result{"type"}{$rp}{$type}++;
                        $result{"nt"}{$rp}{$ref}++;
                    } elsif ($ch =~ /[\+\-]/ && $rp > 5) {
                        my $type = $ref . "_"  . $ch;
                        $result{"type"}{$rp}{$type}++;
                        $result{"nt"}{$rp}{$ref}++;
                    } elsif ($ch eq ".") { # ref match on forward strand
                        my $type = $ref . "_" . $ref;
                        $result{"type"}{$rp}{$type}++;
                        $result{"nt"}{$rp}{$ref}++;
                    }
                }
            } # end for loop
        } # end checkpoint
    } # end while loop
    close MPILEUP;
    return(\%result);
}

sub clean_bases {
    my ($bases) = @_; #e.g.: .$,...
    $bases =~ s/\^.//g; #^X => start of read followed by mapq
    $bases =~ s/\$//g; #$ => end of read

    # indel found
    if ($bases =~ /\+(\d+)/ || $bases =~ /\-(\d+)/) {
        # insertion
        while ($bases =~ /\+(\d+)/) {
            my $size = $1;
            my $string_to_remove;
            for (my $i = 0; $i < $size; $i++) {
                $string_to_remove = $string_to_remove . ".";
            }
            $string_to_remove = '.\+\d+' . $string_to_remove;
            $bases =~ s/$string_to_remove/\+/;
        }
        # deletion
        while ($bases =~ /\-(\d+)/) {
            my $size = $1;
            my $string_to_remove;
            for (my $i=0; $i < $size; $i++) {
                $string_to_remove = $string_to_remove . ".";
            }
            $string_to_remove = '.\-\d+' . $string_to_remove;
            $bases =~ s/$string_to_remove/\-/;
        }
    }
    return $bases;
}
