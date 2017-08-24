#!/usr/bin/env perl

use strict; use warnings;
use Getopt::Long qw(GetOptions);
use Data::Dumper;

###################################################################
######### Setup ###################################################
###################################################################
my $error_message = "Usage: $0 --in_mp <mpileup> --out_ct <counts_tot.tsv> --out_cp <counts_pos.tsv> --out_cc <counts_con.tsv>\n";
my ($in_mp, $out_ct, $out_cp, $out_cc);

GetOptions (
    "in_mp=s" => \$in_mp,
    "out_ct=s" => \$out_ct,
    "out_cp=s" => \$out_cp,
    "out_cc=s" => \$out_cc
) or die $error_message;

if (!($in_mp && $out_ct && $out_cp && $out_cc)) { die $error_message }

my $baseq_min = 30;
my $cov_min = 1;
my $cov_max = 100;
my $soft_masked = 1;
my $context = 3;

###################################################################
######### Main ####################################################
###################################################################
# total
my $con_result = &get_context_count($in_mp, $soft_masked, $cov_min, $cov_max, $baseq_min, $context);
my $pos_result = &get_position_count2($con_result);
my $tot_result = &get_total_count($pos_result);

my $tot_final = &summarise_tot($tot_result);
open(OUT_CT, ">", $out_ct)  or die "Can't open $out_ct for writing";
&print_tot($tot_final, \*OUT_CT);
close(OUT_CT);

# position
my $pos_final = &summarise_pos($pos_result);
open(OUT_CP, ">", $out_cp)  or die "Can't open $out_cp for writing";
&print_pos($pos_final, \*OUT_CP);
close(OUT_CP);

## context
my $con_final = &summarise_con($con_result);
open(OUT_CC, ">", $out_cc)  or die "Can't open $out_cc for writing";
&print_con($con_final, \*OUT_CC);
close(OUT_CC);


###################################################################
######### Functions ###############################################
###################################################################

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

sub print_tot {
    my ($final, $fh) = @_;
    print $fh "Type\tTotal_Type\tTotal\n";
    foreach my $type (sort keys %{$final}) {
        my $total_type = $final->{$type}->{"total_type"};
        my $total = $final->{$type}->{"total"};
        print $fh "$type\t$total_type\t$total\n";
    }
}

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

sub print_con {
    my ($final, $fh) = @_;
    print $fh "Type\tPosition\tContext\tAbsolute\tRelative\n";
    foreach my $type (sort keys %{$final}) {
        my $positions = $final->{$type};
        foreach my $pos (keys %$positions) {
            my $contexts = $positions->{$pos};
            foreach my $con (keys %$contexts) {
                my $abs = $positions->{$pos}->{$con}->{'absolute'};
                my $rel = $positions->{$pos}->{$con}->{'relative'};
                print $fh "$type\t$pos\t$con\t$abs\t$rel\n";
            }
        }
    }
}

sub get_position_count2 {
    my ($href) = @_;
    my (%nt_total, %type_total);

    my $positions = $href->{'type'};
    foreach my $pos (keys %$positions) {
        my $mutations = $positions->{$pos};
        foreach my $type (keys %$mutations) {
            my $contexts = $mutations->{$type};
            foreach my $con (keys %$contexts) {
                $type_total{$pos}{$type} += $mutations->{$type}->{$con};
            }
        }
    }

    $positions = $href->{'nt'};
    foreach my $pos (keys %$positions) {
        my $mutations = $positions->{$pos};
        foreach my $type (keys %$mutations) {
            my $contexts = $mutations->{$type};
            foreach my $con (keys %$contexts) {
                $nt_total{$pos}{$type} += $mutations->{$type}->{$con};
            }
        }
    }
    return({ 'type' => \%type_total,
             'nt'   => \%nt_total});

}

# not used now
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

sub get_nt {
    my ($line) =@_;
    my ($chr, $loc, $ref, $cov, $dbases, $baseq, $mapq, $readpos) = split /\t/, $line;
    return ($ref);
}

### get_context_count($in_mp, $soft_masked, $cov_min, $cov_max, $baseq_min, $context);
sub get_context_count {
    my ($mpileup, $soft_masked, $cov_min, $cov_max, $baseq_min, $context) = @_;
    my %result;
    my $l1;
    my $l2;
    my $l3;

    open (MPILEUP, $mpileup) or die "Can't open $mpileup\n";
    while (my $line = <MPILEUP>) {
        $line =~ s/\r|\n//g;
        $l1 = $l2;
        $l2 = $l3;
        $l3 = $line;
        if ($l1 && $l2 && $l3) {
            my $nt1 = get_nt($l1);
            my $nt2 = get_nt($l2);
            my $nt3 = get_nt($l3);

            my $seq_all = "no_context";
            if ($context == 1){$seq_all = $nt1 . "_base"; }
            if ($context == 2){$seq_all = "base_" . $nt3; }
            if ($context == 3){$seq_all = $nt1 . "_base_" . $nt3; }

            my $seq = uc($seq_all);
            my ($chr, $loc, $ref, $cov, $dbases, $baseq, $mapq, $readpos) = split /\t/, $l2;
            my $cbases = clean_bases($dbases);
            my @tmp = split //, $cbases;
            my @readpos_a = split/,/, $readpos;
            my $length_mutation = @tmp;
            my @baseq_a = split //, $baseq;

            if ($soft_masked == 0) {
                $ref = uc($ref);
            }

            # checkpoint
            if ($cov == $length_mutation && $cov >= $cov_min && $cov <= $cov_max && $ref =~ /[ACTG]/) {
                my @cbase_a = split //, $cbases;

                # loop over cbase
                for (my $i = 0; $i < @cbase_a; $i++) {
                    my $qscore = ord($baseq_a[$i]) - 33;

                    if ($qscore > $baseq_min) {
                        my $ch = $cbase_a[$i];
                        my $rp = $readpos_a[$i];

                        if ($ch =~ /[ATCG]/) {
                            my $type = $ref . "_" . $ch;
                            $result{"type"}{$rp}{$type}{$seq}++;
                            $result{"nt"}{$rp}{$ref}{$seq}++;
                        } elsif ($ch =~ /[\+\-]/ && $rp > 5) {
                            my $type = $ref . "_"  . $ch;
                            $result{"type"}{$rp}{$type}{$seq}++;
                            $result{"nt"}{$rp}{$ref}{$seq}++;
                        } elsif ($ch eq ".") {
                            my $type = $ref."_".$ref;
                            $result{"type"}{$rp}{$type}{$seq}++;
                            $result{"nt"}{$rp}{$ref}{$seq}++;
                        }
                    }
                } # end for loop
            } # end checkpoint
        }
    } # end while loop
    close MPILEUP;
    return(\%result);
}

sub summarise_tot {
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

sub summarise_con {
    my ($result) = @_;
    my %final;
    my $positions = %{$result}{"type"};
    foreach my $pos (keys %$positions) {
        my $mutations = $$positions{$pos};
        foreach my $type (keys %$mutations) {
            my $seqs = $$mutations{$type};
            foreach my $seq (keys %$seqs) {
                my $total_type = $$result{"type"}{$pos}{$type}{$seq};
                $type =~/(.)\_(.)/;
                my $total = $$result{"nt"}{$pos}{$1}{$seq};
                my $rel_value = $total_type / $total;
                $final{$type}{$pos}{$seq}{"relative"} = $rel_value;
                $final{$type}{$pos}{$seq}{"absolute"} = $total_type;
            }
        }
    }
    return(\%final);
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
