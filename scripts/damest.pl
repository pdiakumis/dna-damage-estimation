#!/usr/bin/env perl

use strict; use warnings;
use Getopt::Long qw(GetOptions);

my %RC = (
    'A_C' => 'T_G', 'A_A' => 'T_T', 'A_T' => 'T_A', 'A_G' => 'T_C',
    'T_C' => 'A_G', 'T_A' => 'A_T', 'T_T' => 'A_A', 'T_G' => 'A_C',
    'C_C' => 'G_G', 'C_A' => 'G_T', 'C_T' => 'G_A', 'C_G' => 'G_C',
    'G_C' => 'C_G', 'G_A' => 'C_T', 'G_T' => 'C_A', 'G_G' => 'C_C',
    'A_-' => 'T_-', 'T_-' => 'A_-', 'C_-' => 'G_-', 'G_-' => 'C_-',
    'A_+' => 'T_+', 'T_+' => 'A_+', 'C_+' => 'G_+', 'G_+' => 'C_+'
);


my $error_message = "Usage: $0 --ct1 <counts_tot1.tsv> --ct2 <counts_tot2.tsv> --cp1 <counts_pos1.tsv> --cp2 <counts_pos2.tsv> --out_tot <total.damage> --out_pos <pos.damage> --out_con <con.damage> --id <sample_id>\n";
my ($ct1, $ct2, $cp1, $cp2, $cc1, $cc2, $out_tot, $out_pos, $out_con, $id);

GetOptions (
    "ct1=s" => \$ct1,
    "ct2=s" => \$ct2,
    "cp1=s" => \$cp1,
    "cp2=s" => \$cp2,
    "cc1=s" => \$cc1,
    "cc2=s" => \$cc2,
    "out_tot=s" => \$out_tot,
    "out_pos=s" => \$out_pos,
    "out_con=s" => \$out_con,
    "id=s" => \$id
) or die $error_message;

if (!($ct1 && $ct2 && $cp1 && $cp2 && $out_tot && $out_pos && $id)) { die $error_message }


###################################################################
######### Main ####################################################
###################################################################

# Total
my $counts_tot1 = &read_tot_counts($ct1);
my $counts_tot2 = &read_tot_counts($ct2);

open(OUT_TOT, ">", $out_tot) or die "Can't open $out_tot for writing";
&print_tot_out($counts_tot1, $counts_tot2, \%RC, \*OUT_TOT);
close(OUT_TOT);

# Position
my $counts_pos1 = &read_pos_counts($cp1);
my $counts_pos2 = &read_pos_counts($cp2);

open(OUT_POS, ">", $out_pos) or die "Can't open $out_pos for writing";
&print_pos_out($counts_pos1, $counts_pos2, \%RC, \*OUT_POS);
close(OUT_POS);

# Context
my $counts_con1 = &read_con_counts($cc1);
my $counts_con2 = &read_con_counts($cc2);

open(OUT_CON, ">", $out_con) or die "Can't open $out_con for writing";
&print_con_out($counts_con1, $counts_con2, \*OUT_CON);
close(OUT_CON);

###################################################################
######### Functions ###############################################
###################################################################
sub read_tot_counts {
    my ($count_file) = @_;
    my %final;

    open (COUNTS, $count_file) or die "Can't open $count_file\n";
    while (<COUNTS>) {
        chomp;
        unless (/Total/) {
            my ($type, $total_type, $total) = split /\t/;
            $final{$type}{"total_type"} = $total_type;
            $final{$type}{"total"} = $total;
        }
    }
    close COUNTS;
    return(\%final);
}

sub read_pos_counts {
    my ($count_file) = @_;
    my %final;

    open (COUNTS, $count_file) or die "Can't open $count_file\n";
    while (<COUNTS>) {
        chomp;
        unless (/Absolute/) {
            my ($type, $pos, $abs, $rel) = split /\t/;
            $final{$type}{$pos}{"absolute"} = $abs;
            $final{$type}{$pos}{"relative"} = $rel;
        }
    }
    close COUNTS;
    return(\%final);
}

sub read_con_counts {
    my ($count_file) = @_;
    my %final;

    open (COUNTS, $count_file) or die "Can't open $count_file\n";
    while (<COUNTS>) {
        chomp;
        unless (/Absolute/) {
            my ($type, $pos, $con, $abs, $rel) = split /\t/;
            $final{$type}{$pos}{$con}{"absolute"} = $abs;
            $final{$type}{$pos}{$con}{"relative"} = $rel;
        }
    }
    close COUNTS;
    return(\%final);
}

sub print_pos_out {
    my ($counts1, $counts2, $RC, $fh) = @_;

    foreach my $type (keys %$RC) {
        my $ref = $counts1->{$type};
        foreach my $pos (sort {$a<=>$b} keys %$ref) {
            my $value1 = $counts1->{$type}{$pos}{"relative"};
            my $value2 = $counts2->{$type}{$pos}{"relative"};
            my $abs1 = $counts1->{$type}{$pos}{"absolute"};
            my $abs2 = $counts2->{$type}{$pos}{"absolute"};
            if ($counts1->{$type}{$pos}{"relative"} && $counts2->{$type}{$pos}{"relative"}) {
                print $fh "$id\t$type\tR1\t$value1\t$abs1\t$pos\n";
                print $fh "$id\t$type\tR2\t$value2\t$abs2\t$pos\n";
            }
        }
    }
}

sub print_tot_out {
    my ($counts1, $counts2, $RC, $fh) = @_;
    my %final;

    foreach my $type (keys %$counts1) {
        my $type_rc = $RC->{$type};
        my $couple = join("-", sort($type, $type_rc));

        my $total_type1 = $counts1->{$type}{"total_type"};
        my $total_type2 = $counts2->{$type_rc}{"total_type"};
        my $total1 = $counts1->{$type}{"total"};
        my $total2 = $counts2->{$type_rc}{"total"};
        # If both files have counts for the same type
        if ($counts1->{$type}{"total_type"} && $counts2->{$type}{"total_type"}) {
            my $total_type = $total_type1 + $total_type2;
            my $total = $total1 + $total2;
            my $value = $total_type / $total;
            my $line = "$total_type\t$type\t$id\t$value\t$couple";
            $final{$type}{"value"} = $value;
            $final{$type}{"line"} = $line;
        }
    }
    foreach my $type (sort keys %final) {
        my $type_rc = $RC{$type};
        if ($final{$type} && $final{$type_rc}) {
            my $ratio1 = $final{$type}{"value"} / $final{$type_rc}{"value"};
            my $line = $final{$type}{"line"};
            print $fh "$line\t$ratio1\n";
        }
    }
}

sub print_con_out {
    my ($counts1, $counts2, $fh) = @_;
    my %final;

    foreach my $type (keys %$counts1) {
        my $positions = $counts1->{$type};
        foreach my $pos (sort {$a<=>$b} keys %$positions) {
            my $contexts = $positions->{$pos};
            foreach my $con (keys %$contexts) {
                my $rel1 = $$counts1{$type}{$pos}{$con}{'relative'};
                my $rel2 = $$counts2{$type}{$pos}{$con}{'relative'};
                my $abs1 = $$counts1{$type}{$pos}{$con}{'absolute'};
                my $abs2 = $$counts2{$type}{$pos}{$con}{'absolute'};

                if ($$counts1{$type}{$pos}{$con}{'relative'} && $$counts2{$type}{$pos}{$con}{'relative'}) {
                    print $fh "$id\t$type\tR1\t$rel1\t$pos\t$con\t$abs1\n";
                    print $fh "$id\t$type\tR2\t$rel2\t$pos\t$con\t$abs2\n";
                }
            }
        }
    }
}

###################################################################
######### Documentation ###########################################
###################################################################

=pod

=head1 NAME

My::Module2 -  A module with functions for DNAseq damage assessment

=head1 DESCRIPTION

=head2 C<read_counts()>

Read a file containing total mutation types output by `count_mutations.pl`.

=head3 Input

=over

=item C<$count_file> - (scalar string) Name of counts file

=back

=head3 Output

=over

=item C<%final> - (hash) final hash with counts of different types of
                  mutations, along with counts of reference bases.

=back

=head3 Example
    %final = {   'T_G' => { 'total_type' => 1,   'total' => 93 },
                 'C_C' => { 'total_type' => 33,  'total' => 33 },
                 'A_A' => { 'total_type' => 123, 'total' => 123},
                 'G_T' => { 'total_type' => 2,   'total' => 61 },
                 'G_G' => { 'total_type' => 57,  'total' => 61 },
                 'G_C' => { 'total_type' => 1,   'total' => 61 },
                 'T_T' => { 'total_type' => 92,  'total' => 93 },
                 'G_A' => { 'total_type' => 1,   'total' => 61 }};

=cut
