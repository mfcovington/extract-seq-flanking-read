#!/usr/bin/env perl
# Mike Covington
# created: 2014-01-04
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use List::Util 'sum';

#TODO: Add README
#TODO: Consider other CIGAR score variants?

# Defaults
my $bam_file       = "t/sample-files/sample.bam";
my $ref_fa_file    = "t/sample-files/sample.fa";
my $output_fa_file = "out.fa";
my $samtools_path  = glob "~/installs/bin/samtools";

my $flank_length = 10;
my $fa_width     = 80;

my $options = GetOptions(
    "bam_file=s"       => \$bam_file,
    "ref_fa_file=s"    => \$ref_fa_file,
    "output_fa_file=s" => \$output_fa_file,
    "samtools_path=s"  => \$samtools_path,
    "flank_length=i"   => \$flank_length,
    "fa_width=i"       => \$fa_width,
);

check_options( $samtools_path );

my $read_stats = get_read_info( $bam_file, $samtools_path );
extract_flanking_seqs( $read_stats, $flank_length, $ref_fa_file,
    $samtools_path );
write_to_fasta( $read_stats, $output_fa_file, $fa_width );

exit;

sub check_options {
    my $samtools_path = shift;

    die "Specify correct '--samtools_path'\n"
        unless -e $samtools_path;
}

sub get_read_info {
    my ( $bam_file, $samtools_path ) = @_;

    my %read_stats;

    my $samtools_cmd = "$samtools_path view -X $bam_file";
    open my $bam_fh, "-|", $samtools_cmd;
    while ( my $read = <$bam_fh> ) {
        my ( $read_id, $flag, $seq_id, $pos, $cigar )
            = ( split /\t/, $read )[ 0 .. 3, 5 ];

        next if $flag =~ /u/;

        $read_stats{$read_id}{seq_id} = $seq_id;
        $read_stats{$read_id}{strand} = $flag =~ /r/ ? "-" : "+";
        $read_stats{$read_id}{pos}    = $pos;
        $read_stats{$read_id}{cigar}  = $cigar;
    }
    close $bam_fh;

    return \%read_stats;
}

sub extract_flanking_seqs {
    my ( $read_stats, $flank_length, $ref_fa_file, $samtools_path ) = @_;

    get_positions( $read_stats, $flank_length );
    get_sequences( $read_stats, $ref_fa_file, $samtools_path );
}

sub get_positions {
    my ( $read_stats, $flank_length ) = @_;

    for my $read_id ( keys $read_stats ) {
        my $seq_id = $$read_stats{$read_id}{seq_id};
        my $strand = $$read_stats{$read_id}{strand};
        my $pos    = $$read_stats{$read_id}{pos};
        my $cigar  = $$read_stats{$read_id}{cigar};

        my ( $start, $end );
        if ( $strand eq "+" ) {
            $start = $pos - $flank_length;
            $end   = $pos - 1;
        }
        elsif ( $strand eq "-" ) {
            my $rt_pos = $pos - 1 + cigar_to_length($cigar);
            $start = $rt_pos + 1;
            $end   = $rt_pos + $flank_length;
        }
        else { die "Problem with strand info\n" }

        $$read_stats{$read_id}{start} = $start;
        $$read_stats{$read_id}{end}   = $end;
    }
}

sub cigar_to_length {
    my $cigar = shift;

    my @add = $cigar =~ /(\d+)[MND]/g;
    my $length = sum @add;

    return $length;
}

sub get_sequences {
    my ( $read_stats, $ref_fa_file, $samtools_path ) = @_;

    for my $read_id ( keys $read_stats ) {
        my $seq_id = $$read_stats{$read_id}{seq_id};
        my $strand = $$read_stats{$read_id}{strand};
        my $start  = $$read_stats{$read_id}{start};
        my $end    = $$read_stats{$read_id}{end};

        $$read_stats{$read_id}{flank}
            = extract_fa_seq( $samtools_path, $ref_fa_file, $seq_id, $strand,
            $start, $end );
    }
}

sub extract_fa_seq {    # This subroutine from extract-utr v0.2.1
    my ( $samtools_path, $fa_file, $seqid, $strand, $left_pos, $right_pos )
        = @_;

    my $faidx_cmd =
      defined $left_pos && defined $right_pos
      ? "$samtools_path faidx $fa_file $seqid:$left_pos-$right_pos"
      : "$samtools_path faidx $fa_file $seqid";

    my ( $fa_header, @fa_seq ) = `$faidx_cmd`;
    chomp @fa_seq;

    my $seq = join "", @fa_seq;

    if ( defined $strand && $strand eq '-' ) {
        $seq = reverse $seq;
        $seq =~ tr/ACGTacgt/TGCAtgca/;
    }

    return $seq;
}

sub write_to_fasta {
    my ( $read_stats, $output_fa_file, $fa_width ) = @_;

    my @ids_sorted_by_coords = sort {
        $$read_stats{$a}{seq_id} cmp $$read_stats{$b}{seq_id}
            || $$read_stats{$a}{start} <=> $$read_stats{$b}{start}
    } keys $read_stats;

    open my $output_fa_fh, ">", $output_fa_file;
    for my $read_id ( @ids_sorted_by_coords ) {
        my $seq_id = $$read_stats{$read_id}{seq_id};
        my $flank  = $$read_stats{$read_id}{flank};

        output_fa( $read_id, $flank, $output_fa_fh, $fa_width );
    }
    close $output_fa_fh;
}

sub output_fa {    # This subroutine from extract-utr v0.2.1
    my ( $seqid, $seq, $output_fa_fh, $fa_width ) = @_;

    $fa_width //= 80;
    my @fa_seq = unpack "(A$fa_width)*", $seq;

    say $output_fa_fh ">$seqid";
    say $output_fa_fh $_ for @fa_seq;
}

