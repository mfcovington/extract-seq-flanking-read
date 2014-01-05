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
use Data::Printer; # TEMP

#TODO: Add README

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
extract_flanking_seqs( $read_stats, $flank_length, $samtools_path );
write_to_fasta( $read_stats, $output_fa_file, $fa_width );

p $read_stats;

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
        $read_stats{$read_id}{strand} = $flag =~ /r/ ? "rev" : "fwd";
        $read_stats{$read_id}{pos}    = $pos;
        $read_stats{$read_id}{cigar}  = $cigar;
    }
    close $bam_fh;

    return \%read_stats;
}

sub extract_flanking_seqs {
    my ( $read_stats, $flank_length, $samtools_path ) = @_;


}

sub write_to_fasta {
    my ( $read_stats, $output_fa_file, $fa_width ) = @_;


}
