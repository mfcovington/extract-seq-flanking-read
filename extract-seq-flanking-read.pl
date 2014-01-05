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

get_read_info();
extract_flanking_seqs();

exit;

sub check_options {
    my $samtools_path = shift;

    die "Specify correct '--samtools_path'\n"
        unless -e $samtools_path;
}

sub get_read_info {
    # body...
}

sub extract_flanking_seqs {
    # body...
}
