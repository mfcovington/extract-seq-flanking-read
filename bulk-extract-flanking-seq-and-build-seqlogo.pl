#!/usr/bin/env perl
# Mike Covington
# created: 2014-01-14
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use File::Basename;
use File::Path 'make_path';
use Getopt::Long;

my ( $ref_fa_file, $out_dir );
my $flank_length = 10;

my $options = GetOptions(
    "ref_fa_file=s"  => \$ref_fa_file,
    "out_dir=s"      => \$out_dir,
    "flank_length=i" => \$flank_length,
);

die "Specify reference FASTA file with --ref_fa_file"
    if !defined $ref_fa_file;
die "Specify output directory with --out_dir" if !defined $out_dir;

my @bam_file_list = @ARGV;

make_path $out_dir;

for my $bam_file (@bam_file_list) {
    my ( $fa_out, $logo_out ) = output_file_names( $out_dir, $bam_file );
    extract_flanking_seq( $bam_file, $ref_fa_file, $fa_out );
    build_seq_logo( $fa_out, $logo_out );    # $logo_out not yet implemented
}

exit;

sub output_file_names {
    my ( $out_dir, $bam_file ) = @_;

    my $sample = fileparse( $bam_file, ".bam" );
    return "$out_dir/$sample.fa", "$out_dir/$sample.pdf";
}

sub extract_flanking_seq {
    my ( $bam_file, $ref_fa_file, $fa_out ) = @_;

    my $extract_cmd = <<EOF;
~/git.repos/extract-seq-flanking-read/extract-seq-flanking-read.pl \\
  --bam_file $bam_file \\
  --ref_fa_file $ref_fa_file \\
  --output_fa_file $fa_out \\
  --flank_length $flank_length \\
  --fast
EOF

    system($extract_cmd);
}

sub build_seq_logo {
    my ( $fa_out, $logo_out ) = @_;    # $logo_out not yet implemented

    my $logo_cmd = <<EOF;
~/git.repos/fasta-manipulation/logo-from-fasta.pl \\
  $fa_out
EOF

    system($logo_cmd);
}
