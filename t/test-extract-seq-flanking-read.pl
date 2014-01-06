#!/usr/bin/env perl
# Mike Covington
# created: 2014-01-04
#
# Description:
#
use strict;
use warnings;
use autodie;
use Test::More tests => 1;

my $base_extract_cmd = <<CMD;
../extract-seq-flanking-read.pl \\
  --bam_file sample-files/sample.bam \\
  --ref_fa_file sample-files/sample.fa \\
  --output_fa_file got.fa \\
CMD
my $test_name;
my $extract_cmd;

$test_name = "10bp-upstream";
$extract_cmd = "$base_extract_cmd";
compare_extracted_seq( $extract_cmd, $test_name );

sub compare_extracted_seq {
    my ( $extract_cmd, $test_name ) = @_;

    system($extract_cmd);

    my $expect_file = "sample-files/expect.$test_name.fa";
    open my $expect_fh, "<", $expect_file;
    my @expected = <$expect_fh>;
    close $expect_fh;

    open my $got_fh, "<", "got.fa";
    my @got = <$got_fh>;
    close $got_fh;

    is_deeply( \@got, \@expected, $test_name );

    unlink "got.fa";
}
