#!/usr/bin/perl -w

# Author: Kelsey J.R.P. Byers (kbyers@alum.mit.edu)
# License: MIT License

# Copyright 2024 Kelsey J.R.P. Byers

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# An R script to enable quantitative data processing of GC-EAD data. See the preprint at Byers & Jacobs (2024) Quantitative analysis of gas chromatography-coupled electroantennographic detection (GC-EAD) of plant volatiles by insects. bioRxiv BIORXIV/2024/626223.

use strict;

# Check if the script was run correctly; if not, exit after providing usage
unless($ARGV[1] =~ /^[0-9]/) {
    print "gcfid-peaksabove.pl GC-FID-signal.csv GC-FID-cutoff\n";
    exit;
}

# Processing and opening input and output files for writing
my $infile = $ARGV[0];
open(FIDFILE,"$infile") || die "Couldn't open FID file $infile: $!\n|";
my $cutoff = $ARGV[1];
my $outfile = $infile . ".peaksabove." . $cutoff . ".csv";
open(OUTFILE,">$outfile") || die "Couldn't open output file $outfile: $!\n";

# Start assuming we're not in a peak
my $isapeak = 0;

# Define arrays for data output
my @onsets;
my @offsets;

# Main body of script, processing the input file line by line
while(defined(my $fidline = <FIDFILE>)) {
    my @fiddata=split(/,/,$fidline);
    # If we're already in a peak and above the cutoff, continue
    next if ($isapeak == 1 && $fiddata[1] >= $cutoff);
    # If we're already outside a peak and below the cutoff, continue
    next if ($isapeak == 0 && $fiddata[1] <= $cutoff);
    # If we're outside a peak but above the cutoff, mark a peak onset
    if ($isapeak == 0 && $fiddata[1] >= $cutoff) {
	print "At $fiddata[0], the trace became greater than $cutoff\n";
#	print OUTFILE "$fiddata[0] was the onset of a peak\n"; # for debugging
	push(@onsets,$fiddata[0]);
	$isapeak = 1;
	# If we're inside a peak but below the cutoff, mark a peak offset
    } elsif ($isapeak = 1 && $fiddata[1] <= $cutoff) {
	print "At $fiddata[0], the trace became less than $cutoff\n";
#	print OUTFILE "$fiddata[0] was the offset of a peak\n"; # for debugging
	push(@offsets,$fiddata[0]);
	$isapeak = 0;
    }
}

# Output the data to the screen and to the outfile for later use
print OUTFILE "peak starts: ",join(',',@onsets),"\n";
print OUTFILE "peak ends: ",join(',',@offsets),"\n";
print "peak starts: ",join(',',@onsets),"\n";
print "peak ends: ",join(',',@offsets),"\n";

