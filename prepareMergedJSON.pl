#!/usr/bin/perl

use strict;
use warnings;

use JSON;

my $sample = $ARGV[0];		# e.g. PCSI_0024_Ly_R
my $seqType = $ARGV[1];		# e.g. WG or EX - is there something more descriptive in the seqware meta db we could use?
my $pathToPicardMetricsFile = $ARGV[2];

my $pathToCumhistFile;
my $finalPathToImage;

if (exists $ARGV[3])
{
	$pathToCumhistFile = $ARGV[3];
}
else
{
	$pathToCumhistFile = "null";
}

if (exists $ARGV[4])
{
	$finalPathToImage = $ARGV[4];
}
else
{
	$finalPathToImage = "no image";
}

my $sampleGroup;
my $library = "merged";
my $lastModified;
my %collapsedBasesCovered;

my %picardPercentCollapsed;
my %picardEstimatedLibrarySize;

my $l;
my @f;


# use just the PCSI_0024 as the sample group
if ($sample =~  /^(.*?_.*?)_/)
{
	$sampleGroup = $1;
}
else
{
	warn "Couldn't parse sample [$sample].\n";
}


# get the last modified date of the picard metrics file (as a proxy for the bam file)
$lastModified = (stat($pathToPicardMetricsFile))[9];


# read the cumhist file into a hash
unless ($pathToCumhistFile eq "null")
{
	if (open(CUMHIST, $pathToCumhistFile))
	{
		while ($l = <CUMHIST>)
		{
			chomp $l;
			@f = split (/\t/, $l);
			$collapsedBasesCovered{$f[0]} = $f[1];
		}
		close CUMHIST;
	}
	else
	{
		warn "Couldn't open [$pathToCumhistFile].\n";
		%collapsedBasesCovered = ();
	}
}
else		# no cumhist
{
	%collapsedBasesCovered = ();
}


# read the percent collapsed and estimated library size from the picard metrics file
my $inPicardTable = 0;
if (open (PICARD, $pathToPicardMetricsFile))
{
	while ($l = <PICARD>)
	{
		chomp $l;
		@f = ();
		if (@f = split (/\t/, $l))
		{
			if ($f[0] eq "LIBRARY")
			{
				$inPicardTable = 1;
			}
			elsif ($inPicardTable == 1)
			{
				$picardPercentCollapsed{$f[0]} = $f[7];
				$picardEstimatedLibrarySize{$f[0]} = $f[8];
			}
		}
		else
		{
			$inPicardTable = 0;
		}
	}
}
else
{
	warn "Couldn't open [$pathToPicardMetricsFile].\n";
}


# build and print JSON string

my %jsonHash;

$jsonHash{"sample group"} = $sampleGroup;
$jsonHash{"sample"} = $sample;
$jsonHash{"sequencing type"} = $seqType;
$jsonHash{"library"} = $library;
$jsonHash{"last modified"} = $lastModified;
$jsonHash{"collapsed bases covered"} = \%collapsedBasesCovered;
$jsonHash{"exome histogram image"} = $finalPathToImage;
$jsonHash{"percent collapsed"} = \%picardPercentCollapsed;
$jsonHash{"estimated library size"} = \%picardEstimatedLibrarySize;

my $jsonText = encode_json(\%jsonHash);
print $jsonText;


