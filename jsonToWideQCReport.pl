#!/usr/bin/perl

use strict;
use warnings;
use JSON;

my @jsonFiles = @ARGV;


my $cycleQualCut = 20;
my $cycleMismatchCut = 0.1;
my $cycleIndelCut = 0.01;
my $cycleSoftClipCut = 0.5;

my %jsonHash;

my $line;
my $file;

for $file (@jsonFiles)
{
	warn "Decoding $file\n";
	if (open (FILE, $file))
	{

	if ($line = <FILE>)
	{
		unless ($line eq "\n")
		{
			$jsonHash{$file} = decode_json($line);
		}
	}
	else
	{
		warn "No data found in $file!\n";
	}
	close FILE;
	}
	else
	{
		warn "Couldn't open $file!\n";
	}

}

my $rawReads;
my $rawYield;
my $mapRate;

my $qOver30;
my $qTotal;

my $sum;
my $count;

my $cycleCount;

my $best;
my $worst;
my $worstCycle;

my $errorRate;

my $onTargetRate;
my $estimatedYield;
my $estimatedCoverage;
my $totalYield;
my $totalCoverage;

print "File,Instrument,Run,Lane,Barcode,Library,Insert Mean,Insert Stdev,Raw Reads,Raw Yield,Uniquely Mapped %,R1 % >= q30,R2 % >= q30,R1 Mean Quality,R2 Mean Quality,Best Quality,Worst Quality,# cycles under q$cycleQualCut,R1 Mean Mismatch %,R2 Mean Mismatch %,Best Mismatch %,Worst Mismatch %,# cycles over $cycleMismatchCut % mismatch,R1 Mean Indel %, R2 Mean Indel %,Best Indel %,Worst Indel %,Worst Indel Cycle,# cycles over $cycleIndelCut % indel,R1 Mean Soft Clip %,R2 Mean Soft Clip %,Best Soft Clip %,Worst Soft Clip %,# cycles over $cycleSoftClipCut % soft clip,Total Yield,Total Coverage,Reads/SP,% on Target,Est Yield,Est Coverage,Coverage Target\n";

for my $j (sort keys %jsonHash)
{
	print "$j,";
	print $jsonHash{$j}{"instrument"} . ",";
	print $jsonHash{$j}{"run name"} . ",";
	print $jsonHash{$j}{"lane"} . ",";
	if (exists $jsonHash{$j}{"barcode"})
	{
		print $jsonHash{$j}{"barcode"} . ",";
	}
	else
	{
		print "noIndex,";
	}
	print $jsonHash{$j}{"library"} . ",";

	if ($jsonHash{$j}{"number of ends"} eq "paired end")
	{
		printf "%.2f,", $jsonHash{$j}{"insert mean"};
		printf "%.2f,", $jsonHash{$j}{"insert stdev"};
	}
	else
	{
		print "n/a,n/a,";
	}

	$rawReads = $jsonHash{$j}{"mapped reads"} + $jsonHash{$j}{"unmapped reads"} + $jsonHash{$j}{"qual fail reads"};
	$rawYield = int($rawReads * $jsonHash{$j}{"average read length"});		# cast to int because average length is only from mapped reads (and may result in ugly decimal)
	print $rawReads . ",";
	print $rawYield . ",";

	if ($rawReads > 0)
	{
		$mapRate = ($jsonHash{$j}{"mapped reads"} / $rawReads);
		print $mapRate . ",";
	}
	else
	{
		print "0,";
	}

	# quality calculations
	$qOver30 = 0;
	$qTotal = 0;
	# count first in pair and unpaired reads as read 1
	for my $q (keys %{ $jsonHash{$j}{"read 1 quality histogram"} })
	{
		if ($q >= 30)
		{
			$qOver30 += $jsonHash{$j}{"read 1 quality histogram"}{$q};
		}
		$qTotal += $jsonHash{$j}{"read 1 quality histogram"}{$q};
	}
	for my $q (keys %{ $jsonHash{$j}{"read ? quality histogram"} })
	{
		if ($q >= 30)
		{
			$qOver30 += $jsonHash{$j}{"read ? quality histogram"}{$q};
		}
		$qTotal += $jsonHash{$j}{"read ? quality histogram"}{$q};
	}
	
	if ($qTotal > 0)
	{
		print $qOver30 / $qTotal . ",";
	}
	else
	{
		print "n/a,";
	}

	$qOver30 = 0;
	$qTotal = 0;
	for my $q (keys %{ $jsonHash{$j}{"read 2 quality histogram"} })
	{
		if ($q >= 30)
		{
			$qOver30 += $jsonHash{$j}{"read 2 quality histogram"}{$q};
		}
		$qTotal += $jsonHash{$j}{"read 2 quality histogram"}{$q};
	}

	if ($qTotal > 0)
	{
		print $qOver30 / $qTotal . ",";
	}
	else
	{
		print "n/a,";
	}


	$sum = 0;
	$count = 0;

	# count first in pair and unpaired reads as read 1
	for my $q (keys %{ $jsonHash{$j}{"read 1 quality histogram"} })
	{
		$count += $jsonHash{$j}{"read 1 quality histogram"}{$q};
		$sum += $jsonHash{$j}{"read 1 quality histogram"}{$q} * $q;
	}
	for my $q (keys %{ $jsonHash{$j}{"read ? quality histogram"} })
	{
		$count += $jsonHash{$j}{"read ? quality histogram"}{$q};
		$sum += $jsonHash{$j}{"read ? quality histogram"}{$q} * $q;
	}

	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	$sum = 0;
	$count = 0;
	for my $q (keys %{ $jsonHash{$j}{"read 2 quality histogram"} })
	{
		$count += $jsonHash{$j}{"read 2 quality histogram"}{$q};
		$sum += $jsonHash{$j}{"read 2 quality histogram"}{$q} * $q;
	}

	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	
	$best = 0;
	$worst = 100;
	$cycleCount = 0;

	for my $r (qw(1 2 ?))
	{
		for my $cycle (keys %{ $jsonHash{$j}{"read $r quality by cycle"} })
		{
			if ($jsonHash{$j}{"read $r quality by cycle"}{$cycle} > $best)
			{
				$best = $jsonHash{$j}{"read $r quality by cycle"}{$cycle};
			}
			if ($jsonHash{$j}{"read $r quality by cycle"}{$cycle} < $worst)
			{
				$worst = $jsonHash{$j}{"read $r quality by cycle"}{$cycle};
			}
			if ($jsonHash{$j}{"read $r quality by cycle"}{$cycle} < $cycleQualCut)
			{
				$cycleCount++;
			}
		}
	}

	print "$best,$worst,$cycleCount,";


	# mismatch errors
	$sum = 0;
	$count = 0;

	# count first in pair and unpaired reads as read 1
	for my $cycle (keys %{ $jsonHash{$j}{"read 1 mismatch by cycle"} })
	{
		$count += $jsonHash{$j}{"read 1 aligned by cycle"}{$cycle} + $jsonHash{$j}{"read 1 insertion by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read 1 mismatch by cycle"}{$cycle};
	}
	for my $cycle (keys %{ $jsonHash{$j}{"read ? mismatch by cycle"} })
	{
		$count += $jsonHash{$j}{"read ? aligned by cycle"}{$cycle} + $jsonHash{$j}{"read ? insertion by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read ? mismatch by cycle"}{$cycle};
	}

	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	$sum = 0;
	$count = 0;
	for my $cycle (keys %{ $jsonHash{$j}{"read 2 mismatch by cycle"} })
	{
		$count += $jsonHash{$j}{"read 2 aligned by cycle"}{$cycle} + $jsonHash{$j}{"read 2 insertion by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read 2 mismatch by cycle"}{$cycle};
	}

	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	$best = 100;
	$worst = -100;
	$cycleCount = 0;

	for my $r (qw(1 2 ?))
	{
		for my $cycle (keys %{ $jsonHash{$j}{"read $r mismatch by cycle"} })
		{
			if (($jsonHash{$j}{"read $r aligned by cycle"}{$cycle} + $jsonHash{$j}{"read $r insertion by cycle"}{$cycle}) > 0)
			{
				$errorRate = $jsonHash{$j}{"read $r mismatch by cycle"}{$cycle} / ($jsonHash{$j}{"read $r aligned by cycle"}{$cycle} + $jsonHash{$j}{"read $r insertion by cycle"}{$cycle});
			}
			else
			{
				$errorRate = 0;
			}
			if ($errorRate > $worst)
			{
				$worst = $errorRate;
			}
			if ($errorRate < $best)
			{
				$best = $errorRate;
			}
			if ($errorRate > $cycleMismatchCut)
			{
				$cycleCount++;
			}
		}
	}

	print $best . ",";
	print $worst . ",";
	print $cycleCount . ",";



	# indel errors
	$sum = 0;
	$count = 0;

	# count first in pair and unpaired reads as read 1
	for my $cycle (keys %{ $jsonHash{$j}{"read 1 insertion by cycle"} })
	{
		$count += $jsonHash{$j}{"read 1 aligned by cycle"}{$cycle} + $jsonHash{$j}{"read 1 insertion by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read 1 deletion by cycle"}{$cycle} + $jsonHash{$j}{"read 1 insertion by cycle"}{$cycle};
	}
	for my $cycle (keys %{ $jsonHash{$j}{"read ? insertion by cycle"} })
	{
		$count += $jsonHash{$j}{"read ? aligned by cycle"}{$cycle} + $jsonHash{$j}{"read ? insertion by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read ? deletion by cycle"}{$cycle} + $jsonHash{$j}{"read ? insertion by cycle"}{$cycle};
	}


	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	$sum = 0;
	$count = 0;
	for my $cycle (keys %{ $jsonHash{$j}{"read 2 insertion by cycle"} })
	{
		$count += $jsonHash{$j}{"read 2 aligned by cycle"}{$cycle} + $jsonHash{$j}{"read 2 insertion by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read 2 deletion by cycle"}{$cycle} + $jsonHash{$j}{"read 2 insertion by cycle"}{$cycle};
	}

	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	$best = 100;
	$worst = -100;
	$worstCycle = -1;
	$cycleCount = 0;

	for my $r (qw(1 2 ?))
	{
		for my $cycle (keys %{ $jsonHash{$j}{"read $r insertion by cycle"} })
		{
			if (($jsonHash{$j}{"read $r aligned by cycle"}{$cycle} + $jsonHash{$j}{"read $r insertion by cycle"}{$cycle}) > 0)
			{
				$errorRate = ($jsonHash{$j}{"read $r deletion by cycle"}{$cycle} + $jsonHash{$j}{"read $r insertion by cycle"}{$cycle}) / ($jsonHash{$j}{"read $r aligned by cycle"}{$cycle} + $jsonHash{$j}{"read $r insertion by cycle"}{$cycle});
			}
			else
			{
				$errorRate = 0;
			}
			if ($errorRate > $worst)
			{
				$worst = $errorRate;
				$worstCycle = $cycle;
				if ($r eq "2")
				{
					$worstCycle += $jsonHash{$j}{"read 1 average length"};
				}
			}
			if ($errorRate < $best)
			{
				$best = $errorRate;
			}
			if ($errorRate > $cycleIndelCut)
			{
				$cycleCount++;
			}
		}
	}

	print $best . ",";
	print $worst . ",";
	print $worstCycle . ",";
	print $cycleCount . ",";

	# soft clip errors
	$sum = 0;
	$count = 0;

	# count first in pair and unpaired reads as read 1
	for my $cycle (keys %{ $jsonHash{$j}{"read 1 soft clip by cycle"} })
	{
		$count += $jsonHash{$j}{"read 1 aligned by cycle"}{$cycle} + $jsonHash{$j}{"read 1 insertion by cycle"}{$cycle} + $jsonHash{$j}{"read 1 soft clip by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read 1 soft clip by cycle"}{$cycle};
	}
	for my $cycle (keys %{ $jsonHash{$j}{"read ? soft clip by cycle"} })
	{
		$count += $jsonHash{$j}{"read ? aligned by cycle"}{$cycle} + $jsonHash{$j}{"read ? insertion by cycle"}{$cycle} + $jsonHash{$j}{"read ? soft clip by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read ? soft clip by cycle"}{$cycle};
	}

	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	$sum = 0;
	$count = 0;
	for my $cycle (keys %{ $jsonHash{$j}{"read 2 soft clip by cycle"} })
	{
		$count += $jsonHash{$j}{"read 2 aligned by cycle"}{$cycle} + $jsonHash{$j}{"read 2 insertion by cycle"}{$cycle} + $jsonHash{$j}{"read 2 soft clip by cycle"}{$cycle};
		$sum += $jsonHash{$j}{"read 2 soft clip by cycle"}{$cycle};
	}

	if ($count > 0)
	{
		print $sum / $count . ",";
	}
	else
	{
		print "n/a,";
	}

	$best = 100;
	$worst = -100;
	$cycleCount = 0;

	for my $r (qw(1 2 ?))
	{
		for my $cycle (keys %{ $jsonHash{$j}{"read $r soft clip by cycle"} })
		{
			$errorRate = $jsonHash{$j}{"read $r soft clip by cycle"}{$cycle} / ($jsonHash{$j}{"read $r aligned by cycle"}{$cycle} + $jsonHash{$j}{"read $r insertion by cycle"}{$cycle} + $jsonHash{$j}{"read $r soft clip by cycle"}{$cycle});
			if ($errorRate > $worst)
			{
				$worst = $errorRate;
			}
			if ($errorRate < $best)
			{
				$best = $errorRate;
			}
			if ($errorRate > $cycleSoftClipCut)
			{
				$cycleCount++;
			}
		}
	}

	print $best . ",";
	print $worst . ",";
	print $cycleCount . ",";


	print $jsonHash{$j}{"reads per start point"} . ",";

	if ($jsonHash{$j}{"mapped reads"} > 0)
	{
		$onTargetRate = ($jsonHash{$j}{"reads on target"} / $jsonHash{$j}{"mapped reads"});
	}
	else
	{
		$onTargetRate = "n/a,";
	}
	if (($jsonHash{$j}{"reads per start point"} > 0) and ($onTargetRate ne "n/a,"))
	{
		$estimatedYield = int(($jsonHash{$j}{"aligned bases"} * $onTargetRate) / $jsonHash{$j}{"reads per start point"});
		$totalYield = int(($jsonHash{$j}{"aligned bases"} * $onTargetRate));
	}
	else
	{
		$estimatedYield = "n/a,";
	}
	if ($estimatedYield ne "n/a,")
	{
		$estimatedCoverage = $estimatedYield / $jsonHash{$j}{"target size"};
		$totalCoverage = $totalYield / $jsonHash{$j}{"target size"};
	}
	else
	{
		$estimatedCoverage = "n/a,";
		$totalCoverage = "n/a,";
	}

	print $totalYield . ",";
	print $totalCoverage . ",";

	print $onTargetRate . ",";
	print $estimatedYield . ",";
	print $estimatedCoverage . ",";


	print $jsonHash{$j}{"target file"} . ",";


	print "\n";
}
