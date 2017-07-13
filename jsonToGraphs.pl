#!/usr/bin/perl

use strict;
use warnings;

use JSON;

my $jsonFile = $ARGV[0];		# should add usage and ability to select subset of graphs

my $insertMax = 650;		# should really pass this in
my $title;		# could pass this in too (would be different for runs and for projects)

my %jsonHash;
my $line;

open (FILE, $jsonFile) or die "Couldn't open $jsonFile.\n";

if ($line = <FILE>)
{
	%jsonHash = %{ decode_json($line) };

	unless (-d "${jsonFile}.graphs")
	{
		mkdir "${jsonFile}.graphs";
	}
}
else
{
	warn "No data found in $jsonFile!\n";
}
close FILE;

# draw graphs
my $goodColour = "forestgreen";
my $badColour = "firebrick";
my $midColour  = "goldenrod";
my @otherColours = qw(darkslategray3 darkorchid4 green2 red mediumblue orange grey yellow pink forestgreen goldenrod firebrick);

my $read1max;

my @pieArray;
my @pieLabels = ("on target", "off target", "repeat/low quality", "unmapped");
my @pieColours = ($goodColour, $midColour, $otherColours[0], $badColour);

my %tempHash;
my %histHash;
my @histArray;
my @histLabels;
my @tempArray;
my $histMax;
my @colours;

my $histHighQualCut = 30;
my $histLowQualCut = 20;

my %lineHash;
my @lineX;
my @lineY;
my $qualLineMax = 50;		# seems reasonable for phred scale

my %alignedHash;
my %insertHash;
my %softClipHash;

my $insertStep = 50;
my $insertMean;
my $insertStdev;

my %errorHash;

my $percentCovered;
my $targetSize;
my $basesCoveredDisplayMax = 200;
my $basesCoveredHighCut = 90;
my $basesCoveredLowCut = 80;

my @readArray = ("read 1", "read 2", "read ?");

@pieArray = ();
@histArray = ();
@colours = ();

if (exists $jsonHash{"barcode"})
{
	$title = $jsonHash{"run name"} . " Lane: " . $jsonHash{"lane"} . " Barcode: " . $jsonHash{"barcode"} . "\\n" . $jsonHash{"library"};
}
else
{
	$title = $jsonHash{"run name"} . " Lane: " . $jsonHash{"lane"} . "\\n" . $jsonHash{"library"};
}


push(@pieArray, $jsonHash{"reads on target"});
push(@pieArray, $jsonHash{"mapped reads"} - $jsonHash{"reads on target"});
push(@pieArray, $jsonHash{"qual fail reads"});
push(@pieArray, $jsonHash{"unmapped reads"});

pieGraph("${jsonFile}.graphs", "readPie.png", \@pieArray, \@pieLabels, \@pieColours, "$title Read Breakdown");

@histArray = ();
@histLabels = ();
@colours = ();

for my $r (@readArray)
{
	%histHash = %{ $jsonHash{"$r quality histogram"} };
	@tempArray = sort { $a <=> $b } keys %histHash;
	if ((scalar @tempArray) > 0)
	{
		$histMax = $tempArray[$#tempArray];
		if ($histMax > 0)
		{
			for (my $i = 0; $i <= $histMax; $i++)
			{
				if (exists $histHash{$i})
				{
					push(@histArray, $histHash{$i});
				}
				else
				{
					push(@histArray, 0);
				}

				push (@histLabels, $i);
		
				if ($i < $histLowQualCut)
				{
					push (@colours, $badColour);
				}
				elsif ($i < $histHighQualCut)
				{
					push (@colours, $midColour);
				}
				else
				{
					push (@colours, $goodColour);
				}
			}
		}
	}
}
barGraph("${jsonFile}.graphs", "qualHist.png", \@histArray, \@histLabels, \@colours, "$title Quality Histogram", "Base Quality", "# Bases", "");


# print non-collapsed bases covered histogram

@histArray = ();
@histLabels = ();
@colours = ();

$targetSize = $jsonHash{"target size"};
if (defined $jsonHash{"non collapsed bases covered"})
{
	%histHash = %{ $jsonHash{"non collapsed bases covered"} };
	
	if ((scalar keys %histHash) > 0)
	{
		for (my $i = 1; $i <= $basesCoveredDisplayMax; $i++)
		{
			if (exists $histHash{$i})
			{
				$percentCovered = ($histHash{$i} / $targetSize) * 100;
			}
			else
			{
				$percentCovered = 0;
			}
	
			push(@histArray, $percentCovered);
			push (@histLabels, "\"${i}x\"");
	
			if ($percentCovered < $basesCoveredLowCut)
			{
				push (@colours, $badColour);
			}
			elsif ($percentCovered < $basesCoveredHighCut)
			{
				push (@colours, $midColour);
			}
			else
			{
				push (@colours, $goodColour);
			}
		}
	}
	barGraph("${jsonFile}.graphs", "nonCollapsedCover.png", \@histArray, \@histLabels, \@colours, "$title Non Collapsed Coverage", "Coverage Depth", "% Target Covered", ", ylim=c(0, 100)");
}


# print collapsed bases covered histogram

@histArray = ();
@histLabels = ();
@colours = ();

$targetSize = $jsonHash{"target size"};
if (defined $jsonHash{"collapsed bases covered"})
{
	%histHash = %{ $jsonHash{"collapsed bases covered"} };
	
	if ((scalar keys %histHash) > 0)
	{
		for (my $i = 1; $i <= $basesCoveredDisplayMax; $i++)
		{
			if (exists $histHash{$i})
			{
				$percentCovered = ($histHash{$i} / $targetSize) * 100;
			}
			else
			{
				$percentCovered = 0;
			}
	
			push(@histArray, $percentCovered);
			push (@histLabels, "\"${i}x\"");
	
			if ($percentCovered < $basesCoveredLowCut)
			{
				push (@colours, $badColour);
			}
			elsif ($percentCovered < $basesCoveredHighCut)
			{
				push (@colours, $midColour);
			}
			else
			{
				push (@colours, $goodColour);
			}
		}
	}
	barGraph("${jsonFile}.graphs", "collapsedCover.png", \@histArray, \@histLabels, \@colours, "$title Collapsed Coverage", "Coverage Depth", "% Target Covered", ", ylim=c(0, 100)");
}


# read length histogram

@lineX = ();
@lineY = ();
@colours = ();
%lineHash = ();
for my $r (@readArray)
{
	if (defined $jsonHash{"$r length histogram"})
	{
		for my $l (keys  %{ $jsonHash{"$r length histogram"} })
		{
			if (exists $lineHash{$l})
			{
				$lineHash{$l} += $jsonHash{"$r length histogram"}{$l};
			}
			else
			{
				$lineHash{$l} = $jsonHash{"$r length histogram"}{$l};
			}
		}
	}

	for my $l (sort {$a <=> $b} keys %lineHash)
	{
		push (@lineX, $l);
		push (@lineY, $lineHash{$l});
		push (@colours, $goodColour);
	}
}
lineGraph("${jsonFile}.graphs", "readLength.png", \@lineX, \@lineY, \@colours, "$title Read Length Histogram", "Read Length (bp)", "Number of reads", "");



unless ($jsonHash{"number of ends"} eq "single end")
{
	# insert graph

	@lineX = ();
	@lineY = ();
	@colours = ();
	$insertMean =  $jsonHash{"insert mean"};
	%lineHash = %{ $jsonHash{"insert histogram"} };
	for my $i (sort {$a <=> $b} keys %lineHash)
	{
		if ($i < $insertMax)
		{
			push(@lineX, $i);
			push(@lineY, $lineHash{$i});

			if (($i < ($insertMean - (2 * $insertStep))) or ($i > ($insertMean + (2 * $insertStep))))
			{
				push(@colours, $badColour);
			}
			elsif (($i < ($insertMean - $insertStep)) or ($i > ($insertMean + $insertStep)))
			{
				push(@colours, $midColour);
			}
				else
			{
				push(@colours, $goodColour);
			}
		}
	}
	lineGraph("${jsonFile}.graphs", "insert.png", \@lineX, \@lineY, \@colours, "$title Insert Distribution", "Insert Size (bp)", "Pairs", "");
}

# quality by cycle graph

@lineX = ();
@lineY = ();
@colours = ();
$read1max = 0;
for my $r (@readArray)
{
	if (defined $jsonHash{"$r quality by cycle"})
	{
		%lineHash = %{ $jsonHash{"$r quality by cycle"} };
		for my $i (sort {$a <=> $b} keys %lineHash)
		{
			if ($jsonHash{"number of ends"} eq "single end")
			{
				push(@lineX, $i);
				push(@lineY, $lineHash{$i});
			}
			else
			{
				if ($r eq "read 1")
				{
					push(@lineX, $i);
					push(@lineY, $lineHash{$i});
					$read1max++;
				}
				if ($r eq "read 2")
				{
					push(@lineX, $i + $read1max);
					push(@lineY, $lineHash{$i});
				}
			}

			if ($lineHash{$i} < $histLowQualCut)
			{
				push (@colours, $badColour);
			}
			elsif ($lineHash{$i} < $histHighQualCut)
			{
				push (@colours, $midColour);
			}
			else
			{
				push (@colours, $goodColour);
			}
		}
	}
	if ($r eq "read 1")
	{
		push(@lineX, $read1max);
		push(@lineY, 0);
		push(@colours, "white");
		$read1max++;
	}
}
lineGraph("${jsonFile}.graphs", "qualCycle.png", \@lineX, \@lineY, \@colours, "$title Quality by Cycle", "Cycle", "Mean Quality", ", ylim=c(0, $qualLineMax)");


for my $r (@readArray)
{
	$alignedHash{$r} = \%{ $jsonHash{"$r aligned by cycle"} };
	$insertHash{$r} = \%{ $jsonHash{"$r insertion by cycle"} };
	$softClipHash{$r} = \%{ $jsonHash{"$r soft clip by cycle"} };

}

# mismatch by cycle

@lineX = ();
@lineY = ();
$read1max = 0;
@colours = ();
for my $r (@readArray)
{
	if ((scalar keys %{ $jsonHash{"$r mismatch by cycle"} }) > 0)
	{
		%errorHash = %{ $jsonHash{"$r mismatch by cycle"} };
		for my $i (sort {$a <=> $b} keys %errorHash)
		{
			if ($jsonHash{"number of ends"} eq "single end")
			{
				push(@lineX, $i);
			}
			else
			{
				if ($r eq "read 1")
				{
					push(@lineX, $i);
					$read1max++;
				}
				elsif ($r eq "read 2")
				{
					push(@lineX, $i + $read1max);
				}
			}
			if (($alignedHash{$r}{$i} + $insertHash{$r}{$i}) > 0)
			{
				push(@lineY, (($errorHash{$i} / ($alignedHash{$r}{$i} + $insertHash{$r}{$i})) * 100));
			}
			else
			{
				push(@lineY, 0);
			}
			push(@colours, $badColour);
		}
		if ($r eq "read 1")
		{
			push(@lineX, $read1max);
			push(@lineY, 0);
			push(@colours, "white");
			$read1max++;
		}
	}
}
lineGraph("${jsonFile}.graphs", "misCycle.png", \@lineX, \@lineY, \@colours, "$title Mismatches by Cycle", "Cycle", "% Bases Mismatched", ", ylim=c(0, 10)");

# indel by cycle

@lineX = ();
@lineY = ();
@colours = ();
$read1max = 0;
for my $r (@readArray)
{
	if ((scalar keys %{ $jsonHash{"$r deletion by cycle"} }) > 0)
	{
		%errorHash = %{ $jsonHash{"$r deletion by cycle"} };
		for my $i (sort {$a <=> $b} keys %errorHash)
		{
			if ($jsonHash{"number of ends"} eq "single end")
			{
				push(@lineX, $i);
			}
			else
			{
				if ($r eq "read 1")
				{
					push(@lineX, $i);
					$read1max++;
				}
				elsif ($r eq "read 2")
				{
					push(@lineX, $i + $read1max);
				}
			}
			if (($alignedHash{$r}{$i} + $insertHash{$r}{$i}) > 0)
			{
				push(@lineY, ((($errorHash{$i} + $insertHash{$r}{$i}) / ($alignedHash{$r}{$i} + $insertHash{$r}{$i})) * 100));
			}
			else
			{
				push(@lineY, 0);
			}
			push(@colours, $badColour);
		}
		if ($r eq "read 1")
		{
			push(@lineX, $read1max);
			push(@lineY, 0);
			push(@colours, "white");
			$read1max++;
		}
	}
}
lineGraph("${jsonFile}.graphs", "indelCycle.png", \@lineX, \@lineY, \@colours, "$title Indels by Cycle", "Cycle", "% Bases Inserted/Deleted", ", ylim=c(0, 10)");

# soft clip by cycle

@lineX = ();
@lineY = ();
@colours = ();
$read1max = 0;
for my $r (@readArray)
{
	if ((scalar keys %{ $jsonHash{"$r soft clip by cycle"} }) > 0)
	{
		%errorHash = %{ $jsonHash{"$r soft clip by cycle"} };
		for my $i (sort {$a <=> $b} keys %errorHash)
		{
			if ($jsonHash{"number of ends"} eq "single end")
			{
				push(@lineX, $i);
			}
			else
			{
				if ($r eq "read 1")
				{
					push(@lineX, $i);
					$read1max++;
				}
				elsif ($r eq "read 2")
				{
					push(@lineX, $i + $read1max);
				}
			}
			if (($alignedHash{$r}{$i} + $insertHash{$r}{$i} + $errorHash{$i}) > 0)
			{
				push(@lineY, (($errorHash{$i} / ($alignedHash{$r}{$i} + $errorHash{$i} + $insertHash{$r}{$i})) * 100));
			}
			else
			{
				push(@lineY, 0);
			}
			push(@colours, $badColour);
		}
		if ($r eq "read 1")
		{
			push(@lineX, $read1max);
			push(@lineY, 0);
			push(@colours, "white");
			$read1max++;
		}
	}
}
lineGraph("${jsonFile}.graphs", "softCycle.png", \@lineX, \@lineY, \@colours, "$title Soft Clips by Cycle", "Cycle", "% Bases Soft Clipped", ", ylim=c(0, 100)");


# hard clip by cycle

@lineX = ();
@lineY = ();
@colours = ();
$read1max = 0;
for my $r (@readArray)
{
	if ((scalar keys %{ $jsonHash{"$r hard clip by cycle"} }) > 0)
	{
		%errorHash = %{ $jsonHash{"$r hard clip by cycle"} };
		for my $i (sort {$a <=> $b} keys %errorHash)
		{
			if ($jsonHash{"number of ends"} eq "single end")
			{
				push(@lineX, $i);
			}
			else
			{
				if ($r eq "read 1")
				{
					push(@lineX, $i);
					$read1max++;
				}
				elsif ($r eq "read 2")
				{
					push(@lineX, $i + $read1max);
				}
			}
			if (($alignedHash{$r}{$i} + $insertHash{$r}{$i} + $errorHash{$i}) > 0)
			{
				push(@lineY, (($errorHash{$i} / ($alignedHash{$r}{$i} + $errorHash{$i} + $insertHash{$r}{$i} + $softClipHash{$r}{$i})) * 100));
			}
			else
			{
				push(@lineY, 0);
			}
			push(@colours, $badColour);
		}
		if ($r eq "read 1")
		{
			push(@lineX, $read1max);
			push(@lineY, 0);
			push(@colours, "white");
			$read1max++;
		}
	}
}
lineGraph("${jsonFile}.graphs", "hardCycle.png", \@lineX, \@lineY, \@colours, "$title Hard Clips by Cycle", "Cycle", "% Bases Hard Clipped", ", ylim=c(0, 100)");


sub pieGraph
{
	my $path = $_[0];
	my $name = $_[1];
	my @values = @{ $_[2] };
	my @labels = @{ $_[3] };
	my @colours = @{ $_[4] };
	my $title = $_[5];

	open (RFILE, ">${path}/${name}.Rcode") or die "Couldn't open ${path}/${name}.Rcode.\n";

	print RFILE "slices <- c($values[0]";
	for (my $i = 1; $i < scalar @values; $i++)
	{
		print RFILE ", $values[$i]";
	}
	print RFILE ")\n";

	print RFILE "lbls <- c(\"$labels[0]\"";
	for (my $i = 1; $i < scalar @labels; $i++)
	{
		print RFILE ", \"$labels[$i]\"";
	}
	print RFILE ")\n";

	print RFILE "cols <- c(\"$colours[0]\"";
	for (my $i = 1; $i < scalar @colours; $i++)
	{
		print RFILE ", \"$colours[$i]\"";
	}
	print RFILE ")\n";

	print RFILE "pct <- round(slices/sum(slices)*100)\n";
	print RFILE "lbls <- paste(lbls, pct)\n";
	print RFILE "lbls <- paste(lbls,\"%\",sep=\"\")\n";

	print RFILE "png(filename = \"${path}/${name}\", width = 640, height = 640)\n";
	print RFILE "pie(slices, labels = lbls, col=cols, main=\"$title\", border=NA)\n";
	print RFILE "dev.off()\n";

	close RFILE;

	`Rscript ${path}/${name}.Rcode`;
}
sub lineGraph
{
	my $path = $_[0];
	my $name = $_[1];
	my @xVal = @{ $_[2] };
	my @yVal = @{ $_[3] };
	my @colours = @{ $_[4] };
	my $title = $_[5];
	my $xlab = $_[6];
	my $ylab = $_[7];
	my $additionalParams = $_[8];
	
	open (RFILE, ">${path}/${name}.Rcode") or die "Couldn't open ${path}/${name}.Rcode.\n";

	print RFILE "xvals <- c($xVal[0]";
	for (my $i = 1; $i < scalar @xVal; $i++)
	{
		print RFILE ", $xVal[$i]";
	}
	print RFILE ")\n";

	print RFILE "yvals <- c($yVal[0]";
	for (my $i = 1; $i < scalar @yVal; $i++)
	{
		print RFILE ", $yVal[$i]";
	}
	print RFILE ")\n";

	print RFILE "cols <- c(\"$colours[0]\"";
	for (my $i = 1; $i < scalar @colours; $i++)
	{
		print RFILE ", \"$colours[$i]\"";
	}
	print RFILE ")\n";

	print RFILE "png(filename = \"${path}/${name}\", width = 640, height = 640)\n";
	print RFILE "plot(xvals, yvals, main=\"$title\", type=\"n\", col=\"black\", xlab=\"$xlab\", ylab=\"$ylab\"$additionalParams)\n";
	print RFILE "for(i in 1:(length(yvals)-1))\n{\n";
#	print RFILE "polygon(c(xvals[i], xvals[i], mean(c(xvals[i], xvals[i+1])), mean(c(xvals[i], xvals[i+1]))), c(0, yvals[i], yvals[i], 0), col=cols[i], border=NA)\n";   # c(0, yvals[i], mean(c(yvals[i], yvals[i+1])), 0), col=cols[i], border=NA)\n";
#	print RFILE "polygon(c(mean(c(xvals[i], xvals[i+1])), mean(c(xvals[i], xvals[i+1])), xvals[i+1], xvals[i+1]), c(0, yvals[i+1], yvals[i+1], 0), col=cols[i+1], border=NA)\n";  # c(0, mean(c(yvals[i], yvals[i+1])), yvals[i+1], 0), col=cols[i+1], border=NA)\n";
	print RFILE "polygon(c(xvals[i] - 0.5, xvals[i] - 0.5, xvals[i] + 0.5, xvals[i] + 0.5), c(0, yvals[i], yvals[i], 0), col=cols[i], border=NA)\n";   # c(0, yvals[i], mean(c(yvals[i], yvals[i+1])), 0), col=cols[i], border=NA)\n";
	print RFILE "}\n";
	print RFILE "dev.off()\n";

	close RFILE;
	`Rscript ${path}/${name}.Rcode`;
}
sub barGraph
{
	my $path = $_[0];
	my $name = $_[1];
	my @values = @{ $_[2] };
	my @labels = @{ $_[3] };
	my @colours = @{ $_[4] };
	my $title = $_[5];
	my $xlab = $_[6];
	my $ylab = $_[7];
	my $additionalParams = $_[8];
	
	open (RFILE, ">${path}/${name}.Rcode") or die "Couldn't open ${path}/${name}.Rcode.\n";

	print RFILE "values <- c($values[0]";
	for (my $i = 1; $i < scalar @values; $i++)
	{
		print RFILE ", $values[$i]";
	}
	print RFILE ")\n";

	print RFILE "labels <- c($labels[0]";
	for (my $i = 1; $i < scalar @labels; $i++)
	{
		print RFILE ", $labels[$i]";
	}
	print RFILE ")\n";

	print RFILE "cols <- c(\"$colours[0]\"";
	for (my $i = 1; $i < scalar @colours; $i++)
	{
		print RFILE ", \"$colours[$i]\"";
	}
	print RFILE ")\n";

	print RFILE "png(filename = \"${path}/${name}\", width = 640, height = 640)\n";
	print RFILE "barplot(values, col=cols, names.arg=labels, main=\"$title\", xlab=\"$xlab\", ylab=\"$ylab\", border=cols$additionalParams)\n";
	print RFILE "dev.off()\n";

	close RFILE;
	`Rscript ${path}/${name}.Rcode`;
}
