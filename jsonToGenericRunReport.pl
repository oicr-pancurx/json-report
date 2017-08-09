#!/usr/bin/perl

use strict;
use warnings;

use JSON;
use Cwd 'abs_path';
use Getopt::Std;
use File::Basename;
use vars qw/ %opt /;

my $scriptPath = dirname(abs_path($0));

sub usage
{
        print "\nUsage is jsonToGenericRunReport.pl [options] path/to/*.json\n";
        print "Options are as follows:\n";
		print "\t -c show bases covered data.  Default is no bases covered.\n";
        print "\t-r show read length histogram.  Default is no histogram.\n";
        print "\t-p print all images.  Default is to only show thumbnails (with links).\n";
		print "\t-H show hard clip stats and graph.  Default is off.\n";
        print "\t-h displays this usage message.\n";

        die "\n@_\n\n";
}


my $showReadLengthHist = 0;
my $printAllImages = 0;
my $showHardClip = 0;
my $showBasesCovered = 0;

my $opt_string = "crpHh";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
	usage("Help requested.");
}
if (exists $opt{r})
{
	$showReadLengthHist = 1;
}
if (exists $opt{c})
{
	$showBasesCovered = 1;
}
if (exists $opt{p})
{
	$printAllImages = 1;
}
if (exists $opt{H})
{
	$showHardClip = 1;
}


my @jsonFiles = @ARGV;


my %jsonHash;
my %runList;

my $line;


my $ends = "single end";

for my $file (@jsonFiles)
{
	open (FILE, $file) or die "Couldn't open $file.\n";

	if ($line = <FILE>)
	{
		$jsonHash{$file} = decode_json($line);
		$runList{$jsonHash{$file}{"run name"}}{$jsonHash{$file}{"lane"}}++;

		unless (-d "${file}.graphs")
		{
			mkdir "${file}.graphs";
		}
		if ($jsonHash{$file}{"number of ends"} eq "paired end")
		{
			$ends = "paired end";
		}
	}
	else
	{
		warn "No data found in $file!\n";
	}
}

my $date = `date`;
chomp $date;

my $rawReads;
my $rawYield;
my $insertMean;
my $insertStdev;
my $mapRate;
my $errorRate;
my $softClipRate;
my $hardClipRate;
my $readsPerStartPoint;
my $onTargetRate;
my $estimatedYield;
my $estimatedCoverage;

my %qualCuts;
my $qualCut;

my $totalReads;
my $totalRawYield;
my $totalEstYield;

my $readLength;

my $targetSize;
my $numberOfTargets;
my $basesCovered;
my @coverageXs = qw(1 4 8 15 30 50 100 200);

my $linkName;

unless (-e "sorttable.js")
{
	`ln -s /.mounts/labs/PCSI/production/phoenix-report/sorttable.js`;
}

print "<html>\n<head>\n";
print "<script src=\"./sorttable.js\"></script>\n";
print "<style type=\"text/css\">\n.na { color: #ccc; }\nth, td {\n  padding: 3px !important;\n}\ntable\n{\nborder-collapse:collapse;\n}\n/* Sortable tables */\ntable.sortable thead {\n\tbackground-color:#eee;\n\tcolor:#000000;\n\tfont-weight: bold;\n\tcursor: default;\n}\n</style>\n";
print "</head>\n<body>\n";
print "<p>Generic run report generated on $date.</p>\n";

for my $run (sort keys %runList)
{
	$totalRawYield = 0;
	$totalReads = 0;
	$totalEstYield = 0;

	print "<h1><a name=\"$run\">$run</a></h1>\n";

	print "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";

	print "<th>Lane</th>";
	print "<th>Barcode</th>";
	print "<th>Group ID</th>";
    print "<th>External Name</th>";
	print "<th>Library</th>";

	if ($ends eq "paired end")
	{
		print "<th>Insert Mean (SD)</th>";
	}

	print "<th>Read Length</th>";
	print "<th>Raw Reads</th>";
	print "<th>Raw Yield</th>";
	print "<th>Map %</th>";
	print "<th>Error %</th>";
	print "<th>Soft Clip %</th>";

	if ($showHardClip == 1)
	{
		print "<th>Hard Clip %</th>";
	}
	print "<th>Reads/SP</th>";
	print "<th>% on Target</th>";
	print "<th>Estimated Yield*</th>";
	print "<th>Coverage*</th>";

	print "\n</thead>\n<tbody>\n";

	for my $lane (sort keys %{ $runList{$run} })
	{
		for my $j (keys %jsonHash)
		{
			if (($jsonHash{$j}{"run name"} eq $run) and ($jsonHash{$j}{"lane"} eq $lane))
			{
				$rawReads = $jsonHash{$j}{"mapped reads"} + $jsonHash{$j}{"unmapped reads"} + $jsonHash{$j}{"qual fail reads"};

				$rawYield = int($rawReads * $jsonHash{$j}{"average read length"});		# cast to int because average length is only from mapped reads (and may result in ugly decimal)

				if ($ends eq "paired end")
				{
					$insertMean = sprintf "%.2f", $jsonHash{$j}{"insert mean"};
					$insertStdev =  sprintf "%.2f", $jsonHash{$j}{"insert stdev"};
				}

				if ($jsonHash{$j}{"number of ends"} eq "paired end")
				{
					$readLength = $jsonHash{$j}{"read 1 average length"} . ", " . $jsonHash{$j}{"read 2 average length"};
				}
				else
				{
					$readLength = sprintf "%.2f", $jsonHash{$j}{"read ? average length"};
				}

				$mapRate = "0%";
				if ($rawReads > 0)
				{
					$mapRate = ($jsonHash{$j}{"mapped reads"} / $rawReads) * 100;
					$mapRate = sprintf "%.2f%%", $mapRate;
				}

				$errorRate = "0%";
				$softClipRate = "0%";
				$hardClipRate = "0%";
				if ($jsonHash{$j}{"aligned bases"} > 0)
				{
					$errorRate = (($jsonHash{$j}{"mismatch bases"} + $jsonHash{$j}{"inserted bases"} + $jsonHash{$j}{"deleted bases"}) / $jsonHash{$j}{"aligned bases"}) * 100;
					$errorRate = sprintf "%.2f%%", $errorRate;

					$softClipRate = $jsonHash{$j}{"soft clip bases"} / ($jsonHash{$j}{"aligned bases"} + $jsonHash{$j}{"soft clip bases"} + $jsonHash{$j}{"hard clip bases"}) * 100;
					$softClipRate = sprintf "%.2f%%", $softClipRate;

					$hardClipRate = $jsonHash{$j}{"hard clip bases"} / ($jsonHash{$j}{"aligned bases"} + $jsonHash{$j}{"soft clip bases"} + $jsonHash{$j}{"hard clip bases"}) * 100;
					$hardClipRate = sprintf "%.2f%%", $hardClipRate;
				}

				$readsPerStartPoint = sprintf "%.2f", $jsonHash{$j}{"reads per start point"};

				$onTargetRate = 0;
				if ($rawReads > 0)
				{
					$onTargetRate = ($jsonHash{$j}{"reads on target"} / $jsonHash{$j}{"mapped reads"}) * 100;  # $rawReads) * 100;   # could argue using this either way
				}

				$estimatedYield = int(($jsonHash{$j}{"aligned bases"} * ($onTargetRate / 100)) / $jsonHash{$j}{"reads per start point"});
				
				$estimatedCoverage = $estimatedYield / $jsonHash{$j}{"target size"};
				$estimatedCoverage = sprintf "%.2f", $estimatedCoverage;

				$totalRawYield += $rawYield;
				$totalReads += $rawReads;
				$totalEstYield += $estimatedYield;

				$onTargetRate = sprintf "%.2f%%", $onTargetRate;
				$rawReads =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
				$rawYield =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
				$estimatedYield =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;

				$qualCuts{$jsonHash{$j}{"qual cut"}} = 1;

				if (exists $jsonHash{$j}{barcode})
				{
					$linkName = "${run}_${lane}_$jsonHash{$j}{barcode}";
				}
				else
				{
					$linkName = "${run}_${lane}";
				}

				print "<tr>\n";
				print "<td><a href=\"#$linkName\">$jsonHash{$j}{lane}</a></td>";
				if (exists $jsonHash{$j}{barcode})
				{
					print "<td><a href=\"#$linkName\">$jsonHash{$j}{barcode}</a></td>";
				}
				else
				{
					print "<td>none</td>";
				}
				
	            if (exists $jsonHash{$j}{"group id"})
                {
                    my $group = "<td";
                    if (exists $jsonHash{$j}{"group id description"}) {
                          $group .= " title='" . $jsonHash{$j}{"group id description"} . "'";
                    }
                    print $group . ">$jsonHash{$j}{'group id'}</td>";
                }
                else
                {
                    print "<td class='na'>na</td>";
                }
                if (exists $jsonHash{$j}{"external name"})
                {
                    print "<td>$jsonHash{$j}{'external name'}</td>";
                }
                else
                {
                    print "<td class='na'>na</td>";
                }
				
				print "<td><a href=\"#$linkName\">$jsonHash{$j}{library}</a></td>";
				if ($ends eq "paired end")
				{
					print "<td>$insertMean ($insertStdev)</td>";
				}
				print "<td>$readLength</td>";
				print "<td>$rawReads</td>";
				print "<td>$rawYield</td>";
				print "<td>$mapRate</td>";
				print "<td>$errorRate</td>";
				print "<td>$softClipRate</td>";
				if ($showHardClip == 1)
				{
					print "<td>$hardClipRate</td>";
				}
				print "<td>$readsPerStartPoint</td>";
				print "<td>$onTargetRate</td>";
				print "<td>$estimatedYield</td>";
				print "<td>${estimatedCoverage}x</td>";
				print "\n</tr>\n";

			}
		}
	}

	$totalRawYield =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
	$totalReads =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
	$totalEstYield =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;

	print "</tbody>\n";
	print "<tfoot>\n<tr>\n";
	print "<td>Total</td>";
	if ($ends eq "paired end")
	{
		print "<td></td>";
	}
	print "<td></td><td></td><td></td><td></td><td></td><td>$totalReads</td><td>$totalRawYield</td><td></td><td></td>";
	if ($showHardClip == 1)
	{
		print "<td></td>";
	}
	print "<td></td><td></td><td></td><td>$totalEstYield</td><td></td>\n";
	print "</tr>\n";
	print "</table>\n";


	# handle multiple qual cuts
	my @qualArray;
	if (scalar(keys %qualCuts) == 1)
	{
		@qualArray = keys %qualCuts;
		$qualCut = $qualArray[0];
	}
	else
	{
		$qualCut = "";
		for my $q (keys %qualCuts)
		{
			if ($qualCut eq "")
			{
				$qualCut = $q;
			}
			else
			{
				$qualCut = "$qualCut, $q";
			}
		}
	}

	print "<p>* Estimates exclude unmapped, off target, non-primary or MAPQ < $qualCut reads and use reads per start point to approximate loss due to collapse.</p>\n";


	# print coverage table
	if ($showBasesCovered == 1)
	{
		print "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";

		print "<th colspan=\"5\"></th>";
		print "<th colspan=\"" . scalar @coverageXs . "\">Non-collapsed % bases covered</th>";
		print "<th colspan=\"" . scalar @coverageXs . "\">Collapsed % bases covered</th>";

		print "</tr><tr>";

		print "<th>Lane</th>";
		print "<th>Barcode</th>";
		print "<th>Library</th>";
		print "<th>Target Size (bp)</th>";
		print "<th># Targets</th>";

		for my $i (@coverageXs)
		{
			print "<th>${i}x</th>";
		}
		for my $i (@coverageXs)
		{
			print "<th>${i}x</th>";
		}

		print "\n</tr>\n</thead>\n<tbody>\n";

		for my $lane (sort keys %{ $runList{$run} })
		{
			for my $j (keys %jsonHash)
			{
				if (($jsonHash{$j}{"run name"} eq $run) and ($jsonHash{$j}{"lane"} eq $lane))
				{
					$targetSize = $jsonHash{$j}{"target size"};
					$targetSize =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
					$numberOfTargets = $jsonHash{$j}{"number of targets"};
					$numberOfTargets =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;

					if (exists $jsonHash{$j}{barcode})
					{
						$linkName = "${run}_${lane}_$jsonHash{$j}{barcode}";
					}
					else
					{
						$linkName = "${run}_${lane}";
					}

					print "<tr>\n";
					print "<td><a href=\"#$linkName\">$jsonHash{$j}{lane}</a></td>";
					if (exists $jsonHash{$j}{barcode})
					{
						print "<td><a href=\"#$linkName\">$jsonHash{$j}{barcode}</a></td>";
					}
					else
					{
						print "<td>none</td>";
					}
					print "<td><a href=\"#$linkName\">$jsonHash{$j}{library}</a></td>";
					print "<td>$targetSize</td>";
					print "<td>$numberOfTargets</td>";


					$targetSize = $jsonHash{$j}{"target size"};		# so it isn't formatted with commas
					for my $collapsed ("non collapsed", "collapsed")
					{
						for my $i (@coverageXs)
						{
							if (exists $jsonHash{$j}{"$collapsed bases covered"}{$i})
							{
								$basesCovered = $jsonHash{$j}{"$collapsed bases covered"}{$i};
							}
							else
							{
								$basesCovered = 0;
							}

							$basesCovered = $basesCovered / $targetSize;
							$basesCovered = sprintf "%.2f%%", $basesCovered * 100;

							print "<td>$basesCovered</td>";


						}
					}
				}
			}
		}
	
		print "</tr>\n";
		print "</table>\n";
		print "<p></p>\n";

	}


	# print image thumbnail table

	print "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";
	
	print "<th>Lane</th>";
	print "<th>Barcode</th>";
	print "<th>Library</th>";
	print "<th>Read Breakdown</th>";
	if ($ends eq "paired end")
	{
		print "<th>Insert Distribution</th>";
	}
	if ($showReadLengthHist == 1)
	{
		print "<th>Read length Histogram</th>";
	}
	print "<th>Quality Histogram</th>";
	print "<th>Quality by Cycle</th>";
	print "<th>Mismatch by Cycle</th>";
	print "<th>Indels by Cycle</th>";
	print "<th>Soft Clip by Cycle</th>";
	if ($showHardClip == 1)
	{
		print "<th>Hard Clip by Cycle</th>";
	}

	print "\n</tr>\n</thead>\n<tbody>\n";

	for my $lane (sort keys %{ $runList{$run} })
	{
		for my $j (keys %jsonHash)
		{
			if (($jsonHash{$j}{"run name"} eq $run) and ($jsonHash{$j}{"lane"} eq $lane))
			{
				print "<tr>\n";
				print "<td><a href=\"#$linkName\">$jsonHash{$j}{lane}</a></td>";
				if (exists $jsonHash{$j}{barcode})
				{
					print "<td><a href=\"#$linkName\">$jsonHash{$j}{barcode}</a></td>";
				}
				else
				{
					print "<td>none</td>";
				}
				print "<td><a href=\"#$linkName\">$jsonHash{$j}{library}</a></td>";

				print "<td><a href=\"${j}.graphs/readPie.png\"><img src=\"${j}.graphs/readPie.png\" width=\"100\" height=\"100\"/></a></td>";
				if ($ends eq "paired end")
				{
					print "<td><a href=\"${j}.graphs/insert.png\"><img src=\"${j}.graphs/insert.png\" width=\"100\" height=\"100\"/></a></td>";
				}
				if ($showReadLengthHist == 1)
				{
					print "<td><a href=\"${j}.graphs/readLength.png\"><img src=\"${j}.graphs/readLength.png\" width=\"100\" height=\"100\"/></a></td>";
				}
				print "<td><a href=\"${j}.graphs/qualHist.png\"><img src=\"${j}.graphs/qualHist.png\" width=\"100\" height=\"100\"/></a></td>";
				print "<td><a href=\"${j}.graphs/qualCycle.png\"><img src=\"${j}.graphs/qualCycle.png\" width=\"100\" height=\"100\"/></a></td>";
				print "<td><a href=\"${j}.graphs/misCycle.png\"><img src=\"${j}.graphs/misCycle.png\" width=\"100\" height=\"100\"/></a></td>";
				print "<td><a href=\"${j}.graphs/indelCycle.png\"><img src=\"${j}.graphs/indelCycle.png\" width=\"100\" height=\"100\"/></a></td>";
				print "<td><a href=\"${j}.graphs/softCycle.png\"><img src=\"${j}.graphs/softCycle.png\" width=\"100\" height=\"100\"/></a></td>";
				if ($showHardClip == 1)
				{
					print "<td><a href=\"${j}.graphs/hardCycle.png\"><img src=\"${j}.graphs/hardCycle.png\" width=\"100\" height=\"100\"/></a></td>";
				}
				print "\n</tr>\n";
			}
		}
	}
	print "</tbody>\n</table>\n";
}






# print lane specific info


for my $run (sort keys %runList)
{
	for my $lane (sort keys %{ $runList{$run} })
	{
		for my $j (keys %jsonHash)
		{
			if (($jsonHash{$j}{"run name"} eq $run) and ($jsonHash{$j}{"lane"} eq $lane))
			{
				if (exists $jsonHash{$j}{"barcode"})
				{
					print "<h2><a name=\"${run}_${lane}_$jsonHash{$j}{\"barcode\"}\">$run Lane: $lane Barcode: $jsonHash{$j}{\"barcode\"}</a></h2>\n";
				}
				else
				{
					print "<h2><a name=\"${run}_${lane}\">$run Lane: $lane</a></h2>\n";
				}
				print "<p><a href=\"#$run\">Back to $run.</a></p>\n";
				print "<ul>\n";
				print "<li>Library: $jsonHash{$j}{\"library\"}</li>";
				print "<li>Target file: $jsonHash{$j}{\"target file\"}</li>";
				$targetSize = $jsonHash{$j}{"target size"};
				$targetSize =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
				$numberOfTargets = $jsonHash{$j}{"number of targets"};
				$numberOfTargets =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
				print "<li>Target size: $targetSize bp</li>";
				print "<li>Number of targets: $numberOfTargets</li>";
				if($jsonHash{$j}{"workflow name"}) { print "<li>Workflow Name: $jsonHash{$j}{'workflow name'}</li>"; }
				if($jsonHash{$j}{"workflow version"}) { print "<li>Workflow Version: $jsonHash{$j}{'workflow version'}</li>"; }
				if($jsonHash{$j}{"workflow run creation timestamp"}) { print "<li>Workflow Run Creation Timestamp: $jsonHash{$j}{'workflow run creation timestamp'}</li>"; }
                if($jsonHash{$j}{"bam file creation timestamp"}) { print "<li>Bam File Creation Timestamp: $jsonHash{$j}{'bam file creation timestamp'}</li>"; }
				print "</ul>\n";

				if ($printAllImages == 1)
				{
					print "<img src=\"${j}.graphs/readPie.png\"/>\n";
					print "<img src=\"${j}.graphs/insert.png\"/>\n";
					print "<img src=\"${j}.graphs/qualHist.png\"/>\n";
					print "<img src=\"${j}.graphs/qualCycle.png\"/>\n";
					print "<img src=\"${j}.graphs/misCycle.png\"/>\n";
					print "<img src=\"${j}.graphs/indelCycle.png\"/>\n";
					print "<img src=\"${j}.graphs/softCycle.png\"/>\n";
					print "<img src=\"${j}.graphs/hardCycle.png\"/>\n";
				}
			}
		}
	}
}

print "</body>\n</html>\n";


# draw graphs

my $insertMax = 0;		# find the largest insert mean * 6 insert stdev in the data
for my $j (keys %jsonHash)
{
	if ($insertMax <= ($jsonHash{$j}{"insert mean"} + ($jsonHash{$j}{"insert stdev"} * 6)))
	{
		$insertMax = $jsonHash{$j}{"insert mean"} + ($jsonHash{$j}{"insert stdev"} * 6);
	}
}

my $title;


for my $j (keys %jsonHash)
{
	if (exists $jsonHash{$j}{"barcode"})
	{
		$title = $jsonHash{$j}{"run name"} . " Lane: " . $jsonHash{$j}{"lane"} . " Barcode: " . $jsonHash{$j}{"barcode"} . "\\n" . $jsonHash{$j}{"library"};
	}
	else
	{
		$title = $jsonHash{$j}{"run name"} . " Lane: " . $jsonHash{$j}{"lane"} . "\\n" . $jsonHash{$j}{"library"};
	}

	warn "graphing $j\n";

	`$scriptPath/jsonToGraphs.pl $j`;
}
