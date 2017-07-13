#!/usr/bin/perl

use strict;
use warnings;

use JSON;


use Getopt::Std;
use Cwd 'abs_path';
use File::Basename;
use vars qw/ %opt /;

my $reportPrefix = "";
my $projectName = "";

my $scriptPath = dirname(abs_path($0));

sub usage
{
        print "\nUsage is jsonToGenericProjectReport.pl [options] path/to/*.json\n";
        print "Options are as follows:\n";
		print "\t-P project name.  Default is an attempt to parse from path\n";
		print "\t-n report name prefix.  Default is $reportPrefix.\n";
		print "\t-c show bases covered data.  Default is no bases covered.\n";
        print "\t-r show read length histogram.  Default is no histogram.\n";
        print "\t-p print all images.  Default is to only show thumbnails (with links).\n";
		print "\t-H show hard clip stats and graph.  Default is off.\n";
		print "\t-G don't graph!\n";
        print "\t-h displays this usage message.\n";

        die "\n@_\n\n";
}


my $showReadLengthHist = 0;
my $printAllImages = 0;
my $showHardClip = 0;
my $showBasesCovered = 0;

my $opt_string = "P:n:crpHGh";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
	usage("Help requested.");
}
if (exists $opt{P})
{
	$projectName = $opt{P};
}
if (exists $opt{n})
{
	$reportPrefix = $opt{n};
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

my %jsonOrder;
my %laneCount;

my %jsonGATK;


my $line;

my $gatkGroup;

my $seqType;

for my $file (@jsonFiles)
{
	open (FILE, $file) or die "Couldn't open $file.\n";

	if ($line = <FILE>)
	{
		$jsonHash{$file} = decode_json($line);


		#TEMPORARY: get merged library status from file name
		if ($file =~ /_merged\./)
		{

			$jsonOrder{$jsonHash{$file}{"sample group"}}{$jsonHash{$file}{"sequencing type"}}{$jsonHash{$file}{"sample"}}{"merged"}{$file} = 1;
		}
		elsif ($file =~ /(.*)_(.*?)_gatk\.json/)
		{
			$jsonGATK{$1}{$2}{"gatk"} = $file;
		}
		else
		{
			$jsonOrder{$jsonHash{$file}{"sample group"}}{$jsonHash{$file}{"sequencing type"}}{$jsonHash{$file}{"group id"}}{$jsonHash{$file}{"sample"}}{$jsonHash{$file}{"library"}}{$file} = 1;
			$laneCount{$jsonHash{$file}{"sample group"}}{$jsonHash{$file}{"sequencing type"}}{$jsonHash{$file}{"group id"}}{$jsonHash{$file}{"sample"}}++;
		}


		if (exists $opt{P})
		{
			# already have the project name
		}
		elsif ($file =~ /(.*?)0.*/)
		{
			$projectName = $1;
		}
	}
	else
	{
		warn "No data found in $file!\n";
	}
	close FILE;
}

my $date = `date`;
chomp $date;


unless (-e "sorttable.js")
{
	`ln -s ~rdenroche/sorttable.js`;
}


my $averageReadLength;

my $percentOnTarget;
my $numReadsOnTarget;
my $numBasesOnTarget;

my $readsPerSP;

my $estimatedYield;
my $estimatedCoverage;

my $rawReads;
my $rawYield;
my $ends = "paired end";		# should determine this value for each sample page
my $insertMean;
my $insertStdev;
my $readLength;
my $mapRate;
my $errorRate;
my $softClipRate;
my $hardClipRate;
my $readsPerStartPoint;
my $onTargetRate;

my $totalReads = 0;
my $totalRawYield = 0;
my $totalEstYield = 0;

my $qualCut;
my %qualCuts;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
my @month = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

my $bestDepth = 0;

my $subName;
if ($reportPrefix ne "")
{
	open (HTML, ">${reportPrefix}_${projectName}_report.html") or die "Couldn't open ${reportPrefix}_${projectName}_report.html\n";
	warn "Printing to ${reportPrefix}_${projectName}_report.html\n";
	$subName = "${reportPrefix}_${projectName}";
}
else
{
	open (HTML, ">${projectName}_report.html") or die "Couldn't open ${projectName}_report.html\n";
	warn "Printing to ${projectName}_report.html\n";
	$subName = "$projectName";
}

print HTML "<html>\n<head>\n";
print HTML "<title>$projectName Report</title>\n";
print HTML "<script src=\"./sorttable.js\"></script>\n";
print HTML "<style type=\"text/css\">\nth, td {\n  padding: 3px !important;\n}\ntable\n{\nborder-collapse:collapse;\n}\n/* Sortable tables */\ntable.sortable thead {\n\tbackground-color:#eee;\n\tcolor:#000000;\n\tfont-weight: bold;\n\tcursor: default;\n}\n</style>\n";
print HTML "</head>\n<body>\n";
print HTML "<p>Generic project report generated on $date.</p>\n";

for my $sampleGroup (sort keys %jsonOrder)
{
	print HTML "<h1><a name=\"$sampleGroup\">$sampleGroup</a></h1>\n";
	for my $seqType (sort keys %{ $jsonOrder{$sampleGroup} })
	{
		print HTML "<h2>$seqType</h2>\n";
		print HTML "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";
		print HTML "<th>Sample</th>";
		print HTML "<th># Lanes</th>";
		print HTML "<th>Last Modified</th>";
        print HTML "<th>Group ID</th>";
        print HTML "<th>External Name</th>";
       	print HTML "<th>Avg Read Length</th>";
		print HTML "<th>Error %</th>";
		print HTML "<th>Soft Clip %</th>";
		print HTML "<th>% Reads on Target</th>";
		print HTML "<th># Reads on Target</th>";
		print HTML "<th># Bases on Target</th>";
		print HTML "<th>Reads/SP</th>";
		print HTML "<th>Coverage</th>";
		if (($seqType eq "exome") or ($seqType eq "EX"))
		{
			print HTML "<th>% at 8x</th>";
			print HTML "<th>90% covered at</th>";
			print HTML "<th>Graph</th>";
		}

		print HTML "\n</thead>\n<tbody>\n";


		for my $groupId (sort keys %{ $jsonOrder{$sampleGroup}{$seqType} })
		{
        for my $sample (sort keys %{ $jsonOrder{$sampleGroup}{$seqType}{$groupId} })
        {
			print HTML "<tr>\n";
			print HTML "<td><a href=\"${subName}_${sample}_${seqType}_${groupId}_report.html\">$sample</a></td>";
			
			print HTML "<td><a href=\"${subName}_${sample}_${seqType}_${groupId}_report.html\">$laneCount{$sampleGroup}{$seqType}{$groupId}{$sample}</a></td>";




			if (exists $jsonOrder{$sampleGroup}{$seqType}{$groupId}{$sample}{"merged"})
			{
				for my $j (keys %{ $jsonOrder{$sampleGroup}{$seqType}{$groupId}{$sample}{"merged"}})
				{
					($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($jsonHash{$j}{"last modified"});
					$year += 1900;
					print HTML "<td>$month[$mon] $mday, $year</td>";

					if ($jsonHash{$j}{"number of ends"} eq "single end")
					{
						$averageReadLength = $jsonHash{$j}{"read ? average length"};
						$numBasesOnTarget = $jsonHash{$j}{"reads on target"} * $jsonHash{$j}{"read ? average length"};

						if ($jsonHash{$j}{"read ? average length"} =~ m/^\d+\.\d+$/)
						{
							$averageReadLength = sprintf "%.2f", $jsonHash{$j}{"read ? average length"};
						}
					}
					else
					{
						$averageReadLength = $jsonHash{$j}{"read 1 average length"} . ", " . $jsonHash{$j}{"read 2 average length"};
						$numBasesOnTarget = ($jsonHash{$j}{"reads on target"} * ($jsonHash{$j}{"read 1 average length"} + $jsonHash{$j}{"read 2 average length"})) / 2;

						if (($jsonHash{$j}{"read 1 average length"} =~ m/^\d+\.\d+$/) or ($jsonHash{$j}{"read 2 average length"} =~ m/^\d+\.\d+$/))
						{
							$averageReadLength = sprintf "%.2f", $jsonHash{$j}{"read 1 average length"};
							$averageReadLength = "$averageReadLength," . sprintf "%.2f", $jsonHash{$j}{"read 2 average length"};
						}
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

					$numReadsOnTarget = $jsonHash{$j}{"reads on target"};

					$percentOnTarget = 0;
					if ($jsonHash{$j}{"mapped reads"} > 0)
					{
						$percentOnTarget = ($jsonHash{$j}{"reads on target"} / $jsonHash{$j}{"mapped reads"}) * 100;
					}

					$readsPerSP = sprintf "%.2f", $jsonHash{$j}{"reads per start point"};
					
					$estimatedYield = int($jsonHash{$j}{"aligned bases"} * ($percentOnTarget / 100));
					

					$estimatedCoverage = $estimatedYield / $jsonHash{$j}{"target size"};
					$estimatedCoverage = sprintf "%.2f", $estimatedCoverage;

					$percentOnTarget = sprintf "%.2f%%", $percentOnTarget;
					$numReadsOnTarget =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
					$numBasesOnTarget = int($numBasesOnTarget);
					$numBasesOnTarget =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;

                    if (exists $jsonHash{$j}{"group id"})
                    {
                        my $group = "<td";
                        if (exists $jsonHash{$j}{"group id description"}) {
                              $group .= " title='" . $jsonHash{$j}{"group id description"} . "'";
                        }
                        print HTML $group . ">$jsonHash{$j}{'group id'}</td>";
                    }
                    else
                    {
                        print HTML "<td class='na'>na</td>";
                    }
                    if (exists $jsonHash{$j}{"external name"})
                    {
                        print HTML "<td>$jsonHash{$j}{'external name'}</td>";
                    }
                    else
                    {
                        print HTML "<td class='na'>na</td>";
                    }
					print HTML "<td>$averageReadLength</td>";
					print HTML "<td>$errorRate</td>";
					print HTML "<td>$softClipRate</td>";
					print HTML "<td>$percentOnTarget</td>";
					print HTML "<td>$numReadsOnTarget</td>";
					print HTML "<td>$numBasesOnTarget</td>";
					print HTML "<td>$readsPerSP</td>";
					print HTML "<td>${estimatedCoverage}x</td>";

					if (($seqType eq "exome") or ($seqType eq "EX"))
					{
						if (exists $jsonHash{$j}{"collapsed bases covered"}{8})
						{
							print HTML "<td>$jsonHash{$j}{\"collapsed bases covered\"}{8}%</td>";
							$bestDepth = 0;
							for my $i (sort { $a <=> $b } keys %{ $jsonHash{$j}{"collapsed bases covered"} })
							{
								if ($jsonHash{$j}{"collapsed bases covered"}{$i} >= 90)
								{
									$bestDepth = $i;
								}
							}
							if ($bestDepth == 50)
							{
								$bestDepth = ">50";
							}
						
							print HTML "<td>${bestDepth}x</td>";
							print HTML "<td><a href=\"../$jsonHash{$j}{\"exome histogram image\"}\"><img src=\"../$jsonHash{$j}{\"exome histogram image\"}\" width=\"27\" height=\"20\"/></a></td>";
						}
						else
						{
							print HTML "<td colspan=\"3\">Bed coverage not found</td>";
						}
						
					}
				}
			}
			else
			{
				if (($seqType eq "exome") or ($seqType eq "EX"))
				{
					print HTML "<td colspan=\"12\">Merge has not been generated</td>";
				}
				else
				{
					print HTML "<td colspan=\"11\">Merge has not been generated</td>";
				}
			}
			print HTML "</tr>";
		}
		}
		print HTML "</tbody>\n</table>\n";

		if (exists $jsonGATK{$sampleGroup}{$seqType}{"gatk"})
		{
			print HTML "<p>";
			print HTML "<table border=\"1\" class=\"sortable\"><thead><tr>\n";
			print HTML "<th>Super Merge Name</th>";
			print HTML "<th>GATK Complete</th>";
			print HTML "<th>Pass Filter SNP Candidates</th>";
			print HTML "<th>Indel Candidates</th>";
			print HTML "\n</tr></thead><tbody><tr>\n";

			if (exists $jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"sample"})
			{
				print HTML "<td>$jsonHash{$jsonGATK{$sampleGroup}{$seqType}{\"gatk\"}}{\"sample\"}</td>";
			}
			else
			{
				print HTML "<td>Incomplete</td>";
			}


			if (defined $jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"snp date"})
			{
				if ($jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"snp date"} > 0)
				{
					($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"snp date"});
					$year += 1900;
					print HTML "<td>$month[$mon] $mday, $year</td>";
				}
				else
				{
					print HTML "<td>n/a</td>";
				}
			}
			else
			{
				print HTML "<td>n/a</td>";
			}

			$jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"pass snp count"} =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
			print HTML "<td>" . $jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"pass snp count"} . "</td>";

			$jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"indel count"} =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
			print HTML "<td>" . $jsonHash{$jsonGATK{$sampleGroup}{$seqType}{"gatk"}}{"indel count"} . "</td>";

			print HTML "\n</tr></tbody></table>\n";
			print HTML "<p>";

		}
		else
		{
			print HTML "<p>Super merge GATK calls have not been generated.</p>\n";
		}
	}
}


print HTML "</body>\n</html>\n";
close HTML;

# generate sub pages for each sample



for my $sampleGroup (sort keys %jsonOrder)
{
	for my $seqType (sort keys %{ $jsonOrder{$sampleGroup} })
	{
	 
		for my $groupId (sort keys %{ $jsonOrder{$sampleGroup}{$seqType} })
		{
        for my $sample (sort keys %{ $jsonOrder{$sampleGroup}{$seqType}{$groupId} })
        {

			$totalReads = 0;
			$totalRawYield = 0;
			$totalEstYield = 0;

			open (HTML, ">${subName}_${sample}_${seqType}_${groupId}_report.html") or die "Couldn't open ${subName}_${sample}_${groupId}_report.html.\n";

			print HTML "<html>\n<head>\n";
			print HTML "<title>$sample $seqType Report</title>\n";
			print HTML "<script src=\"./sorttable.js\"></script>\n";
			print HTML "<style type=\"text/css\">\n.na { color: #ccc; }\nth, td {\n  padding: 3px !important;\n}\ntable\n{\nborder-collapse:collapse;\n}\n/* Sortable tables */\ntable.sortable thead {\n\tbackground-color:#eee;\n\tcolor:#000000;\n\tfont-weight: bold;\n\tcursor: default;\n}\n</style>\n";
			print HTML "</head>\n<body>\n";
			print HTML "<p>Generic sample report generated on $date.</p>\n";
			print HTML "<p>Return to <a href=\"${subName}_report.html\">$projectName Report</a>.<p>\n";

			print HTML "<h1>$sample $seqType</h1>\n";

			print HTML "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";


			print HTML "<th>Run Name</th>";
			print HTML "<th>Group ID</th>";
            print HTML "<th>External Name</th>";
			print HTML "<th>Lane</th>";
		    print HTML "<th>Barcode</th>";
		    print HTML "<th>Library</th>";

		    if ($ends eq "paired end")
		    {
		        print HTML "<th>Insert Mean (SD)</th>";
		    }

		    print HTML "<th>Read Length</th>";
		    print HTML "<th>Raw Reads</th>";
		    print HTML "<th>Raw Yield</th>";
		    print HTML "<th>Map %</th>";
		    print HTML "<th>Error %</th>";
		    print HTML "<th>Soft Clip %</th>";

		    if ($showHardClip == 1)
			{
				print HTML "<th>Hard Clip %</th>";
			}
			print HTML "<th>Reads/SP</th>";
			print HTML "<th>% on Target</th>";
			print HTML "<th>Estimated Yield*</th>";
			print HTML "<th>Coverage*</th>";


			print HTML "\n</thead>\n<tbody>\n";

			for my $lib (sort keys %{ $jsonOrder{$sampleGroup}{$seqType}{$groupId}{$sample} })
			{
				unless ($lib eq "merged")
				{
					for my $j (keys %{ $jsonOrder{$sampleGroup}{$seqType}{$groupId}{$sample}{$lib}})
					{

		                $rawReads = $jsonHash{$j}{"mapped reads"} + $jsonHash{$j}{"unmapped reads"} + $jsonHash{$j}{"qual fail reads"};
		
		                $rawYield = int($rawReads * $jsonHash{$j}{"average read length"});      # cast to int because average length is only from mapped reads (and may result in ugly decimal)
		
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
	                if ($jsonHash{$j}{"mapped reads"} > 0)
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
	
	                print HTML "<tr>\n";
					print HTML "<td>" . $jsonHash{$j}{"run name"} . "</td>";
					if (exists $jsonHash{$j}{"group id"})
                    {
                        my $group = "<td";
                        if (exists $jsonHash{$j}{"group id description"}) {
                              $group .= " title='" . $jsonHash{$j}{"group id description"} . "'";
                        }
                        print HTML $group . ">$jsonHash{$j}{'group id'}</td>";
                    }
                    else
                    {
                        print HTML "<td class='na'>na</td>";
                    }
                    if (exists $jsonHash{$j}{"external name"})
                    {
                        print HTML "<td>$jsonHash{$j}{'external name'}</td>";
                    }
                    else
                    {
                        print HTML "<td class='na'>na</td>";
                    }
	                print HTML "<td>$jsonHash{$j}{lane}</td>";
	                if (exists $jsonHash{$j}{barcode})
	                {
	                    print HTML "<td>$jsonHash{$j}{barcode}</td>";
	                }
	                else
	                {
	                    print HTML "<td>none</td>";
	                }
	                print HTML "<td>$jsonHash{$j}{library}</td>";
	                if ($ends eq "paired end")
	                {
	                    print HTML "<td>$insertMean ($insertStdev)</td>";
	                }
	                print HTML "<td>$readLength</td>";
	                print HTML "<td>$rawReads</td>";
	                print HTML "<td>$rawYield</td>";
	                print HTML "<td>$mapRate</td>";
	                print HTML "<td>$errorRate</td>";
	                print HTML "<td>$softClipRate</td>";
	                if ($showHardClip == 1)
		                {
		                    print HTML "<td>$hardClipRate</td>";
		                }
		                print HTML "<td>$readsPerStartPoint</td>";
		                print HTML "<td>$onTargetRate</td>";
		                print HTML "<td>$estimatedYield</td>";
		                print HTML "<td>${estimatedCoverage}x</td>";
		                print HTML "\n</tr>\n";
					}
	            }
	        }
	
		    $totalRawYield =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
		    $totalReads =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
		    $totalEstYield =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
		
		    print HTML "</tbody>\n";
		    print HTML "<tfoot>\n<tr>\n";
		    print HTML "<td>Total</td>";
		    if ($ends eq "paired end")
		    {
		        print HTML "<td></td>";
		    }
		    print HTML "<td></td><td></td><td></td><td></td><td></td><td></td><td>$totalReads</td><td>$totalRawYield</td><td></td><td></td>";
		    if ($showHardClip == 1)
		    {
		        print HTML "<td></td>";
		    }
		    print HTML "<td></td><td></td><td></td><td>$totalEstYield</td><td></td>\n";
		    print HTML "</tr>\n";
		    print HTML "</table>\n";
		
		
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
		
		    print HTML "<p>* Estimates exclude unmapped, off target, non-primary or MAPQ < $qualCut reads and use reads per start point to approximate loss due to collapse.</p>\n";

    # print image thumbnail table

    print HTML "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";

	print HTML "<th>Run name</th>";
    print HTML "<th>Lane</th>";
    print HTML "<th>Barcode</th>";
    print HTML "<th>Library</th>";
    print HTML "<th>Read Breakdown</th>";
    if ($ends eq "paired end")
    {
        print HTML "<th>Insert Distribution</th>";
    }
    if ($showReadLengthHist == 1)
    {
        print HTML "<th>Read length Histogram</th>";
    }
    print HTML "<th>Quality Histogram</th>";
    print HTML "<th>Quality by Cycle</th>";
    print HTML "<th>Mismatch by Cycle</th>";
    print HTML "<th>Indels by Cycle</th>";
    print HTML "<th>Soft Clip by Cycle</th>";
    if ($showHardClip == 1)
    {
        print HTML "<th>Hard Clip by Cycle</th>";
    }

    print HTML "\n</tr>\n</thead>\n<tbody>\n";

			for my $lib (sort keys %{ $jsonOrder{$sampleGroup}{$seqType}{$groupId}{$sample} })
			{
				unless ($lib eq "merged")
				{
					for my $j (keys %{ $jsonOrder{$sampleGroup}{$seqType}{$groupId}{$sample}{$lib}})
					{
                print HTML "<tr>\n";
				print HTML "<td>" . $jsonHash{$j}{"run name"} . "</td>";
                print HTML "<td>$jsonHash{$j}{lane}</td>";
                if (exists $jsonHash{$j}{"barcode"})
                {
                    print HTML "<td>$jsonHash{$j}{barcode}</td>";
                }
                else
                {
                    print HTML "<td>none</td>";
                }
                print HTML "<td>$jsonHash{$j}{library}</td>";

                print HTML "<td><a href=\"${j}.graphs/readPie.png\"><img src=\"${j}.graphs/readPie.png\" width=\"100\" height=\"100\"/></a></td>";
                if ($ends eq "paired end")
                {
                    print HTML "<td><a href=\"${j}.graphs/insert.png\"><img src=\"${j}.graphs/insert.png\" width=\"100\" height=\"100\"/></a></td>";
                }
                if ($showReadLengthHist == 1)
                {
                    print HTML "<td><a href=\"${j}.graphs/readLength.png\"><img src=\"${j}.graphs/readLength.png\" width=\"100\" height=\"100\"/></a></td>";
                }
                print HTML "<td><a href=\"${j}.graphs/qualHist.png\"><img src=\"${j}.graphs/qualHist.png\" width=\"100\" height=\"100\"/></a></td>";
                print HTML "<td><a href=\"${j}.graphs/qualCycle.png\"><img src=\"${j}.graphs/qualCycle.png\" width=\"100\" height=\"100\"/></a></td>";
                print HTML "<td><a href=\"${j}.graphs/misCycle.png\"><img src=\"${j}.graphs/misCycle.png\" width=\"100\" height=\"100\"/></a></td>";
                print HTML "<td><a href=\"${j}.graphs/indelCycle.png\"><img src=\"${j}.graphs/indelCycle.png\" width=\"100\" height=\"100\"/></a></td>";
                print HTML "<td><a href=\"${j}.graphs/softCycle.png\"><img src=\"${j}.graphs/softCycle.png\" width=\"100\" height=\"100\"/></a></td>";
                if ($showHardClip == 1)
                {
                    print HTML "<td><a href=\"${j}.graphs/hardCycle.png\"><img src=\"${j}.graphs/hardCycle.png\" width=\"100\" height=\"100\"/></a></td>";
                }
                print HTML "\n</tr>\n";
            }
        }
    }
    print HTML "</tbody>\n</table>\n";


			close HTML;
		}
	}
	}
}

my $title;
unless (exists $opt{G})
{
	for my $j (keys %jsonHash)
	{
		unless ($j =~ /gatk/)
		{
			unless (-d "${j}.graphs")
			{
				mkdir "${j}.graphs";
			}
		
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
	}
}

