#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Std;
use vars qw/ %opt /;

use JSON;

my $scriptLocation = "~rdenroche/svn/seq_prod_bio/pipeline/jsonReport/";

my $wgsTargetBed = "/oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed";
my $exomeTargetBed = $wgsTargetBed;

my $l;
my @f;

my %summaryHash;		# $summaryHash{sample}{seqtype}{sampletype/gatk} = 1;

sub usage
{
    print "\nUsage is runSamStatsProject.pl [options].\n";
	print "NOTE: should be run from the project/analysis dir (i.e. /oicr/data/archive/projects/PCSI/analysis/)\n";
    print "Options are as follows:\n";
	print "\t-t target.bed to calculate *exome* reads on target and coverage against.  Default is $wgsTargetBed .  wgs will be calculated against the default.\n";
	print "\t-v paths have versions in them (i.e. .../novoalign/V2.07.09/hg19_random/...).  Default is paths do not have versions.\n";
    print "\t-a regenerate all json files.  Default is to only generate new json files if a json file doesn't exist, or if the json's timestamp is older than the bam's timestamp.\n";
	print "\t-h displays this usage.\n";

    die "\n@_\n\n";
}

my $opt_string = "t:avh";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
    usage("Help requested.");
}

my $generateAll = 0;
if (exists $opt{a})
{
	$generateAll = 1;
}

if (exists $opt{t})
{
	$exomeTargetBed = $opt{t};
}

my $laneBamFiles;
my @laneBamFiles;

my $mergedBamFiles;
my @mergedBamFiles;

my $superMergeGATKs;
my @superMergeGATKs;

if (exists $opt{v})
{
	# assuming the old version tree structure!  (i.e. .../novoalign/V2.07.09/hg19_random/...)
	
	$laneBamFiles = `ls */*/*/novoalign/*/*/libraries/*/bam/*.bam`;
	chomp $laneBamFiles;
	@laneBamFiles = split(/\n/, $laneBamFiles);
	
	$mergedBamFiles = `ls */*/*/novoalign/*/*/collapsed/merge/*.bam`;
	chomp $mergedBamFiles;
	@mergedBamFiles = split(/\n/, $mergedBamFiles);

	$superMergeGATKs = `ls -d */*superMerge/*/novoalign/*/*/gatk/gatk*/`;
	chomp $superMergeGATKs;
	@superMergeGATKs = split(/\n/, $superMergeGATKs);
}
else
{
	$laneBamFiles = `ls */*/*/novoalign/*/libraries/*/bam/*.bam`;
	chomp $laneBamFiles;
	@laneBamFiles = split(/\n/, $laneBamFiles);

	$mergedBamFiles = `ls */*/*/novoalign/*/collapsed/merge/*.bam`;
	chomp $mergedBamFiles;
	@mergedBamFiles = split(/\n/, $mergedBamFiles);

	$superMergeGATKs = `ls -d */*superMerge/*/novoalign/*/gatk/gatk*/`;
	chomp $superMergeGATKs;
	@superMergeGATKs = split(/\n/, $superMergeGATKs);
}


unless (-d "./jsonReport")
{
	mkdir "./jsonReport";
}

my $project;

my $cwd = `pwd`;
chomp $cwd;
if ($cwd =~ /.*projects\/(.*?)\//)
{
	$project = $1;
}
elsif ($cwd =~ /.*cpcgene.*/)
{
	$project = "CPC-GENE";
}
elsif ($cwd =~ /.*cpc-gene.*/)
{
	$project = "CPC-GENE";
}
elsif ($cwd =~ /\/\.mounts\/sata\/bas021\/seqprodbio\/archive\/projects\/CPC-GENE_bas021/)
{
	$project = "CPC-GENE";
}
elsif ($cwd =~ /.*HepatoCarFr.*/)
{
	$project = "HepatoCarcinomaFrance";
}
elsif ($cwd =~ /.*DynamicsOfChemotherapyResponse.*/)
{
	$project = "DynamicsOfChemotherapyResponse";
}
elsif ($cwd =~ /.*CPCG.*/)
{
	$project = "CPC-GENE";
}
else
{
	die "Can't parse $cwd\n";
}

my $needToJson;

my $sampleGroup;
my $sample;
my $library;

my $basename;

my $runName;
my $lane;
my $ends;
my $instrument;

my $seqType;

my $path;

my %jsonHash;
my $jsonText;

my $jsonOutFile;

open (GOJSON, ">./jsonReport/goJson") or die "Couldn't open ./jsonReport/goJson\n";
open (NOBAM, ">./jsonReport/brokenBamLinks.txt") or die "Couldn't open ./jsonReport/brokenBamLinks.txt\n";

my $holdList = "";

for my $bam (@laneBamFiles)
{
	%jsonHash = ();
	unless (-e $bam)
	{
        warn "$bam is a broken link, writing to ./jsonReport/brokenBamLinks.txt\n";
		print NOBAM "$bam\n";
	}

	if ($bam =~ /(.*?)\/(.*?)\/(.*?)\/.*libraries\/(.*?)\/bam\/(.*)/)
	{
		$sampleGroup = $1;
		$sample = $2;
		$seqType = $3;
		$library = $4;
		$basename = $5;
		print $basename . "\n";
	}
	else
	{
        die "couldn't parse $bam\n";
	}
	if ($basename =~ /(.*)_s_(.)_sequence_(.*)_novoalign.*\.bam/)
	{
		$runName = $1;
		$lane = $2;
		$ends = $3;
	}
	# SWID_61719_AOE_0001_Ly_R_PE_414_WG_120126_h801_0067_AD0PN5ACXX_NoIndex_L001_R1_001.fastq.gz.annotated.bam
	elsif ($basename =~ /SWID_.*?_.*?_.*?_.*?_.*?_.*?_.*?_(.*)_.*_L00(.)_/)
	{
		$runName = $1;
		$lane = $2;
		$ends = "don't use this";
	}
	elsif ($basename =~ /PMRA_.*?_.*?_.*?_.*?_.*?_.*?_(.*)_.*?_L00(.)\.bam/)		# if this is a common format you might expand it to not be specifically for PMRA
	{
		$runName = $1;
		$lane = $2;
		$ends = "don't use this";
	}
	elsif ($basename =~ /OVCA.*?_sequence_(.*?)_novoalign\.sam\.sorted\.bam/)		# Fouad has weird names
	{
		$runName = "J_Brenton_sample";
		$lane = 1;
		$ends = $1;
	}
	elsif ($basename =~ /(.*)_(.*)_(.*)_CHC_.*\.bam/)		# bams from the HepatoCarFr project look like this: 61UW3AAXX_5_B009YW2_CHC_0018_nn_P.bam
	{
		$runName = $1;
		$lane = $2;
#		$barcode = $3;		# don't need barcode
	}
	elsif ($basename =~ /^(.*?_.*?_.*?_.*?)_(.)_/)
	{
		$runName = $1;
		$lane = $2;
	}
	else
	{
		#die "couldn't parse $basename\n";
		warn "couldn't parse $basename\n";
	}
	if ($runName =~ /.*?_(.*?)_/)
	{
		$instrument = $1;
	}
	elsif ($project eq "HepatoCarcinomaFrance")
	{
		$instrument = "Fr";
	}
	else
	{
		die "couldn't parse $runName\n";
	}

	if ($bam =~ /(.*)\/$basename/)
	{
		$path = $1;
	}



	$jsonHash{"sample group"} = $sampleGroup;
	$jsonHash{"sample"} = $sample;
	$jsonHash{"run name"} = $runName;
	$jsonHash{"instrument"} = $instrument;
	$jsonHash{"lane"} = $lane;
	$jsonHash{"library"} = $library;
	$jsonHash{"sequencing type"} = $seqType;
	$jsonHash{"last modified"} = (stat($bam))[9];

	$jsonText = encode_json(\%jsonHash);
#	$jsonText =~ s/\"/\\\"/g;		don't need to escape since we're putting the command in a file now

	$jsonOutFile = "${sample}_${seqType}_${runName}_${lane}_${library}_lane.json";

	$needToJson = 0;

	unless (-e "$path/$jsonOutFile")
	{
		$needToJson = 1;
	}
	elsif ((stat($bam))[9] > (stat("$path/$jsonOutFile"))[9])		# if the bam file is newer than the json file (does this look at the bam file itself, or just the link?)
	{
		$needToJson = 1;
	}
	elsif ((-s "$path/$jsonOutFile") == 0)
	{
		$needToJson = 1;
	}

	$path = "../$path";		# since we'll be down in the jsonReport folder
	if (($needToJson == 1) or ($generateAll == 1))
	{
		open (COMMAND, ">./jsonReport/$jsonOutFile.command") or die "Couldn't open >./jsonReport/$jsonOutFile.command\n";

		if ($seqType eq "exome")
		{
			print COMMAND "module load samtools; samtools view ../$bam | $scriptLocation/samStats.pl -r $exomeTargetBed -j '$jsonText' > $path/$jsonOutFile\n";
		}
		else
		{
			print COMMAND "module load samtools; samtools view ../$bam | $scriptLocation/samStats.pl -r $wgsTargetBed -j '$jsonText' > $path/$jsonOutFile\n";
		}
		close COMMAND;

		print GOJSON "qsub -cwd -b y -V -N stat${lane}_${runName} -e ${jsonOutFile}.log -o ${jsonOutFile}.log \"bash ./$jsonOutFile.command\"\n";
		$holdList = "stat${lane}_${runName},$holdList";
	}

	`ln -s $path/$jsonOutFile ./jsonReport/$jsonOutFile`;
}


my $exomeHistFile;
my $exomeHistImage;
my %exomeHist;

for my $bam (@mergedBamFiles)
{
	%jsonHash = ();

    unless (-e $bam)
    {
        warn "$bam is a broken link, writing to ./jsonReport/brokenBamLinks.txt\n";
		print NOBAM "$bam\n";
    }

    if ($bam =~ /^(.*?)\/(.*?)\/(.*?)\/.*collapsed\/merge\/(.*)/)
    {
        $sampleGroup = $1;
        $sample = $2;
		$seqType = $3;
        $basename = $4;
    }
    else
    {
        die "couldn't parse $bam\n";
    }

    if ($bam =~ /(.*)\/$basename/)
    {
        $path = $1;
    }

	$summaryHash{$seqType}{$sampleGroup}{$sample} = 1;

	if ($seqType eq "exome")
	{
		if ($bam =~ /^(.*)\/merge/)
		{
			$exomeHistFile = "$1/coverage/$sample.cumhist.dat";
			$exomeHistImage = "$1/coverage/$sample.coverage.png";
		}
		else
		{
			die "couldn't parse $bam for exome coverage stats\n";
		}

		if (-e $exomeHistFile)
		{
			if ((stat($bam))[9] > (stat($exomeHistFile))[9])
			{
				warn "exome coverage stats are out of date for $bam\n";
			}

			open (CUMHIST, "$exomeHistFile") or die "Couldn't open $exomeHistFile.\n";
			while ($l = <CUMHIST>)
			{
				chomp $l;
				@f = split (/\t/, $l);
				$exomeHist{$f[0]} = $f[1];
			}
			close CUMHIST;
		}
		else
		{
			%exomeHist = ();
			$exomeHistImage = "no image";
			warn "no bed coverage stats for $bam\n";
		}
	}

    $jsonHash{"sample group"} = $sampleGroup;
    $jsonHash{"sample"} = $sample;
	$jsonHash{"sequencing type"} = $seqType;
	$jsonHash{"library"} = "merged";
	$jsonHash{"last modified"} = (stat($bam))[9];
	$jsonHash{"collapsed bases covered"} = \%exomeHist;
	$jsonHash{"exome histogram image"} = $exomeHistImage;

    $jsonText = encode_json(\%jsonHash);
#    $jsonText =~ s/\"/\\\"/g;		writing to a command file now, so no need to escape

    $jsonOutFile = "${sample}_${seqType}_merged.json";

	$needToJson = 0;

	unless (-e "$path/$jsonOutFile")
	{
		$needToJson = 1;
	}
	elsif ((stat($bam))[9] > (stat("$path/$jsonOutFile"))[9])		# if the bam file is newer than the json file (does this look at the bam file itself, or just the link?)
	{
		$needToJson = 1;
	}
	elsif ((-s "$path/$jsonOutFile") == 0)
	{
		$needToJson = 1;
	}


	$path = "../$path";		# since we'll be down in the jsonReport folder
	if (($needToJson == 1) or ($generateAll == 1))
	{
		open (COMMAND, ">./jsonReport/$jsonOutFile.command") or die "Couldn't open >./jsonReport/$jsonOutFile.command\n";
		if ($seqType eq "exome")
		{
			print COMMAND "module load samtools; samtools view ../$bam | $scriptLocation/samStats.pl -r $exomeTargetBed -j '$jsonText' > $path/$jsonOutFile\n";
		}
		else
		{
			print COMMAND "module load samtools; samtools view ../$bam | $scriptLocation/samStats.pl -r $wgsTargetBed -j '$jsonText' > $path/$jsonOutFile\n";
		}
		close COMMAND;

		print GOJSON "qsub -cwd -b y -V -N stat${sample} -e ${jsonOutFile}.log -o ${jsonOutFile}.log \"bash ./$jsonOutFile.command\"\n";
		$holdList = "stat${sample},$holdList";
	}

    `ln -s $path/$jsonOutFile ./jsonReport/$jsonOutFile`;
}

my $gatkFileName;

for my $gatkDir (@superMergeGATKs)
{
	
	if ($gatkDir =~ /^(.*?)\/(.*?)\/(.*?)\/.*gatk\//)
	{
		$sampleGroup = $1;
		$sample = $2;
		$seqType = $3;
	}
	else
	{
		die "Couldn't parse $gatkDir.\n";
	}

	$summaryHash{$seqType}{$sampleGroup}{"gatk"} = 1;

	$gatkFileName = "${sampleGroup}_${seqType}_gatk.json";

	print GOJSON "qsub -cwd -b y -V -o $gatkFileName.log -e $gatkFileName.log -N stat_gatk$sampleGroup \"$scriptLocation/gatkStats.pl ../$gatkDir > ../$gatkDir/$gatkFileName\"\n";

	`ln -s ../$gatkDir/$gatkFileName ./jsonReport/$gatkFileName`;
}

my $user = `whoami`;
chomp $user;

unless ($holdList eq "")
{
	print GOJSON "qsub -cwd -b y -V -N jsonReport$project -hold_jid $holdList \"~rdenroche/jsonToGenericProjectReport.pl -P $project *.json\"\n";
	print GOJSON "qsub -cwd -b y -V -N jsonReport$project -hold_jid $holdList \"~rdenroche/jsonToShortProjectReport.pl -G -P $project -n short *.json\"\n";
}
else # no need to hold_jid
{
	print GOJSON "qsub -cwd -b y -V -N jsonReport$project \"~rdenroche/jsonToGenericProjectReport.pl -P $project *.json\"\n";
	print GOJSON "qsub -cwd -b y -V -N jsonReport$project \"~rdenroche/jsonToShortProjectReport.pl -G -P $project -n short *.json\"\n";
}
print GOJSON "qsub -cwd -b y -e mail.log -o mail.log -N jsonEmail$project -hold_jid jsonReport$project \"mail -s ${project}_report_is_done ${user}\@oicr.on.ca < mail.txt\"\n";

close GOJSON;

my $htmlPath;
if ($cwd =~ /.*archive\/(.*)/)
{
	$htmlPath = $1;
}
elsif ($cwd =~ /.*analysis\/(.*)/)
{
	$htmlPath = "projects/${project}/analysis\/$1";
}
else
{
	$htmlPath = "projects/${project}/analysis";
}

open (MAIL, ">./jsonReport/mail.txt") or die "Couldn't open ./jsonReport/mail.txt\n";
print MAIL "Report is here: http://hn1.hpc.oicr.on.ca/archive/${htmlPath}/jsonReport/${project}_report.html\n";
close MAIL;


for my $seq (sort keys %summaryHash)
{
	print "\n$seq\n";
	for my $samp (sort keys %{ $summaryHash{$seq} })
	{
		print "$samp,";
		for my $type (qw(R P X C))
		{
			if (exists $summaryHash{$seq}{$samp}{"${samp}${type}"})
			{
				print "$type";
			}
			print ",";
		}

		if (exists $summaryHash{$seq}{$samp}{"gatk"})
		{
			print "gatk";
		}
		print ",\n";

	}
}
