#!/usr/bin/perl
#

use strict;
use warnings;

use Getopt::Std;
use Cwd 'abs_path';
use vars qw/ %opt /;

use JSON;

my $scriptPath = dirname(abs_path($0));
my $bamFileDefault = "./bamPaths.txt";

my $bedFileDefault = "/oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed";

sub usage
{
	print "\nUsage is runSamStatsIllumina.pl [options].\n";
	print "Options are as follows:\n";
	print "\t-b is the name of a file containing the name and path to all the bam files, one per line.  Default is $bamFileDefault\n";
	print "\t-t target.bed.  Default is $bedFileDefault.\n";
	print "\t\tNimblegen is /oicr/data/reference/genomes/homo_sapiens/Nimblegen/2.1M_Human_Exome_Annotation_21191/hg19/080904_ccds_exome_rebalfocus_hg19/processed/2.1M_Human_Exome.bed.txt\n";
	print "\t\tAgilent is /oicr/data/genomes/homo_sapiens/Agilent/SureSelect_Whole_Exome_ICGC_Sanger/GRCh37hg19/sanger.exons.bed.hg19\n";
	print "\t\tTruSeq is /oicr/data/reference/genomes/homo_sapiens_mc/TruSeq/TruSeq-Exome-Targeted-Regions-BED-file\n";
	print "\t-T turns on interactive target entry mode.  Default is off.  Overrides -r.\n";
	print "\t-h displays this usage message.\n";

	die "\n@_\n\n";
}

my $runName;
my $flowcell;
my $instrument;
my $bin;

my $bedFile = $bedFileDefault;

my $bamFile = $bamFileDefault;

my $readId;
my $line;
my @fields;
my $doDemultiplex = 0;
my $sampDirFile;

my $lane;

my $interactiveTargetMode = 0;

my $opt_string = "b:t:Th";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
	usage("Help requested.");
}

if (exists $opt{b})
{
	$bamFile = $opt{b};
}

if (exists $opt{t})
{
	$bedFile = $opt{t};
}
if (exists $opt{T})
{
	$interactiveTargetMode = 1;
}




# parse instrument, run name and flowcell from cwd

my $path = `pwd`;
chomp $path;

if ($path =~ /.*\/(.*?)\/(.*?)$/)
{
	$instrument = $1;
	$runName = $2;
}
else
{
	die "Couldn't parse $path\n";
}

my @splitLib;

# get file names and parse lane, index from *.bam
my %bamHash;
open (FILE, $bamFile) or usage("Couldn't open $bamFile!\n");
while ($line = <FILE>)
{
	chomp $line;

	if ($line =~ /.*Project_(.*?)\/Sample_(.*?)\/.*?_(.*?)_L00(.)_.*_novoalign\.sam\.sorted\.bam/)
	{
		$bamHash{$line}{"project"} = $1;
		$bamHash{$line}{"library"} = $2;
		$bamHash{$line}{"barcode"} = $3;
		$bamHash{$line}{"lane"} = $4;
	}
	else
	{
		die "Couldn't parse $line.\n";
	}

	@splitLib = split(/_/, $bamHash{$line}{"library"});
	$bamHash{$line}{"sample"} = $splitLib[0] . $splitLib[1] . $splitLib[3];
}



# determine target files
my %projectHash;

if ($interactiveTargetMode ==0)
{
	for my $ius (keys %bamHash)
	{
		$bamHash{$ius}{"target"} = $bedFile;
	}
}
else
{
	# build project hash
	for my $ius (keys %bamHash)
	{
		if (exists $projectHash{$bamHash{$ius}{"project"}})
		{
			$projectHash{$bamHash{$ius}{"project"}}++;
		}
		else
		{
			$projectHash{$bamHash{$ius}{"project"}} = 1;
		}
	}

	warn "\nEntering targets in interactive mode.\n";
	warn "Common targets are:\n";
	warn "\thg19 whole genome: /oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed\n";
	warn "\tNimblegen: /oicr/data/reference/genomes/homo_sapiens/Nimblegen/2.1M_Human_Exome_Annotation_21191/hg19/080904_ccds_exome_rebalfocus_hg19/processed/2.1M_Human_Exome.bed.txt\n";
	warn "\tAgilent: /oicr/data/genomes/homo_sapiens/Agilent/SureSelect_Whole_Exome_ICGC_Sanger/GRCh37hg19/sanger.exons.bed.hg19\n";
	warn "\tTruSeq: /oicr/data/reference/genomes/homo_sapiens_mc/TruSeq/TruSeq-Exome-Targeted-Regions-BED-file\n";

	warn "\n";
	for my $p (keys %projectHash)
	{
		warn "Please enter a target for the $projectHash{$p} bam files from project $p:\n";
		$bedFile = <>;
		chomp $bedFile;
		$projectHash{$p} = $bedFile;
	}
	warn "\nTargets aquired.\n";

	for my $ius (keys %bamHash)
	{
		$bamHash{$ius}{"target"} = $projectHash{$bamHash{$ius}{"project"}};
	}
}




# print commands

my $fileName;
my $sequenceFiles;
my %jsonHash;
my $jsonText;

my $jsonOutFile;
my $hold_jid = "";

unless (-d "./jsonReport")
{
	mkdir "./jsonReport";
}

open (COMMANDS, ">./jsonReport/goJson") or die "Couldn't open ./jsonReport/goJson\n";

for my $ius (sort keys %bamHash)
{
	if ($ius =~ /(.*)\/.*?/)
	{
		$path = $1;
	}

	%jsonHash = ();

	$jsonHash{"run name"} = $runName;
	$jsonHash{"instrument"} = $instrument;
	if (exists $bamHash{$ius}{"barcode"})
	{
		$jsonHash{"barcode"} = $bamHash{$ius}{"barcode"};
	}
	$jsonHash{"library"} = $bamHash{$ius}{"library"};
	$jsonHash{"sample"} = $bamHash{$ius}{"sample"};
	$jsonHash{"lane"} = $bamHash{$ius}{"lane"};

	$jsonText = encode_json(\%jsonHash);
	$jsonText =~ s/\"/\\\"/g;		# escaping \s for make safe qsub

	if (exists $bamHash{$ius}{"barcode"})
	{
		$jsonOutFile = "$jsonHash{\"run name\"}_$jsonHash{\"lane\"}_$jsonHash{\"barcode\"}_$jsonHash{\"library\"}_novoalign.json";
	}
	else
	{
		$jsonOutFile = "$jsonHash{\"run name\"}_$jsonHash{\"lane\"}_$jsonHash{\"library\"}_novoalign.json";
	}

	print COMMANDS "qsub -cwd -b y -V -N stat$jsonHash{\"lane\"}_$jsonHash{\"run name\"} -e ${jsonOutFile}.log -o ${jsonOutFile}.log \"module load samtools; samtools view ../$ius | $scriptPath/samStats.pl -r $bamHash{$ius}{\"target\"} -j '$jsonText' > ../$path/$jsonOutFile\"\n";

	$hold_jid = "stat$jsonHash{\"lane\"}_$jsonHash{\"run name\"},$hold_jid";

	`cd ./jsonReport; ln -s ../$path/$jsonOutFile .; cd ..`;
}

my $user = `whoami`;
chomp $user;
my $cwd = `pwd`;
chomp $cwd;
if ($cwd =~ /.*archive\/(.*)/)
{
	$cwd = $1;
}

print COMMANDS "qsub -cwd -b y -V -N jsonReport$runName -hold_jid $hold_jid \"$scriptPath/jsonToGenericRunReport.pl *.json > ${runName}_report.html\"\n";
print COMMANDS "qsub -cwd -b y -e mail.log -o mail.log -N jsonEmail$runName -hold_jid jsonReport$runName \"mail -s ${runName}_report_is_done ${user}\@oicr.on.ca < mail.txt\"\n";

close COMMANDS;

open (MAIL, ">./jsonReport/mail.txt") or die "Couldn't open ./jsonReport/mail.txt\n";
print MAIL "Report is here: http://hn1.hpc.oicr.on.ca/archive/${cwd}/jsonReport/${runName}_report.html\n";
close MAIL;

