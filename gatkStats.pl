#!/usr/bin/perl

use strict;
use warnings;

use JSON;

my $gatkDir = $ARGV[0];

my $snpVcf = `ls $gatkDir/*.snps.raw.filtered.vcf`;
chomp $snpVcf;

my $indelVcf = `ls $gatkDir/*.indels.raw.filtered.vcf`;
chomp $indelVcf;


my %jsonHash;
my $count;



$jsonHash{"snp date"} = (stat($snpVcf))[9];
$jsonHash{"indel date"} = (stat($indelVcf))[9];

$count = `grep -v "^#" $snpVcf | wc -l`;
chomp $count;
$jsonHash{"snp count"} = $count;

$count = `grep PASS $snpVcf | wc -l`;
chomp $count;
$jsonHash{"pass snp count"} = $count;

$count = `grep -v "^#" $indelVcf | wc -l`;
chomp $count;
$jsonHash{"indel count"} = $count;

$count = `grep PASS $indelVcf | wc -l`;
chomp $count;
$jsonHash{"pass indel count"} = $count;

warn $snpVcf . "\n";
if ($snpVcf =~ /\/\/(.*?)\.bam\.realigned\.recal\.bam\.snps\.raw\.filtered\.vcf/)
{
	$jsonHash{"sample"} = $1;
}


print encode_json(\%jsonHash);





