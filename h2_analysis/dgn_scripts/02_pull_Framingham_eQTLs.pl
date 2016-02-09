#!/usr/bin/perl
use warnings;
use strict;

###pull Framingham eQTL SNP lists at various thresholds: p<0.0001 (entire file), FDR<0.05
###For now, only include hapmap2 SNPs because that's what we selected from the DGN imputation


my $framing = '/group/im-lab/nas40t2/haky/Data/Transcriptome/eQTL-results/Framingham/eqtl-gene-1000g-peer-validate_c1-23_c28.csv.gz';
#my $framing = 'test.csv.gz';
my $hapmap = '/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list'; ##from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz

open(OUTALL, ">perl.Framingham_eqtl-gene_p0.0001_hapmapSnpsCEU.rsIDlist");
open(OUTFDR, ">perl.Framingham_eqtl-gene_fdr0.05_hapmapSnpsCEU.rsIDlist");

open(HAP, "$hapmap");
my %haplist;
while(<HAP>){
    chomp;
    my ($snp) = split(/\n/);
    $haplist{$snp} = 1;
}

system("zcat $framing > tmp.framing");

open(A, "tmp.framing");

while(<A>){
    chomp;
    my @array = split(/,/);
    my $snp = $array[3];
    if($snp eq 'Rs_ID'){
        next;
    }
    elsif(defined($haplist{$snp})){
        print OUTALL "$snp\n";
        my $logfdr = $array[22];
        my $fdr = 10**$logfdr;
        if($fdr < 0.05){
            print OUTFDR "$snp\n";
	}
    }
}
    

system("rm tmp.framing");
