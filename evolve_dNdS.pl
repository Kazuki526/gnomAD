#!/usr/bin/perl
use strict;
use warnings;


my $genome = "/Volumes/areca42TB/GRCh37-lite/GRCh37-lite.fa";
-e $genome or die "ERROR::$genome not exist\n";
my $file = "maf/ensembl_transcripts/ENST_list.tsv";
my %enst =();
open(ELS,"$file") or die "ERROR::$file is not exist\n";
<ELS>;
while(<ELS>){
		chomp;
		my @line=split(/\t/,);
		my ($enst_id,$chr) = ("","");
		if($line[0] =~ /_chr([0-9XY]+)\./){$chr=$1;}
		if($line[3] =~ /^(ENST\d+)\.\d+$/){$enst_id=$1;}else{die "ERROR::what ENST?\n$_\n";next;}
		$enst{$enst_id}{maf_id}=$line[3];
		$enst{$enst_id}{line}  =join("\t",@line[1..5]);
}
close ELS;

print "read gff\n";
my $gff = "/Volumes/areca42TB/GRCh37-lite/ensembl.GRCh37.gff3.gz";
open(GFF,"gunzip -c $gff|");
my $step=0;
while(<GFF>){
		if($_=~/^#/){$step=0;next;
		}elsif($step==-1){next;}
		chomp;
		my @line = split(/\t/,);
		if($step==0){
				if($line[2] eq "gene"){$step=1;}else{$step=-1;}
		}elsif($step==1){
				if($line[2] !~ "transcript"){$step=2;next;}
				if($line[8] =~ /^ID=transcript:(ENST\d+);/){
						if(defined$enst{$1}){
								$enst{$1}{strand}=$line[6];
								$step=3;
						}else{$step=2;}
				}else{die "ERROR::next of gene is not transcript\n$_\n";}
		}elsif($step==2){
				if($line[2] !~ "transcript"){next;}
				if($line[8] =~ /^ID=transcript:(ENST\d+);/){
						if(defined$enst{$1}){
								$enst{$1}{strand}=$line[6];
								$step=3;
						}else{$step=2;}
				}
		}elsif($step==3){
				if(($line[2] =~ "transcript")){
						if($line[8] =~ /^ID=transcript:(ENST\d+);/){
								if(defined$enst{$1}){
										$enst{$1}{strand}=$line[6];
										$step=3;
								}else{$step=2;}
						}else{die "ERROR::what gff line? transcript??\n$_\n";}
				}elsif($line[2] eq "exon"){
						if($line[8] =~ /^Parent=transcript:(ENST\d+);/){$enst{$1}{exon}++;
						}else{die "ERROR::what gff line? exon??\n$_\n";}
				}elsif($line[2] eq "CDS"){
						if($line[8] =~ /;Parent=transcript:(ENST\d+);/){
								$enst{$1}{cds} .= "$line[0]:$line[3]-$line[4],";
						}else{die "ERROR::what gff line? CDS??\n$_\n";}
				}
		}elsif($line[2] !~ /UTR/){die "ERROR::$_\n";}
}
close GFF;

print "calculate number of site\n";
my %nos = &calculate_codon(); #nomber of site
open(OUT,">maf/ensembl_transcripts/ensembl_enst_dNdS.tsv");
print OUT "chr\tSYMBOL\tGene\tENST\tCANONICAL\tn\tsynonymous_site\tnonsynonymous_site\ttruncating_site\tsplice_site\tfilter\n";
foreach my $enst_id(sort (keys %enst)){
		if(!defined$enst{$enst_id}{strand}){print "WARNNING::$enst_id:$enst{$enst_id}{strand}:$enst{$enst_id}{cds} has no cds info!!\n";next;}
		my @cds = split(/,/,$enst{$enst_id}{cds});
		my ($cds_seq,$chr)=("","");
		foreach my $region (@cds){
				my($start,$end);
				if($region =~ /^([0-9XY]+):(\d+)-(\d+)$/){($chr,$start,$end)=($1,$2,$3);}
				my $fa = `samtools faidx $genome $region`;
				$fa=~s/^>[0-9XY]+:\d+-\d+\n//;$fa=~s/\n//g;
				$cds_seq.=$fa;
		}
		if(length($cds_seq) %3 !=0){
				print "WARNNING::length of $enst_id is not a multiple of 3\n";
				print OUT "$chr\t$enst{$enst_id}{line}\t0\t0\t0\t0\tnot3\n";
		}else{
				if($enst{$enst_id}{strand} eq "-"){
						$cds_seq =~ tr/ACGTacgt/TGCAtgca/;
						$cds_seq = reverse($cds_seq);
				}
				my $splice = ($enst{$enst_id}{exon}-1)*4;
				my ($syn,$nonsyn,$trunc,$error) = &seq2nos($cds_seq,"$enst_id$enst{$enst_id}{strand}");
				print OUT "$chr\t$enst{$enst_id}{line}\t$syn\t$nonsyn\t$trunc\t$splice\t$error\n";
		}
}
close OUT;



sub calculate_codon( ){
		my(%codon2prot) = (
			'TCA' => 'S',    # Serine
			'TCC' => 'S',    # Serine
			'TCG' => 'S',    # Serine
			'TCT' => 'S',    # Serine
			'TTC' => 'F',    # Phenylalanine
			'TTT' => 'F',    # Phenylalanine
			'TTA' => 'L',    # Leucine
			'TTG' => 'L',    # Leucine
			'TAC' => 'Y',    # Tyrosine
			'TAT' => 'Y',    # Tyrosine
			'TAA' => '*',    # Stop
			'TAG' => '*',    # Stop
			'TGC' => 'C',    # Cysteine
			'TGT' => 'C',    # Cysteine
			'TGA' => '*',    # Stop
			'TGG' => 'W',    # Tryptophan
			'CTA' => 'L',    # Leucine
			'CTC' => 'L',    # Leucine
			'CTG' => 'L',    # Leucine
			'CTT' => 'L',    # Leucine
			'CCA' => 'P',    # Proline
			'CCC' => 'P',    # Proline
			'CCG' => 'P',    # Proline
			'CCT' => 'P',    # Proline
			'CAC' => 'H',    # Histidine
			'CAT' => 'H',    # Histidine
			'CAA' => 'Q',    # Glutamine
			'CAG' => 'Q',    # Glutamine
			'CGA' => 'R',    # Arginine
			'CGC' => 'R',    # Arginine
			'CGG' => 'R',    # Arginine
			'CGT' => 'R',    # Arginine
			'ATA' => 'I',    # Isoleucine
			'ATC' => 'I',    # Isoleucine
			'ATT' => 'I',    # Isoleucine
			'ATG' => 'M',    # Methionine
			'ACA' => 'T',    # Threonine
			'ACC' => 'T',    # Threonine
			'ACG' => 'T',    # Threonine
			'ACT' => 'T',    # Threonine
			'AAC' => 'N',    # Asparagine
			'AAT' => 'N',    # Asparagine
			'AAA' => 'K',    # Lysine
			'AAG' => 'K',    # Lysine
			'AGC' => 'S',    # Serine
			'AGT' => 'S',    # Serine
			'AGA' => 'R',    # Arginine
			'AGG' => 'R',    # Arginine
			'GTA' => 'V',    # Valine
			'GTC' => 'V',    # Valine
			'GTG' => 'V',    # Valine
			'GTT' => 'V',    # Valine
			'GCA' => 'A',    # Alanine
			'GCC' => 'A',    # Alanine
			'GCG' => 'A',    # Alanine
			'GCT' => 'A',    # Alanine
			'GAC' => 'D',    # Aspartic Acid
			'GAT' => 'D',    # Aspartic Acid
			'GAA' => 'E',    # Glutamic Acid
			'GAG' => 'E',    # Glutamic Acid
			'GGA' => 'G',    # Glycine
			'GGC' => 'G',    # Glycine
			'GGG' => 'G',    # Glycine
			'GGT' => 'G',    # Glycine
						) ;
		my $nos =(); #number of site
		foreach my $codon (keys %codon2prot){
				my @nuc = qw(A T G C);
				my ($syn,$nsyn,$trun)=(0,0,0);
				for(my$posi=0;$posi<3;$posi++){
						foreach my $nuc(@nuc){
								my $mut_codon=$codon;
								substr($mut_codon,$posi,1,$nuc);
								if($codon2prot{$codon} eq $codon2prot{$mut_codon}){$syn++;
								}elsif(($codon2prot{$codon} eq '*')||($codon2prot{$mut_codon} eq '*')){$trun++;
								}else{$nsyn++;}
						}
						$syn--;
				}
				$nos{$codon}="$syn,$nsyn,$trun";
		}
		return(%nos);
}

sub seq2nos( $ $ ){
		my $seq = uc $_[0];
		my $error="pass";
		my($syn,$nsyn,$trun,$start_miss)=(0,0,0,0);
		if(substr($seq,0,3) eq "ATG"){
				$trun+=9;
		}else{
				$start_miss=9;$error="start";
				print "WARNING:$_[1] is not start with ATG\n";
		} #Consequence is coding_sequence_variant, so donot count
		my $stop=substr($seq,-3);		
		if(($stop ne "TAA")&&($stop ne "TAG")&&($stop ne "TGA")){
				$error="stop";
				print "WARNING:$_[1] is not stop with stop codon\n";
		}
		for(my $posi=3;$posi<length($seq);$posi+=3){
				my $codon = substr($seq,$posi,3);
				my($s,$n,$t)= split(/,/,$nos{$codon}) or die "ERROR:$codon\n";
				($syn,$nsyn,$trun)=($syn+$s,$nsyn+$n,$trun+$t);
		}
		if($syn+$nsyn+$trun+$start_miss != length($seq)*3){die "ERROR:of calculate number of site\n$seq\nsyn:$syn nonsyn:$nsyn truncating:$trun\n";}
		return($syn/3,$nsyn/3,$trun/3,$error);
}
