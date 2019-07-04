#!/usr/bin/perl
use warnings;
use strict;

my $gnomAD_version = "file/gnomad.exomes.r2.1.1.sites";

# this script execute in gnomAD exon raw file (by chromosome) dir
my @ls = `ls|grep $gnomAD_version`;
if(scalar(@ls)!= 24){die "ERROR::doing on wrong dir?\n";}

#make extract dir
mkdir "maf";
mkdir "maf/all_maf";
mkdir "maf/non_cancer_maf";
mkdir "maf/control_maf";

#vep extract colum
my @vep_focal = qw(SYMBOL Consequence IMPACT Gene HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons
				   STRAND CANONICAL SIFT PolyPhen LoF LoF_filter);
my @race = qw(afr sas amr eas nfe fin asj oth);
my $ac_col="AC\tAN\tnhomalt";
foreach my $race (@race){
		$ac_col .= "\tAC_$race\tAN_$race\tnhomalt_$race";
}
my @ac_col = split(/\t/,$ac_col);
my @chr = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);



use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(5);
foreach my $chr (@chr){
#fork start
		print "start chr$chr\n";
		$pm->start and next;
		&print_extracted_file("$chr");
#fork end 
		$pm->finish;
}
$pm->wait_all_children;
exit;

sub print_extracted_file( $ ){
		my $chr = $_[0];
		my $infile = "$gnomAD_version.$chr.vcf.bgz";
		open(VCF,"gunzip -c $infile|");
		open(ALL,"|gzip -c >all_maf/gnomAD_chr$chr.maf.gz");
		open(CAN,"|gzip -c >non_cancer_maf/non_cancer_chr$chr.maf.gz");
		open(CONT,"|gzip -c >control_maf/control_chr$chr.maf.gz");
		print ALL "chr\tposi\tref\talt\tfilter\t". join("\t",@vep_focal) ."\trf_tp_probability\tInbreedingCoeff\t$ac_col\n";
		print CAN "chr\tposi\tref\talt\tfilter\t". join("\t",@vep_focal) ."\trf_tp_probability\tInbreedingCoeff\t$ac_col\n";
		print CONT "chr\tposi\tref\talt\tfilter\t". join("\t",@vep_focal) ."\trf_tp_probability\tInbreedingCoeff\t$ac_col\n";
		my %vepcol=();
		while(<VCF>){
				if($_ =~ /^#/){
						if($_ =~ /^##INFO=<ID=vep/){
								chomp;
								%vepcol = &vep2colum($_);
						}
						next;
				}
				my @line = split(/\t/,);
				if($line[6] =~/AC0/){next;}
				my $vepout = &vepout($line[7],\%vepcol);
				if($vepout eq ""){next;}
				my $outbase = "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[6]$vepout";
				my %info = &info2hash($line[7]);
				my $allout ="$info{rf_tp_probability}\t$info{InbreedingCoeff}";
				foreach my $out_info (@ac_col){
						$allout .= "\t$info{$out_info}";
				}
				print ALL "$outbase\t$allout\n";
				if($info{non_cancer_AC}!=0){
						my $noncanout ="$info{rf_tp_probability}\t$info{InbreedingCoeff}";
						foreach my $out_info(@ac_col){
								my $out_info_nc=$info{"non_cancer_$out_info"};
								$noncanout .= "\t$out_info_nc";
						}
						print CAN "$outbase\t$noncanout\n";
				}
				if($info{controls_AC}!=0){
						my $contout ="$info{rf_tp_probability}\t$info{InbreedingCoeff}";
						foreach my $out_info(@ac_col){
								my $out_info_co=$info{"controls_$out_info"};
								$contout .= "\t$out_info_co";
						}
						print CONT "$outbase\t$contout\n";
				}
		}
		close VCF;
		close ALL;
		close CAN;
		close CONT;
}



sub vep2colum( $ ){
		my $info = $_[0];
		my $vepinfo="";
		if($info =~ /^##INFO=.*Format: (Allele.*LoF_info)\">$/){$vepinfo=$1;}else{die "ERROR::INFO vep is not matched\n";}
		my @vepinfo = split(/\|/,$vepinfo);
		my %vepcol=();
		for(my $i=0;$i < scalar(@vepinfo);$i++){
				$vepcol{$vepinfo[$i]}=$i;
		}
		foreach my $vep_focal(@vep_focal){
				if(!defined $vepcol{$vep_focal}){die "ERROR::not exist $vep_focal on INFO of vep\n";}
		}
		return(%vepcol);
}

sub vepout( $ $ ){
		my $info = $_[0];
		my %vepcol = %{$_[1]};
		my $vep_text;
		if($info =~ /;vep=([^;]+)$/){$vep_text=$1;}
		my @vep_text = split(/,/,$vep_text);
		my @vep = split(/\|/,$vep_text[0]);
		my $out="";
		if(($vep[$vepcol{IMPACT}] eq "MODIFIER") || ($vep[$vepcol{BIOTYPE}] ne "protein_coding")){
				return($out);
		}else{
				foreach my $vep_focal(@vep_focal){
						if(!defined $vep[$vepcol{$vep_focal}]){$out .="\t";
						}else{$out .= "\t$vep[$vepcol{$vep_focal}]";}
				}
				return($out);
		}
}
sub info2hash( $ ){
		my $info = $_[0];
		my @info = split(/;/,$info);
		my %out =();
		for(my $i=0;$i<scalar(@info);$i++){
				if($i<6){
						if($info[$i] =~ /^([^=]+)=(.+)$/){$out{$1}=$2;}else{die "ERROR::what info?? $info[$i]\n";}
				}else{
						if($info[$i] =~ /^([^=]+)=(\d+)$/){$out{$1}=$2;}
				}
		}
		return(%out);
}
		




sub info_check_all( $ ){
		my $info=$_[0];
		my @info = split(/;/,$info);
		my $out="";
		if($info[2] =~ /^rf_tp_probability=(.+)$/){$out.="$1\t";}else{die "ERROR::$info[2] is not rf_tp_probability?\n";}
		if($info[4] =~ /^InbreedingCoeff=(.+)$/){$out.="$1\t";}else{die "ERROR::$info[4] is not InbreedingCoeff?\n";}
		if($info[0] =~ /^AC=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[0] is not AC?\n";}
		if($info[1] =~ /^AN=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[1] is not AN?\n";}
		if($info[151] =~ /^nhomalt=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[151] is not nhomalt?\n";}
#afr
		if($info[70] =~ /^AC_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[70] is not AC_afr?\n";}
		if($info[71] =~ /^AN_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[71] is not AN_afr?\n";}
		if($info[72] =~ /^nhomalt_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[72] is not nhomalt_afr?\n";}
#sas
		if($info[103] =~ /^AC_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[103] is not AC_sas?\n";}
		if($info[104] =~ /^AN_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[104] is not AN_sas?\n";}
		if($info[105] =~ /^nhomalt_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[105] is not nhomalt_sas?\n";}
#amr
		if($info[133] =~ /^AC_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[133] is not AC_amr?\n";}
		if($info[134] =~ /^AN_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[134] is not AN_amr?\n";}
		if($info[135] =~ /^nhomalt_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[135] is not nhomalt_amr?\n";}
#eas
		if($info[148] =~ /^AC_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[148] is not AC_eas?\n";}
		if($info[149] =~ /^AN_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[149] is not AN_eas?\n";}
		if($info[150] =~ /^nhomalt_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[150] is not nhomalt_eas?\n";}
#nfe
		if($info[400] =~ /^AC_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[400] is not AC_nfe?\n";}
		if($info[401] =~ /^AN_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[401] is not AN_nfe?\n";}
		if($info[402] =~ /^nhomalt_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[402] is not nhomalt_nfe?\n";}
#fin
		if($info[425] =~ /^AC_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[425] is not AC_fin?\n";}
		if($info[426] =~ /^AN_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[426] is not AN_fin?\n";}
		if($info[427] =~ /^nhomalt_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[427] is not nhomalt_fin?\n";}
#asj
		if($info[473] =~ /^AC_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[473] is not AC_asj?\n";}
		if($info[474] =~ /^AN_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[474] is not AN_asj?\n";}
		if($info[475] =~ /^nhomalt_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[475] is not nhomalt_asj?\n";}
#oth
		if($info[497] =~ /^AC_oth=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[497] is not AC_oth?\n";}
		if($info[498] =~ /^AN_oth=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[498] is not AN_oth?\n";}
		if($info[499] =~ /^nhomalt_oth=(\d+)$/){$out.="$1";}else{die "ERROR::$info[499] is not nhomalt_oth?\n";}
		return($out);
}

sub info_check_noncancer( $ ){
		my $info=$_[0];
		my @info = split(/;/,$info);
		my $out="";
		if($info[2] =~ /^rf_tp_probability=(.+)$/){$out.="$1\t";}else{die "ERROR::$info[2] is not rf_tp_probability?\n";}
		if($info[4] =~ /^InbreedingCoeff=(.+)$/){$out.="$1\t";}else{die "ERROR::$info[4] is not InbreedingCoeff?\n";}
		if($info[485] =~ /^non_cancer_AC=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[485] is not non_cancer_AC?\n";}
		if($info[486] =~ /^non_cancer_AN=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[486] is not non_cancer_AN?\n";}
		if($info[487] =~ /^non_cancer_nhomalt=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[487] is not non_cancer_nhomalt?\n";}
#afr
		if($info[255] =~ /^non_cancer_AC_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[255] is not non_cancer_AC_afr?\n";}
		if($info[256] =~ /^non_cancer_AN_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[256] is not non_cancer_AN_afr?\n";}
		if($info[257] =~ /^non_cancer_nhomalt_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[257] is not non_cancer_nhomalt_afr?\n";}
#sas
		if($info[561] =~ /^non_cancer_AC_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[561] is not non_cancer_AC_sas?\n";}
		if($info[562] =~ /^non_cancer_AN_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[562] is not non_cancer_AN_sas?\n";}
		if($info[563] =~ /^non_cancer_nhomalt_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[563] is not non_cancer_nhomalt_sas?\n";}
#amr
		if($info[282] =~ /^non_cancer_AC_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[282] is not non_cancer_AC_amr?\n";}
		if($info[283] =~ /^non_cancer_AN_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[283] is not non_cancer_AN_amr?\n";}
		if($info[284] =~ /^non_cancer_nhomalt_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[284] is not non_cancer_nhomalt_amr?\n";}
#eas
		if($info[162] =~ /^non_cancer_AC_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[162] is not non_cancer_AC_eas?\n";}
		if($info[163] =~ /^non_cancer_AN_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[163] is not non_cancer_AN_eas?\n";}
		if($info[164] =~ /^non_cancer_nhomalt_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[164] is not non_cancer_nhomalt_eas?\n";}
#nfe
		if($info[464] =~ /^non_cancer_AC_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[464] is not non_cancer_AC_nfe?\n";}
		if($info[465] =~ /^non_cancer_AN_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[465] is not non_cancer_AN_nfe?\n";}
		if($info[466] =~ /^non_cancer_nhomalt_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[466] is not non_cancer_nhomalt_nfe?\n";}
#fin
		if($info[506] =~ /^non_cancer_AC_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[506] is not non_cancer_AC_fin?\n";}
		if($info[506] =~ /^non_cancer_AN_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[507] is not non_cancer_AN_fin?\n";}
		if($info[508] =~ /^non_cancer_nhomalt_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[508] is not non_cancer_nhomalt_fin?\n";}
#asj
		if($info[207] =~ /^non_cancer_AC_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[207] is not non_cancer_AC_asj?\n";}
		if($info[208] =~ /^non_cancer_AN_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[208] is not non_cancer_AN_asj?\n";}
		if($info[209] =~ /^non_cancer_nhomalt_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[209] is not non_cancer_nhomalt_asj?\n";}
#oth
		if($info[345] =~ /^non_cancer_AC_oth=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[345] is not non_cancer_AC_oth?\n";}
		if($info[346] =~ /^non_cancer_AN_oth=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[346] is not non_cancer_AN_oth?\n";}
		if($info[347] =~ /^non_cancer_nhomalt_oth=(\d+)$/){$out.="$1";}else{die "ERROR::$info[347] is not non_cancer_nhomalt_oth?\n";}
		return($out);
}

sub info_check_cont( $ ){
		my $info=$_[0];
		my @info = split(/;/,$info);
		my $out="";
		if($info[2] =~ /^rf_tp_probability=(.+)$/){$out.="$1\t";}else{die "ERROR::$info[2] is not rf_tp_probability?\n";}
		if($info[4] =~ /^InbreedingCoeff=(.+)$/){$out.="$1\t";}else{die "ERROR::$info[4] is not InbreedingCoeff?\n";}
		if($info[582] =~ /^controls_AC=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[582] is not controls_AC?\n";}
		if($info[583] =~ /^controls_AN=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[583] is not controls_AN?\n";}
		if($info[584] =~ /^controls_nhomalt=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[584] is not controls_nhomalt?\n";}
#afr
		if($info[109] =~ /^controls_AC_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[109] is not controls_AC_afr?\n";}
		if($info[110] =~ /^controls_AN_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[110] is not controls_AN_afr?\n";}
		if($info[111] =~ /^controls_nhomalt_afr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[111] is not controls_nhomalt_afr?\n";}
#sas
		if($info[382] =~ /^controls_AC_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[382] is not controls_AC_sas?\n";}
		if($info[383] =~ /^controls_AN_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[383] is not controls_AN_sas?\n";}
		if($info[384] =~ /^controls_nhomalt_sas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[384] is not controls_nhomalt_sas?\n";}
#amr
		if($info[446] =~ /^controls_AC_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[446] is not controls_AC_amr?\n";}
		if($info[447] =~ /^controls_AN_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[447] is not controls_AN_amr?\n";}
		if($info[448] =~ /^controls_nhomalt_amr=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[448] is not controls_nhomalt_amr?\n";}
#eas
		if($info[373] =~ /^controls_AC_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[373] is not controls_AC_eas?\n";}
		if($info[374] =~ /^controls_AN_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[374] is not controls_AN_eas?\n";}
		if($info[375] =~ /^controls_nhomalt_eas=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[375] is not controls_nhomalt_eas?\n";}
#nfe
		if($info[213] =~ /^controls_AC_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[213] is not controls_AC_nfe?\n";}
		if($info[214] =~ /^controls_AN_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[214] is not controls_AN_nfe?\n";}
		if($info[215] =~ /^controls_nhomalt_nfe=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[215] is not controls_nhomalt_nfe?\n";}
#fin
		if($info[306] =~ /^controls_AC_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[306] is not controls_AC_fin?\n";}
		if($info[307] =~ /^controls_AN_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[307] is not controls_AN_fin?\n";}
		if($info[308] =~ /^controls_nhomalt_fin=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[308] is not controls_nhomalt_fin?\n";}
#asj
		if($info[219] =~ /^controls_AC_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[219] is not controls_AC_asj?\n";}
		if($info[220] =~ /^controls_AN_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[220] is not controls_AN_asj?\n";}
		if($info[221] =~ /^controls_nhomalt_asj=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[221] is not controls_nhomalt_asj?\n";}
#oth
		if($info[394] =~ /^controls_AC_oth=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[394] is not controls_AC_oth?\n";}
		if($info[395] =~ /^controls_AN_oth=(\d+)$/){$out.="$1\t";}else{die "ERROR::$info[395] is not controls_AN_oth?\n";}
		if($info[396] =~ /^controls_nhomalt_oth=(\d+)$/){$out.="$1";}else{die "ERROR::$info[396] is not controls_nhomalt_oth?\n";}
		return($out);
}
