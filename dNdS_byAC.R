library(tidyverse)
library(pipeR)
setwd('/Volumes/areca42TB/gnomAD/')
write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}
site_num = read_tsv("maf/ensembl_transcripts/ensembl_enst_dNdS.tsv")
classify_consequence = function(.data) {
  mutate(.data,
         mutype= dplyr::recode(Consequence,
                               `incomplete_terminal_codon_variant&coding_sequence_variant`='truncating',
                               missense_variant = 'nonsynonymous',
                               stop_gained = 'truncating',
                               stop_lost = 'truncating',
                               stop_retained_variant = 'synonymous',
                               synonymous_variant = 'synonymous',
                               splice_acceptor_variant = 'splice',
                               splice_donor_variant = 'splice',
                               start_lost = 'truncating',
                               `missense_variant&splice_region_variant`='nonsynonymous',
                               `splice_acceptor_variant&stop_gained`='out',
                               `splice_acceptor_variant&synonymous_variant`='out',
                               `splice_donor_variant&missense_variant`='out',
                               `splice_region_variant&3_prime_UTR_variant`='out',
                               `splice_region_variant&5_prime_UTR_variant`='out',
                               `splice_region_variant&coding_sequence_variant`='out',
                               `splice_region_variant&intron_variant`='out',
                               `splice_region_variant&stop_retained_variant`='synonymous',
                               `splice_region_variant&synonymous_variant`='synonymous',
                               `start_lost&splice_region_variant`='truncating',
                               `stop_gained&splice_region_variant`='truncating',
                               `stop_gained&start_lost`='truncating',
                               `stop_lost&splice_region_variant`='truncating'))
}
##############################################  all cases( gnomAD ) ###########################################
files = system("ls maf/all_maf",intern = T)
read_maf = function(.file){
  .file = paste0("maf/all_maf/",.file)
  read_tsv(.file) %>>%
    filter(str_length(ref)==1,str_length(alt)==1) %>>%
    tidyr::separate(HGVSc,into = c("ENST","nuc_change"),sep=":")%>>%
    mutate(chr=as.character(chr))
}
all_maf = tibble(file=files) %>>%
  mutate(maf = purrr::map(file,~read_maf(.)))%>>%
  unnest()
all_maf =classify_consequence(all_maf)%>>%filter(mutype != "out")%>>%dplyr::select(-file)

##nonfinish european
nfe_enst_list =all_maf %>>%
  filter(AC_nfe>0,AN_nfe>max(AN_nfe)*0.95,filter=="PASS")%>>%
  dplyr::select(chr,posi,ref,alt,SYMBOL,ENST,mutype,AC_nfe,AN_nfe,nhomalt_nfe) %>>%
  left_join(site_num)%>>%filter(filter=="pass")%>>%
  dplyr::select(-n)%>>%
  count(SYMBOL,ENST,synonymous_site,nonsynonymous_site,truncating_site,splice_site)
all_site_num = tibble(
  synonymous = sum(nfe_enst_list$synonymous_site),
  nonsynonymous = sum(nfe_enst_list$nonsynonymous_site),
  truncating = sum(nfe_enst_list$truncating_site),
  splice = sum(nfe_enst_list$splice_site)
)%>>%tidyr::gather(mutype,num_of_site)

mu_num_by_AC = all_maf %>>%
  filter(AC_nfe>0,AC_nfe<50,AN_nfe>max(AN_nfe)*0.95,filter=="PASS") %>>%
  inner_join(nfe_enst_list) %>>%
  count(mutype,AC_nfe) %>>%
  left_join(all_site_num) %>>%
  dplyr::rename(n=nn)%>>%
  mutate(dN=n/num_of_site)

#site_freq_spectra = 
mu_num_by_AC %>>%
  ggplot()+
  geom_bar(aes(x=AC_nfe,y=n),stat = "identity")+
  facet_wrap(~ mutype, scales="free")+
  scale_y_log10()
#dNdS plot
mu_num_by_AC %>>%
  filter(mutype!="synonymous")%>>%
  left_join(mu_num_by_AC%>>%
              filter(mutype=="synonymous")%>>%
              dplyr::select(AC_nfe,dN)%>>%dplyr::rename(dS=dN))%>>%
  mutate(dNdS=dN/dS)%>>%
  filter(dNdS<1)%>>%
  ggplot()+
  geom_point(aes(x=AC_nfe,y=dNdS))+
  facet_wrap(~mutype)+
  xlab(paste0("control AC_nfe (AN_nfe=",max(all_maf$AN_nfe)," )"))

### bar plot ###
AF_group = tibble(AF_start=c(0, 1, 3, 5, 10, 30, 100, 1000),
                  AF_end  =c(1, 3, 5,10, 30,100,1000,50000))
AF_grouping = function(AF){
  ifelse(AF<10**-5,0,
         ifelse(AF<2*10**-5,1,
                ifelse(AF<5*10**-5,3,
                       ifelse(AF<10**-4,5,
                              ifelse(AF<3*10**-4,10,
                                     ifelse(AF<10**-3,30,
                                            ifelse(AF<10**-2,100,1000)))))))
}
mut_AF_group =all_maf %>>%
  filter(AC_nfe>0,AN_nfe>max(AN_nfe)*0.95,filter=="PASS") %>>%
  inner_join(nfe_enst_list) %>>%
  mutate(AF=AC_nfe/AN_nfe) %>>%mutate(AF_start=AF_grouping(AF))%>>%
  count(mutype,AF_start)%>>%
  left_join(all_site_num) %>>%
  mutate(dN=nn/num_of_site)

mut_AF_group %>>%
  filter(mutype != "synonymous") %>>%
  left_join(mut_AF_group%>>%
              filter(mutype=="synonymous")%>>%
              dplyr::select(AF_start,dN)%>>%dplyr::rename(dS=dN))%>>%
  mutate(dNdS=dN/dS)%>>%
  left_join(AF_group)%>>%
  mutate(AF_start=AF_start/10**5,AF_end=AF_end/10**5)%>>%
  ggplot()+
  geom_rect(aes(xmin=AF_start,xmax=AF_end,ymin=0,ymax=dNdS))+
  facet_wrap(~mutype)+
  scale_x_log10()+
  labs(x="AF",y="dNdS")


###### nonsynonymous split by PolyPhen & SIFT #######
#PolyPhen
all_maf %>>%
  filter(AC_nfe>0,AN_nfe>max(AN_nfe)*0.95,filter=="PASS",mutype=="nonsynonymous") %>>%
  inner_join(nfe_enst_list)%>>%
  mutate(AF=AC_nfe/AN_nfe) %>>%mutate(AF_start=AF_grouping(AF))%>>%
  mutate(PolyPhen = str_extract(PolyPhen,"^[A-Za-z_]+")) %>>%
  count(PolyPhen,AF_start)%>>%
  mutate(dN=nn/sum(nfe_enst_list$nonsynonymous_site))%>>%
  filter(!is.na(PolyPhen),PolyPhen!="unknown")%>>%
  left_join(mut_AF_group%>>%
              filter(mutype=="synonymous")%>>%
              dplyr::select(AF_start,dN)%>>%dplyr::rename(dS=dN))%>>%
  mutate(dNdS=dN/dS)%>>%
  left_join(AF_group)%>>%
  mutate(AF_start=AF_start/10**5,AF_end=AF_end/10**5)%>>%
  ggplot()+
  geom_rect(aes(xmin=AF_start,xmax=AF_end,ymin=0,ymax=dNdS))+
  facet_wrap(~PolyPhen)+
  scale_x_log10()+
  labs(x="AF",y="dNdS")

#number of benign, possibly_damaging or probably_damaging siteが計算できないため無理


 #################################################### control only ###############################
files_cont = system("ls maf/control_maf",intern = T)
read_maf = function(.file){
  .file = paste0("maf/control_maf/",.file)
  read_tsv(.file) %>>%
    filter(str_length(ref)==1,str_length(alt)==1) %>>%
    tidyr::separate(HGVSc,into = c("ENST","nuc_change"),sep=":")%>>%
    mutate(chr=as.character(chr))
}
cont_maf = tibble(file=files_cont) %>>%
  mutate(maf = purrr::map(file,~read_maf(.)))%>>%
  unnest()
cont_maf =classify_consequence(cont_maf)%>>%filter(mutype != "out")%>>%dplyr::select(-file)

##nonfinish european
mu_num_by_AC_cont = cont_maf %>>%
  filter(AC_nfe>0,AC<50,AN_nfe>max(AN_nfe)*0.95,filter=="PASS") %>>%
  inner_join(nfe_enst_list) %>>%
  count(mutype,AC_nfe) %>>%
  left_join(all_site_num) %>>%
  dplyr::rename(n=nn)%>>%
  mutate(dN=n/num_of_site)
#dNdS plot
mu_num_by_AC_cont %>>%
  filter(mutype!="synonymous")%>>%
  left_join(mu_num_by_AC_cont%>>%
              filter(mutype=="synonymous")%>>%
              dplyr::select(AC_nfe,dN)%>>%dplyr::rename(dS=dN))%>>%
  mutate(dNdS=dN/dS)%>>%
  filter(dNdS<1)%>>%
  ggplot()+
  geom_point(aes(x=AC_nfe,y=dNdS))+
  facet_wrap(~mutype)+
  xlab(paste0("control AC_nfe (AN_nfe=",max(cont_maf$AN_nfe)," )"))
