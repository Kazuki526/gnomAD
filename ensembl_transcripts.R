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


files = system("ls maf/all_maf",intern = T)
pick_enst = function(.file){
  .file = paste0("maf/all_maf/",.file)
  maf = read_tsv(.file) %>>%
    filter(str_length(ref)==1,str_length(alt)==1) %>>%
    tidyr::separate(HGVSc,into = c("ENST","nuc_change"),sep=":")%>>%
    count(SYMBOL,Gene,ENST,CANONICAL)
}
enst_list = tibble(file=files) %>>%
  mutate(enst_list = purrr::map(file,~pick_enst(.)))%>>%
  unnest()
write_df(enst_list,"maf/ensembl_transcripts/ENST_list.tsv")
