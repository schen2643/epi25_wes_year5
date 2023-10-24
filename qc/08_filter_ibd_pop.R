setwd("~/qc_files")

library(tidyverse)
library(data.table)

# add arguments
pop <- "eur"

# Read in data
f.ibd <- paste0("gsutil cat gs://epi25/year5/qc_files/08_ibd_",pop,".tsv")
ibd <- fread(f.ibd, sep='\t', stringsAsFactors=TRUE, header=TRUE, data.table=FALSE)
head(ibd)
ibd <- ibd[,-c(3:4)]

# Limit to IBD > 0.2 (related pairs)
pihat <- 0.2
ibd.rel <- ibd %>% 
  filter(PI_HAT > pihat
         )

# Add relatedness info
ibd.rel$inferred.rel <- "Other relatedness"
ibd.rel$inferred.rel <- ifelse(ibd.rel$ibd_Z0<0.11 & ibd.rel$ibd_Z1<0.11, 
                               "Duplicates/MZ-twins", ibd.rel$inferred.rel)
ibd.rel$inferred.rel <- ifelse(ibd.rel$ibd_Z0<0.1 & ibd.rel$ibd_Z1>0.75,
                               "Parent-offspring", ibd.rel$inferred.rel)
ibd.rel$inferred.rel <- ifelse(ibd.rel$ibd_Z0>0.07 & ibd.rel$ibd_Z0<0.375 &
                                 ibd.rel$ibd_Z1>0.25 & ibd.rel$ibd_Z1<0.75,
                               "Siblings", ibd.rel$inferred.rel)
ibd.rel$inferred.rel = ifelse(ibd.rel$PI_HAT <= pihat,
                              "Unrelated", ibd.rel$inferred.rel)
table(ibd.rel$inferred.rel)

# Add case-control info
phe <- fread("../misc_data/epi25_year1-5_pheno_n72763_20210912.tsv", header=T, sep="\t", stringsAsFactors = F, data.table=F)
ibd.rel$i.cc <- phe[match(ibd.rel$i, phe$sample),"case_control"]
ibd.rel$j.cc <- phe[match(ibd.rel$j, phe$sample),"case_control"]
ibd.rel$i.epi <- phe[match(ibd.rel$i, phe$sample),"epilepsy_type"]
ibd.rel$j.epi <- phe[match(ibd.rel$j, phe$sample),"epilepsy_type"]

ibd.rel$i.cohort <- phe[match(ibd.rel$i, phe$sample),"cohort"]
ibd.rel$j.cohort <- phe[match(ibd.rel$j, phe$sample),"cohort"]
table(ibd.rel$i.cohort)
table(ibd.rel$j.cohort)

samples <- names(sort(table(unlist(ibd.rel[,c('i', 'j')])), decreasing=TRUE))
# sorted by occurrence of sample IDs
# number of unique sample IDs: length(unique(unlist(ibd.rel[,c('i', 'j')])))
removed <- character()

#####
# Before IBD pruning:
# Examine the case-control status of these related pairs/samples
phe.rel <- phe %>% filter(sample %in% samples) %>% select(sample, case_control)
table(phe.rel[match(samples, phe.rel$sample),"case_control"]) 

for (row in 1:nrow(ibd.rel)) {
  
  i_sample <- as.character(ibd.rel[row, 'i'])
  j_sample <- as.character(ibd.rel[row, 'j'])
  
  i_index <- match(i_sample, samples)
  j_index <- match(j_sample, samples)
  
  if (is.na(i_index) | is.na(j_index)) { next }
  
  if (i_index <= j_index) {
    removed <- c(removed, i_sample)
    samples <- samples[samples != i_sample]
  } else {
    removed <- c(removed, j_sample)
    samples <- samples[samples != j_sample]
  } 
  
}
length(removed)
length(samples) 

#####
# After IBD pruning:
# To keep
table(phe.rel[match(samples, phe.rel$sample),"case_control"])  

# To remove
table(phe.rel[match(removed, phe.rel$sample),"case_control"]) 
write.table(removed, paste0("08_ibd_",pop,"_sample.remove.list"), quote=F, col.names=F, row.names=F)

names(ibd.rel) <- c("sample1","sample2","Z0","Z1","Z2","PI_HAT","ibs0","ibs1","ibs2",
                    "inferred.rel",
                    "sample1.cc","sample2.cc","sample1.epi","sample2.epi",
                    "sample1.cohort","sample2.cohort")
write.table(ibd.rel, paste0("08_ibd_",pop,"_sample_inferred_relatedness_pihat02.tsv"),
            quote=F, col.names=T, row.names=F, sep="\t")

phe.rel.rm <- phe[match(removed, phe$sample),c("sample","case_control","epilepsy_type","cohort")]
write.table(phe.rel.rm, paste0("08_ibd_",pop,"_sample_removed_pheno.tsv"), quote=F, col.names=T, row.names=F,
            sep="\t")

head(phe.rel.rm)
table(phe.rel.rm$case_control)
table(phe.rel.rm$epilepsy_type)

