setwd("~/qc_files")

library(tidyverse)
library(data.table)
library(ggsci)
library(RColorBrewer)
library(colorRamps)

pop <- "eur"

# Read in data
f.sqc <- paste0("gsutil cat gs://epi25/wes/qc_files/10_",pop,"_finalQC_sample.tsv")
dat.sqc <- fread(f.sqc, stringsAsFactors=FALSE, sep='\t', header=TRUE, data.table=FALSE)
dim(dat.sqc)
names(dat.sqc)[7] <- 'bait_set'

table(dat.sqc$cohort)
table(dat.sqc$epi25_site)
cohort_names <- names(table(dat.sqc$cohort))
epi25_site_names <- names(table(dat.sqc$epi25_site))[-1]
cohort_names[!(cohort_names %in% epi25_site_names)]
dat.sqc$cohort <- factor(dat.sqc$cohort, levels=c(names(table(dat.sqc$epi25_site))[-1],
                                                  cohort_names[!(cohort_names %in% epi25_site_names)]))
names(table(dat.sqc$cohort))
# dat.sqc$cohort <- factor(dat.sqc$cohort, levels=c(levels(dat.sqc$cohort)[c(1:34,36:53,59:60,
#                                                                            35,54:58,61:66)]))
# dat.sqc$cohort <- recode(dat.sqc$cohort, 
#                          `ALSPAC` = "ALSPAC_controls",
#                          `CCDG_MGB_Biobank` = "CCDG_controls",
#                          `Dalio Dutch Posthuma` = "Dutch_controls",
#                          `Dalio USA Smoller & Yolken` = "USA_controls",
#                          `Dalio_German_Reif` = "German_controls",
#                          `Epi PME` = "PME_cohort",
#                          `Epi4k` = "Epi4K_cohort",
#                          `GPC Caucasian` = "GPC_Caucasian_controls",
#                          `GPC Latino` = "GPC_Latino_controls",
#                          `IBC NIDDK CCDG Y2-3` = "NIDDK1_controls",
#                          `IBC NIDDK CCDG Y4` = "NIDDK2_controls",
#                          `IBD_FINRISK` = "FINRISK_controls",
#                          `MIGen_ExS_Leicester` = "MIGen_Leicester_controls",
#                          `UK/IRL Edinburgh + McQuillan` = "UK/IRL_controls",
#                          `Spalletta_controls` =   "Italian_controls",
#                          `Cheung_controls` = "HK_controls"
# )

batch=names(table(dat.sqc$batch))
dat.sqc$batch <- factor(dat.sqc$batch, levels=rev(c(sort(batch[grep("Epi25",batch)]),sort(batch[grep("Epi25",batch, invert=T)]))))
dat.sqc$case_control <- factor(dat.sqc$case_control, levels=rev(c("case","control")))
dat.sqc$case_control <- recode(dat.sqc$case_control,`case` = "Epilepsy", `control` = "Control")
dat.sqc$cohort <- factor(dat.sqc$cohort, levels=rev(levels(dat.sqc$cohort)))
dat.sqc$bait_set <- recode(dat.sqc$bait_set,`twist` = "Twist", `ice` = "Illumina")


# Setting colors:
col.cc <- c(brewer.pal(9,"YlOrRd")[5], "#5b3396")
col.bait <- pal_npg("nrc")(9)[c(3,1)]

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colorCount_cohort <- length(unique(dat.sqc$cohort)) # number of levels
col.cohort <- rev(getPalette(colorCount_cohort))

colorCount_batch <- length(unique(dat.sqc$batch)) # number of levels
col.batch <- rev(getPalette(colorCount_batch))

# Setting other variables
plots_dir <- "~/plots/"


####
# Examine qc metrics by cohorts/all cohorts
# two cohorts have only one sample (ITAUM3 and GPC_Latino_controls)
dat.sqc.cohort <- dat.sqc %>% 
  group_by(cohort) %>%
  summarise(
    # rTiTv
    avg_rTiTv = mean(rTiTv),
    sd_rTiTv = sd(rTiTv),
    upper_rTiTv = avg_rTiTv + 4*sd_rTiTv,
    lower_rTiTv = avg_rTiTv - 4*sd_rTiTv,
    nOutlier.rTiTv = sum(rTiTv > upper_rTiTv | rTiTv < lower_rTiTv, na.rm=T),
    # rHetHom
    avg_rHetHom = mean(rHetHom),sd_rHetHom = sd(rHetHom),
    upper_rHetHom = avg_rHetHom + 4*sd_rHetHom,
    lower_rHetHom = avg_rHetHom - 4*sd_rHetHom,
    nOutlier.rHetHom = sum(rHetHom > upper_rHetHom | rHetHom < lower_rHetHom, na.rm=T),
    # rInsertDelete
    avg_rInDel = mean(rIndel),sd_rInDel = sd(rIndel),
    upper_rInDel = avg_rInDel + 4*sd_rInDel,
    lower_rInDel = avg_rInDel - 4*sd_rInDel,
    nOutlier.rInDel = sum(rIndel > upper_rInDel | rIndel < lower_rInDel, na.rm=T)
  )
sum(dat.sqc.cohort$nOutlier.rTiTv)
sum(dat.sqc.cohort$nOutlier.rHetHom)
sum(dat.sqc.cohort$nOutlier.rInDel)


# Identify outliers samples by each metric and remove them
# rTiTv
outlier.rTiTv = character()
for(cohort in dat.sqc.cohort$cohort[dat.sqc.cohort$nOutlier.rTiTv >0]){
  print(cohort)
  # dat.sqc.tmp = dat.sqc[which(dat.sqc$epi25_cohort %in% cohort | dat.sqc$cohort %in% cohort),]
  dat.sqc.tmp = dat.sqc[which(dat.sqc$cohort %in% cohort),]
  outlier.index = which(
    dat.sqc.tmp$rTiTv < dat.sqc.cohort$lower_rTiTv[dat.sqc.cohort$cohort==cohort] |
      dat.sqc.tmp$rTiTv > dat.sqc.cohort$upper_rTiTv[dat.sqc.cohort$cohort==cohort]
  )
  outlier.rTiTv = c(outlier.rTiTv, dat.sqc.tmp$s[outlier.index])
  rm(dat.sqc.tmp)
}
table(dat.sqc[dat.sqc$s %in% outlier.rTiTv, "case_control"])

# rHetHom
outlier.rHetHom = character()
for(cohort in dat.sqc.cohort$cohort[dat.sqc.cohort$nOutlier.rHetHom >0]){
  print(cohort)
  dat.sqc.tmp = dat.sqc[which(dat.sqc$cohort %in% cohort),]
  outlier.index = which(
    dat.sqc.tmp$rHetHom < dat.sqc.cohort$lower_rHetHom[dat.sqc.cohort$cohort==cohort] | 
      dat.sqc.tmp$rHetHom > dat.sqc.cohort$upper_rHetHom[dat.sqc.cohort$cohort==cohort]
  )
  outlier.rHetHom = c(outlier.rHetHom, dat.sqc.tmp$s[outlier.index])
  rm(dat.sqc.tmp)
}
table(dat.sqc[dat.sqc$s %in% outlier.rHetHom, "case_control"])

# rInDel
outlier.rInDel = character()
for(cohort in dat.sqc.cohort$cohort[dat.sqc.cohort$nOutlier.rInDel >0]){
  print(cohort)
  dat.sqc.tmp = dat.sqc[which(dat.sqc$cohort %in% cohort),]
  outlier.index = which(
    dat.sqc.tmp$rIndel < dat.sqc.cohort$lower_rInDel[dat.sqc.cohort$cohort==cohort] | 
      dat.sqc.tmp$rIndel > dat.sqc.cohort$upper_rInDel[dat.sqc.cohort$cohort==cohort]
  )
  outlier.rInDel = c(outlier.rInDel, dat.sqc.tmp$s[outlier.index])
  rm(dat.sqc.tmp)
}
table(dat.sqc[dat.sqc$s %in% outlier.rInDel, "case_control"])

# Take intersection
outlier_samples <- unique(c(outlier.rTiTv, outlier.rHetHom, outlier.rInDel))
length(outlier_samples) #44
table(dat.sqc[dat.sqc$s %in% outlier_samples, "case_control"])


# Plots of each sample QC metric

##################
# 1. Call Rate
##################

summary(dat.sqc$call_rate)
table(dat.sqc[dat.sqc$call_rate < 0.98, "case_control"])

callrate_thresh <- 0.98
# Overall
ggplot(dat.sqc, aes(call_rate)) +
  geom_histogram(alpha=0.6, binwidth=0.001, size=0.4, color="black", fill="#4DBBD5B2") +
  xlab("Call rate") +
  ylab("Frequency") +
  # coord_cartesian(xlim=c(0.98,1)) +
  ggtitle("Final sample QC") +
  geom_vline(xintercept=callrate_thresh, linetype="dashed") +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_callRate.pdf"), width=7, height=4)


ggplot(dat.sqc, aes(cohort, call_rate)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=cohort), alpha=0.6, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Call rate",
       x="Sample cohort") +
  geom_hline(yintercept=callrate_thresh, linetype="dashed") +
  scale_color_manual(values = col.cohort) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_callRate_by_cohort.pdf"), width=9.5, height=8)


# By case-control
ggplot(dat.sqc, aes(case_control, call_rate)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=case_control), show.legend=FALSE, alpha=0.8) + 
  labs(title="Final sample QC",
       y="Call rate",
       x="Case/Control",
       color="Case/Control") +
  geom_hline(yintercept=callrate_thresh, linetype="dashed") +
  scale_color_manual(values = col.cc) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_callRate_by_casectrl.pdf"), width=5.5, height=4)


###################
# 2. TiTv ratio
###################

# Overall
summary(dat.sqc$rTiTv)
ggplot(dat.sqc, aes(rTiTv)) +
  geom_histogram(alpha=0.6, binwidth=0.005, color="black", fill="#4DBBD5B2", size=0.3) +
  xlab("Transition/Transversion ratio") +
  ylab("Frequency") +
  ggtitle("Final sample QC") +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rTiTv.pdf"), width=7, height=4)

# By case-control
ggplot(dat.sqc, aes(case_control, rTiTv)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=case_control), alpha=0.7, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Transition/Transversion ratio",
       x="Case/Control") +
  scale_color_manual(values = col.cc) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rTiTv_by_casectrl.pdf"), width=6.5, height=5)

# By all cohorts
rTiTv.outlier.color <- ifelse(dat.sqc$s %in% outlier.rTiTv, "red", "turquoise4")
ggplot(dat.sqc, aes(cohort, rTiTv)) +
  geom_boxplot(outlier.shape=NA, coef=0, color='grey50', fill='grey71',
               show.legend=FALSE, alpha=0.7, width=0.6) +
  geom_jitter(
    width=0.15, size=0.7, color=rTiTv.outlier.color, alpha=0.6, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Transition/Transversion ratio",
       x="Sample cohort") +
  # ylim(3.1, 3.53) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rTiTv_by_cohort.pdf"), width=8, height=10)


###################
#  3. Het/Hom ratio
###################
# Overall
summary(dat.sqc$rHetHom)
ggplot(dat.sqc, aes(rHetHom)) +
  geom_histogram(alpha=0.6, binwidth=0.015, color="black", fill="#4DBBD5B2", size=0.3) +
  xlab("Het/Hom ratio") +
  ylab("Frequency") +
  ggtitle("Final sample QC") +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rHetHom.pdf"), width=7, height=4)

# By case-control
ggplot(dat.sqc, aes(case_control, rHetHom)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=case_control), alpha=0.7, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Het/Hom ratio",
       x="Case/Control") +
  scale_color_manual(values = col.cc) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rHetHom_by_casectrl.pdf"), width=6.5, height=5)

# By all cohorts
rHetHom.outlier.color <- ifelse(dat.sqc$s %in% outlier.rHetHom, "red", "turquoise4")
ggplot(dat.sqc, aes(cohort, rHetHom)) +
  geom_boxplot(outlier.shape=NA, coef=0, color='grey50', fill='grey71',
               show.legend=FALSE, alpha=0.7, width=0.6) +
  geom_jitter(
    width=0.15, size=0.7, color=rHetHom.outlier.color, alpha=0.6, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Het/Hom ratio",
       x="All cohorts") +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rHetHom_by_cohort.pdf"), width=8, height=10)


###################
# 4. Insertion/Deletion ratio
###################
# Overall
summary(dat.sqc$rIndel)
ggplot(dat.sqc, aes(rIndel)) +
  geom_histogram(alpha=0.6, binwidth=0.015, color="black", fill="#4DBBD5B2", size=0.3) +
  xlab("Insertion/Deletion ratio") +
  ylab("Frequency") +
  ggtitle("Final sample QC") +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rInDel.pdf"), width=7, height=4)


# By case-control
ggplot(dat.sqc, aes(case_control, rIndel)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=case_control), alpha=0.7, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Insertion/Deletion ratio",
       x="Case/Control") +
  scale_color_manual(values = col.cc) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rInDel_by_casectrl.pdf"), width=6.5, height=5)


# By all cohorts
rInDel.outlier.color <- ifelse(dat.sqc$s %in% outlier.rInDel, "red", "turquoise4")
ggplot(dat.sqc, aes(cohort, rIndel)) +
  geom_boxplot(outlier.shape=NA, coef=0, color='grey50', fill='grey71',
               show.legend=FALSE, alpha=0.7, width=0.6) +
  geom_jitter(width=0.15, size=0.7, color=rInDel.outlier.color, alpha=0.6, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Insertion/Deletion ratio",
       x="All cohorts") +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_rInDel_by_cohort.pdf"), width=8, height=10)


################################
#  5. Singleton synonymous count
################################

### wirte nSYNsingleton to for later use as a covar
write.table(dat.sqc[,c("s", "n_syn_singleton")], paste0("10_",pop,"_nSYNsingleton.tsv"), col.names = T, quote=F, row.names = F)

# Overall
summary(dat.sqc$n_syn_singleton)
ggplot(dat.sqc, aes(n_syn_singleton)) +
  geom_histogram(alpha=0.6, binwidth=1, color="black", fill="#4DBBD5B2", size=0.3) +
  xlab("Number of synonymous singletons called per individual") +
  ylab("Frequency") +
  ggtitle("Final sample QC") +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_nSYNsingleton.pdf"), width=7, height=4)

# By case-control
ggplot(dat.sqc, aes(case_control, n_syn_singleton)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=case_control), alpha=0.7, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Number of synonymous singletons called per individual",
       x="Case/Control") +
  scale_color_manual(values = col.cc) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_nSYNsingleton_by_casectrl.pdf"), width=6.5, height=5)


# By all cohorts
ggplot(dat.sqc, aes(cohort, n_syn_singleton)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=cohort), alpha=0.6, show.legend=FALSE) + 
  labs(title="Final sample QC",
       y="Number of synonymous singletons called per individual",
       x="Sample cohort") +
  scale_color_manual(values = col.cohort) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_nSYNsingleton_by_cohort.pdf"), width=9.5, height=10)

# By all cohorts, ordered by mean value of syn. singleton count
cohort_avg.nSYNsingleton <- dat.sqc %>%
  group_by(cohort) %>%
  summarise(avg.nSYNsingleton = mean(n_syn_singleton))
cohort_avg.nSYNsingleton <- 
  cohort_avg.nSYNsingleton[order(cohort_avg.nSYNsingleton$avg.nSYNsingleton),]

dat.sqc$cohort <- factor(dat.sqc$cohort, levels=rev(cohort_avg.nSYNsingleton$cohort))
levels(dat.sqc$cohort)

# get nSYNsingleton mean of epi25 cohort and controls
table(dat.sqc$cohort, dat.sqc$epi25_site=="")
# by cohort
epi25_mean <- mean(dat.sqc$n_syn_singleton[dat.sqc$epi25_site!=""])
control_mean <- mean(dat.sqc$n_syn_singleton[dat.sqc$epi25_site==""])
print(epi25_mean)

levels(dat.sqc$cohort)
# EUR
# 1:5, 6, 7:24, 25, 26:30, 31, 32:35, 36, 37:40, 41, 42:43, 44:45, 46:49, 50, 51:59, 60, 61, 62:63, 64:66
cohort.color = c(rep("#5b3396",5),rep("#FD8D3C",1),
                 rep("#5b3396",18),rep("#FD8D3C",1),
                 rep("#5b3396",5),rep("#FD8D3C",1),
                 rep("#5b3396",4),rep("#FD8D3C",1),
                 rep("#5b3396",4),rep("#FD8D3C",1),
                 rep("#5b3396",2),rep("#FD8D3C",2),
                 rep("#5b3396",4),rep("#FD8D3C",1),
                 rep("#5b3396",9),rep("#FD8D3C",1),
                 rep("#5b3396",1),rep("#FD8D3C",2),
                 rep("#5b3396",3))
cohort.color = c(rep("#FD8D3C",13),rep("#5b3396",1),
                 rep("#FD8D3C",5))
# AFR
# 1:15, 16, 17:22, 23:24, 25:27, 28, 29, 30, 31, 32, 33:40
cohort.color = c(rep("#5b3396",15),rep("#FD8D3C",1),
                 rep("#5b3396",6),rep("#FD8D3C",2),
                 rep("#5b3396",3),rep("#FD8D3C",1),
                 rep("#5b3396",1),rep("#FD8D3C",1),
                 rep("#5b3396",1),rep("#FD8D3C",1),
                 rep("#5b3396",8))

# EAS
# 1:21, 22:24, 25:32, 33, 34:37, 38:39, 40
cohort.color = c(rep("#5b3396",21),rep("#FD8D3C",3),
                 rep("#5b3396",8),rep("#FD8D3C",1),
                 rep("#5b3396",4),rep("#FD8D3C",2),
                 rep("#5b3396",1))
# SAS
# 1:8, 9, 10:11, 12, 13:21, 22, 23:26, 27, 28, 29, 30:31, 32, 33:35, 36, 37:42, 43, 44
cohort.color = c(rep("#5b3396",8),rep("#FD8D3C",1),
                 rep("#5b3396",2),rep("#FD8D3C",1),
                 rep("#5b3396",9),rep("#FD8D3C",1),
                 rep("#5b3396",4),rep("#FD8D3C",1),
                 rep("#5b3396",1),rep("#FD8D3C",1),
                 rep("#5b3396",2),rep("#FD8D3C",1),
                 rep("#5b3396",3),rep("#FD8D3C",1),
                 rep("#5b3396",6),rep("#FD8D3C",1),
                 rep("#5b3396",1))

# FIN
# 1:7, 8, 9, 10:12, 13:14, 15:16, 17:19, 20, 21:22, 23, 24
cohort.color = c(rep("#5b3396",7),rep("#FD8D3C",1),
                 rep("#5b3396",1),rep("#FD8D3C",3),
                 rep("#5b3396",2),rep("#FD8D3C",2),
                 rep("#5b3396",3),rep("#FD8D3C",1),
                 rep("#5b3396",2),rep("#FD8D3C",1),
                 rep("#5b3396",1))
# AMR
# 1:7, 8:9, 10:19, 20, 21:23, 24, 25, 26, 27:30, 31, 32:33, 34, 35:36, 37:38, 39:40, 41, 42:45
cohort.color = c(rep("#5b3396",7),rep("#FD8D3C",2),
                 rep("#5b3396",10),rep("#FD8D3C",1),
                 rep("#5b3396",3),rep("#FD8D3C",1),
                 rep("#5b3396",1),rep("#FD8D3C",1),
                 rep("#5b3396",4),rep("#FD8D3C",1),
                 rep("#5b3396",2),rep("#FD8D3C",1),
                 rep("#5b3396",2),rep("#FD8D3C",2),
                 rep("#5b3396",2),rep("#FD8D3C",1),
                 rep("#5b3396",4))

ggplot(dat.sqc, aes(cohort, n_syn_singleton)) +
  geom_hline(yintercept=c(epi25_mean, control_mean), color=c("#5b3396","#FD8D3C"), linetype="dashed") +
  geom_boxplot(aes(fill=cohort),color='grey30',fill=cohort.color,
               outlier.size = 0.5,outlier.shape=1,
               alpha=0.6, width=0.6, show.legend=TRUE) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=1, color="white", fill="white") + 
  # geom_jitter(
  # width=0.15, size=0.7, aes(color=all_cohort), alpha=0.3, show.legend=FALSE) + 
  # ylim(0,60) +
  labs(title="Final sample QC",
       y="Synonymous singleton counts per individual",
       x="") +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_nSYNsingleton_by_cohort_ordered.pdf"), width=7, height=10)


# colored by case-control status rather than cohort
ggplot(dat.sqc, aes(cohort, n_syn_singleton)) +
  geom_boxplot(outlier.size=-1, coef=0, color='grey50', fill='grey95', show.legend=FALSE)+
  geom_jitter(width=0.2, size=1, aes(color=case_control), alpha=0.3, show.legend=TRUE) + 
  labs(title="Final sample QC",
       y="Number of synonymous singletons called per individual",
       x="Sample cohort") +
  scale_color_manual(values = col.cc) +
  coord_flip() +
  theme_bw()
ggsave(paste0("../plots/10_",pop,"_finalQC_sample_nSYNsingleton_by_cohort_ordered_colbyCC.pdf"), width=7, height=8)


# Check/Remove outliers
table(dat.sqc$cohort, dat.sqc$epi25_site=="")
# * also check/remove w/i cohort outliers
dat.sqc.nsyn.cohort <- dat.sqc %>% 
  group_by(cohort) %>%
  summarise(
    # n_syn_singleton
    avg_n_syn_singleton = mean(n_syn_singleton),
    sd_n_syn_singleton = sd(n_syn_singleton),
    upper_n_syn_singleton = avg_n_syn_singleton + 4*sd_n_syn_singleton,
    lower_n_syn_singleton = avg_n_syn_singleton - 4*sd_n_syn_singleton,
    nOutlier.n_syn_singleton = sum(n_syn_singleton > upper_n_syn_singleton | n_syn_singleton < lower_n_syn_singleton, na.rm=T)
  )
sum(dat.sqc.nsyn.cohort$nOutlier.n_syn_singleton)

outlier.n_syn_singleton = character()
for(cohort in dat.sqc.nsyn.cohort$cohort[dat.sqc.nsyn.cohort$nOutlier.n_syn_singleton >0]){
  print(cohort)
  # dat.sqc.tmp = dat.sqc[which(dat.sqc$epi25_cohort %in% cohort | dat.sqc$cohort %in% cohort),]
  dat.sqc.tmp = dat.sqc[which(dat.sqc$cohort %in% cohort),]
  outlier.index = which(
    dat.sqc.tmp$n_syn_singleton < dat.sqc.nsyn.cohort$lower_n_syn_singleton[dat.sqc.nsyn.cohort$cohort==cohort] |
      dat.sqc.tmp$n_syn_singleton > dat.sqc.nsyn.cohort$upper_n_syn_singleton[dat.sqc.nsyn.cohort$cohort==cohort]
  )
  outlier.n_syn_singleton = c(outlier.n_syn_singleton, dat.sqc.tmp$s[outlier.index])
  rm(dat.sqc.tmp)
}
table(dat.sqc[dat.sqc$s %in% outlier.n_syn_singleton, "case_control"])


####
# Write out sample QC lists
outlier.callRate <- dat.sqc$s[dat.sqc$call_rate < 0.98]
outlier_samples <- unique(c(outlier.rTiTv, outlier.rHetHom, outlier.rInDel,outlier.callRate))
write.table(outlier_samples, paste0("10_",pop,"_outlier_sample.remove.list"), col.names = F, quote=F, row.names = F)

# manually check/add nSYNsingleton outliers by pop
#######
# EUR #
#######
# Remove CYPCYP GPC_Latino_controls TURBZU LEBABM TURIBU (yr4)
# Remove CYPCYP TURBZU LEBABM TURIBU ITAUM3 Italian_controls
# nSYNsingleton_outlier_samples <- dat.sqc %>% filter(cohort=="CYPCYP" | cohort=="TURBZU" | cohort=="LEBABM" |
#                                                       cohort=="TURIBU" | cohort=="ITAUM3" | 
#                                                       cohort=="Italian_controls") %>% select(s)
# * further Remove PME_cohort ITAUBG ITAUMR ITAUMC ITAIGI ITAICB
nSYNsingleton_outlier_samples <- dat.sqc %>% filter(cohort=="CYPCYP" | cohort=="TURBZU" | cohort=="LEBABM" |
                                                      cohort=="TURIBU" | cohort=="ITAUM3" | 
                                                      cohort=="Italian_controls" | 
                                                      cohort=="PME_cohort" | cohort=="ITAUBG" | cohort=="ITAUMR" |  
                                                      cohort=="ITAUMC" | cohort=="ITAIGI" | cohort=="ITAICB")

#######
# AFR #
#######
nSYNsingleton_outlier_samples <- dat.sqc %>% filter(cohort=="DEUUKB" |  cohort=="AUSRMB" | 
                                                      cohort=="AUSAUS" | cohort=="ITAUMC" | cohort=="BELULB" | cohort=="FRALYU" | 
                                                      cohort=="BELATW")
#######
# EAS #
#######
nSYNsingleton_outlier_samples <- dat.sqc %>% filter(cohort=="NZLUTO" |  cohort=="USAMON" | 
                                                      cohort=="DEUUGS" | cohort=="ITAIGI" | cohort=="USAMSS")

#######
# SAS #
#######
nSYNsingleton_outlier_samples <- dat.sqc %>% filter(cohort=="FRALYU" |  cohort=="CHEUBB" |  cohort=="LEBABM" | 
                                                      cohort=="DEUUKL" | cohort=="UK/IRL Edinburgh + McQuillan" | 
                                                      cohort=="GBRSWU" |  cohort=="USAEGP" | cohort=="USAMGH" | cohort=="BELULB")

#######
# FIN #
#######
# only keep the three "FIN" cohorts
nSYNsingleton_outlier_samples <- dat.sqc %>% filter(cohort!="FINUVH" &  cohort!="FINKPH" &  cohort!="IBD_FINRISK" )

#######
# AMR #
#######
nSYNsingleton_outlier_samples <- dat.sqc %>% filter(cohort=="LEBABM" |  cohort=="ITAIGI" |  cohort=="ITAUMR" |  
                                                      cohort=="DEUUTB" | cohort=="BRAUSP" | cohort=="FRALYU" )


###
# write out per metirc
write.table(outlier.rTiTv, paste0("10_",pop,"_rTiTv_outlier_sample.list"), col.names = F, quote=F, row.names = F)
write.table(outlier.rHetHom, paste0("10_",pop,"_rHetHom_outlier_sample.list"), col.names = F, quote=F, row.names = F)
write.table(outlier.rInDel, paste0("10_",pop,"_rInDel_outlier_sample.list"), col.names = F, quote=F, row.names = F)
# write.table(outlier.n_syn_singleton, paste0("10_",pop,"_n_syn_singleton_outlier_sample.list"), col.names = F, quote=F, row.names = F)
nSYNsingleton_outlier_samples <- unique(c(nSYNsingleton_outlier_samples$s, outlier.n_syn_singleton))
table(dat.sqc[dat.sqc$s %in% nSYNsingleton_outlier_samples, "case_control"])
write.table(nSYNsingleton_outlier_samples, paste0("10_",pop,"_nSYNsingleton_outlier_sample.list"), col.names = F, quote=F, row.names = F)


