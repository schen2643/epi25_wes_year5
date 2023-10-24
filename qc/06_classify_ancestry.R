setwd("~/epi25/wes/qc_files")

library(tidyverse)
library(data.table)

pca.1kg <- fread("05_pca_w_1kg.tsv", header=T, sep="\t")
epi25_pheno <- fread("../annotations/epi25_pheno.tsv", header=T, sep="\t", stringsAsFactors = F)
onekG_pheno <- fread("../misc_data/samples_1kg_30x_GRCh38.tsv", header=T, sep="\t")
pca.1kg <- merge(pca.1kg, epi25_pheno, by="sample", all.x=T)
pca.1kg <- merge(pca.1kg, onekG_pheno[,c("sample","pop","super_pop")], all.x=T)
pca.1kg <- pca.1kg[pca.1kg$sample %in% c(epi25_pheno$sample, onekG_pheno$sample)]

pca.1kg$super_pop <- ifelse(is.na(pca.1kg$super_pop), pca.1kg$case_control, pca.1kg$super_pop)
pca.1kg$pop <- ifelse(is.na(pca.1kg$pop), pca.1kg$case_control, pca.1kg$pop)
pca.1kg$case_control <- ifelse(is.na(pca.1kg$case_control), "1KG", pca.1kg$case_control)
pca.1kg$cohort <- ifelse(is.na(pca.1kg$cohort), "1KG", pca.1kg$cohort)
pca.1kg$epilepsy_type <- ifelse(is.na(pca.1kg$epilepsy_type), "1KG", pca.1kg$epilepsy_type)

cohort_names <- names(table(pca.1kg$cohort))
epi25_site_names <- names(table(pca.1kg$epi25_site))[-1]
cohort_names[!(cohort_names %in% epi25_site_names)]
pca.1kg$cohort <- factor(pca.1kg$cohort, levels=c("1KG", names(table(pca.1kg$epi25_site))[-1],
                                                  cohort_names[!(cohort_names %in% epi25_site_names)][-1]))

pca.1kg$super_pop = factor(pca.1kg$super_pop,
                           levels=c("control","case","AFR","AMR","EAS","EUR","SAS"))
pca.1kg$super_pop <- recode(pca.1kg$super_pop, `case` = "Epilepsy", `control` = "Control")
pca.1kg$case_control <- recode(pca.1kg$case_control, `case` = "Epilepsy", `control` = "Control")

# RF function to predict ancestry using PCs:
pop_forest <- function(training_data, data, ntree=100, seed=42, pcs=1:6) {
  set.seed(seed)
  form = formula(paste('as.factor(known_pop) ~', paste0('PC', pcs, collapse = ' + ')))
  forest = randomForest(form,
                        data = training_data,
                        importance = T,
                        ntree = ntree)
  print(forest)
  fit_data = data.frame(predict(forest, data, type='prob'), sample = data$sample)
  fit_data %>%
    gather(predicted_pop, probability, -sample) %>%
    group_by(sample) %>%
    slice(which.max(probability))
}

# Prep data
# Explanation:
# `df.train` and `df.test` are training and testing data. training data has to have a 
# column `known_pop` and `PC1` to `PC10` or so. `df.test` is expected to have a 
# column `sample` which is just the sample ID, and also PC columns.
df.train <- pca.1kg %>% 
  filter(case_control=='1KG') %>%
  select(super_pop, PC1:PC10) %>%
  rename(known_pop = super_pop)
df.train$known_pop <- as.character(df.train$known_pop)
df.test <- pca.1kg %>% 
  filter(case_control!='1KG') %>%
  select(sample, PC1:PC10)
df.test$sample = as.character(df.test$sample)

# Make prediction
df.pred <- as.data.frame(pop_forest(training_data = df.train, data = df.test))
df.pred.90 <- df.pred %>% filter(probability >= 0.9)
pca.1kg.pred <- merge(pca.1kg, df.pred, by="sample")
pca.1kg.pred.90 <- subset(pca.1kg.pred, probability >= 0.9)
# Predict Finns
df.train <- filter(pca.1kg, case_control=='1KG') %>%
  mutate(Finnish = pop == 'FIN') %>%
  select(c(Finnish, PC1, PC2))
df.test <- filter(pca.1kg, case_control!='1KG') %>%
  select(c(sample, PC1, PC2))
nb <- naiveBayes(Finnish~., data=df.train)

pca.1kg.pred.fin <- mutate(df.test, predicted_fin = predict(nb, df.test[,-1])) %>%
  mutate(predicted_fin = as.factor(ifelse(predicted_fin, 'Finnish', 'non-Finnish'))) %>%
  select(sample, predicted_fin) %>%
  inner_join(pca.1kg, by='sample')


# Write out predicted ancestries
predAFR_samples <- pca.1kg.pred.90 %>% filter(predicted_pop=="AFR") %>% select(sample)
predAMR_samples <- pca.1kg.pred.90 %>% filter(predicted_pop=="AMR") %>% select(sample)
predEAS_samples <- pca.1kg.pred.90 %>% filter(predicted_pop=="EAS") %>% select(sample)
predEUR_samples <- pca.1kg.pred.90 %>% filter(predicted_pop=="EUR") %>% select(sample)
predSAS_samples <- pca.1kg.pred.90 %>% filter(predicted_pop=="SAS") %>% select(sample)
predFIN_samples <- pca.1kg.pred.fin %>% filter(predicted_fin=="Finnish") %>% select(sample)

write.table(predAFR_samples, "06_predicted_AFR_sample.list", col.names = F, quote=F, row.names = F)
write.table(predAMR_samples, "06_predicted_AMR_sample.list", col.names = F, quote=F, row.names = F)
write.table(predEAS_samples, "06_predicted_EAS_sample.list", col.names = F, quote=F, row.names = F)
write.table(predEUR_samples, "06_predicted_EUR_sample.list", col.names = F, quote=F, row.names = F)
write.table(predSAS_samples, "06_predicted_SAS_sample.list", col.names = F, quote=F, row.names = F)
write.table(predFIN_samples, "06_predicted_FIN_sample.list", col.names = F, quote=F, row.names = F)

write.table(pca.1kg.pred,"06_pca_w_1kg_RFpredicted_pop.tsv", col.names=T, row.names=F, quote=F, 
            sep="\t")

