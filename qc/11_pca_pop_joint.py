import hail as hl
hl.init(log='/tmp/hail_pca_final.log')

from pprint import pprint
# Input
MT = 'gs://epi25/wes/genetic_data/epi25_cc_wes_hg38_alqc_split_gqc.mt'

PHENO_TABLE = 'gs://epi25/wes/annotations/samples_annot.ht'
SITES_TABLE = 'gs://epi25/wes/annotations/sites_annot.ht'
HIGH_LD_REGIONS = 'gs://epi25/wes/annotations/high-LD-regions-hg38-GRCh38.txt'

IMPUTED_SEX_TABLE = 'gs://epi25/wes/qc_files/04_imputed_sex.tsv'
MISSING_SEX_LIST = 'gs://epi25/wes/qc_files/04_missing_impSex_sample.remove.list'
DISCORDANT_SEX_LIST = 'gs://epi25/wes/qc_files/04_discordant_sex_sample.remove.list'

VARIANT_PREQC_LIST = 'gs://epi25/wes/qc_files/01_prefilter_variant.remove.list'
# combined lists across pops
JOINT_SAMPLE_LIST = 'gs://epi25/wes/qc_files/06_predicted_pops_sample.list'
JOINT_IBD_REMOVE_LIST = 'gs://epi25/wes/qc_files/08_ibd_pops_sample.remove.list'
JOINT_VARIANT_FINALQC_LIST = 'gs://epi25/wes/qc_files/09_pops_finalQC_variant.remove.list'
JOINT_SAMPLE_FINALQC_LIST = 'gs://epi25/wes/qc_files/10_pops_outlier_sample.remove.list'
JOINT_NSYNAC1_OUTLIER_LIST = 'gs://epi25/wes/qc_files/10_pops_nSYNsingleton_outlier_sample.list'
JOINT_POP_PCA_OUTLIER_LIST = 'gs://epi25/wes/qc_files/07_pops_pca_outlier_sample.list'

# Output
JOINT_PCA_SCORES_TABLE = 'gs://epi25/wes/qc_files/11_joint_pca.tsv'


############
print("")
print("Read in data")
mt = hl.read_matrix_table(MT)
ht_pop = hl.import_table(JOINT_SAMPLE_LIST, impute=True, no_header=True)
ht_pop = ht_pop.annotate(s = ht_pop.f0).key_by('s')

ht_site = hl.read_table(SITES_TABLE)
ht_pheno = hl.read_table(PHENO_TABLE)
ht_impsex = hl.import_table(IMPUTED_SEX_TABLE, impute=True, key="s")

ht_discordSex = hl.import_table(DISCORDANT_SEX_LIST, impute=True, no_header=True)
ht_discordSex = ht_discordSex.annotate(s = ht_discordSex.f0).key_by('s')
ht_missingSex = hl.import_table(MISSING_SEX_LIST, impute=True, no_header=True)
ht_missingSex = ht_missingSex.annotate(s = ht_missingSex.f0).key_by('s')
ht_ibd = hl.import_table(JOINT_IBD_REMOVE_LIST, impute=True, no_header=True)
ht_ibd = ht_ibd.annotate(s = ht_ibd.f0).key_by('s')

ht_vqc1 = (hl.import_table(VARIANT_PREQC_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')
ht_vqc2 = (hl.import_table(JOINT_VARIANT_FINALQC_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')

ht_sqc1 = hl.import_table(JOINT_SAMPLE_FINALQC_LIST, impute=True, no_header=True)
ht_sqc1 = ht_sqc1.annotate(s = ht_sqc1.f0).key_by('s')
ht_sqc2 = hl.import_table(JOINT_NSYNAC1_OUTLIER_LIST, impute=True, no_header=True)
ht_sqc2 = ht_sqc2.annotate(s = ht_sqc2.f0).key_by('s')
ht_sqc3 = hl.import_table(JOINT_POP_PCA_OUTLIER_LIST, impute=True, no_header=True)
ht_sqc3 = ht_sqc3.annotate(s = ht_sqc3.f0).key_by('s')

ht_highLD = hl.import_bed(HIGH_LD_REGIONS, reference_genome='GRCh38')


print("")
print("Annotate with phenotype info and imputed sex")
mt = mt.annotate_cols(impsex = ht_impsex[mt.s], pheno_tmp=ht_pheno[mt.s])
mt = mt.annotate_cols(pheno = mt.pheno_tmp.pheno)
mt = mt.drop('pheno_tmp')

print("Annotate with variant annotations")
mt = mt.annotate_rows(site = ht_site[mt.row_key])
mt.describe()


print("")
print("Retain non-Finnish, unrelated, sex checked JOINT samples")
mt = mt.filter_cols(hl.is_defined(ht_pop[mt.s]))
mt = mt.filter_cols(~hl.is_defined(ht_fin[mt.s]))
mt = mt.filter_cols(hl.is_defined(ht_discordSex[mt.s]), keep=False)
mt = mt.filter_cols(hl.is_defined(ht_missingSex[mt.s]), keep=False)
mt = mt.filter_cols(hl.is_defined(ht_ibd[mt.s]), keep=False)


print("")
print("Remove samples failing final QC")
mt = mt.filter_cols(~hl.is_defined(ht_sqc1[mt.s]))
mt = mt.filter_cols(~hl.is_defined(ht_sqc2[mt.s]))
mt = mt.filter_cols(~hl.is_defined(ht_sqc3[mt.s]))


print("")
print("Remove variants that fail initial and final QC")
mt = mt.filter_rows(hl.is_defined(ht_vqc1[mt.locus, mt.alleles]), keep=False)
mt = mt.filter_rows(hl.is_defined(ht_vqc2[mt.locus, mt.alleles]), keep=False)


print("")
print("DROP EPI4k cohort!")
mt = mt.filter_cols(mt.pheno.cohort=="Epi4k", keep=False)


# LD pruning
print("")
print("Remove variants in high LD region")
mt = mt.filter_rows(hl.is_defined(ht_highLD[mt.locus]), keep=False)
print("") # from Duncan
print("Keep in_x_nonpar or in_autosome_or_par")
mt = mt.filter_rows(mt.locus.in_x_nonpar() | mt.locus.in_autosome_or_par())
print("")
print("Perform variant QC")
mt = hl.variant_qc(mt)
print("")
print("Keep common (MAF>1%) variants for pruning")
mt = mt.filter_rows(mt.locus.in_autosome() & (mt.variant_qc.AF[1] >= 0.01) & (mt.variant_qc.AF[1] <= 0.99) & (mt.variant_qc.call_rate >= 0.98)).persist()
print("")
print("Perform LD pruning")
ht_pruned = hl.ld_prune(mt.GT, r2=0.005, bp_window_size=10000000)
mt_pruned = mt.filter_rows(hl.is_defined(ht_pruned[mt.row_key]))


print("")
print("Perform PCA")
eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt_pruned.GT, k=20, compute_loadings=True)


print("")
print("Export JOINT sample scores and loadings")
scores.export(JOINT_PCA_SCORES_TABLE)

