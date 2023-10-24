import hail as hl
hl.init(log='/tmp/epi25_hail_final_sqc.log')

from pprint import pprint
# Input
MT = 'gs://epi25/wes/genetic_data/epi25_cc_wes_hg38_alqc_split_gqc.mt'

PHENO_TABLE = 'gs://epi25/wes/annotations/samples_annot.ht'
SITES_TABLE = 'gs://epi25/wes/annotations/sites_annot.ht'

IMPUTED_SEX_TABLE = 'gs://epi25/wes/qc_files/04_imputed_sex.tsv'
MISSING_SEX_LIST = 'gs://epi25/wes/qc_files/04_missing_impSex_sample.remove.list'
DISCORDANT_SEX_LIST = 'gs://epi25/wes/qc_files/04_discordant_sex_sample.remove.list'

POP = 'EUR'
POP_SAMPLE_LIST = 'gs://epi25/wes/qc_files/06_predicted_{0}_sample.list'.format(POP)
POP_IBD_REMOVE_LIST = 'gs://epi25/wes/qc_files/08_ibd_{0}_sample.remove.list'.format(POP.lower())
FIN_SAMPLE_LIST = 'gs://epi25/wes/qc_files/07_predicted_FIN_sample.list'

VARIANT_PREQC_LIST = 'gs://epi25/wes/qc_files/01_prefilter_variant.remove.list'
POP_VARIANT_FINALQC_LIST = 'gs://epi25/wes/qc_files/09_{0}_finalQC_variant.remove.list'.format(POP.lower())


# Output
POP_SAMPLE_FINALQC_TABLE = 'gs://epi25/wes/qc_files/10_{0}_finalQC_sample.tsv'.format(POP.lower())



print("")
print("Read in data")
mt = hl.read_matrix_table(MT)
ht_pop = hl.import_table(POP_SAMPLE_LIST, impute=True, no_header=True)
ht_pop = ht_pop.annotate(s = ht_pop.f0).key_by('s')
ht_fin = hl.import_table(FIN_SAMPLE_LIST, impute=True, no_header=True)
ht_fin = ht_fin.annotate(s = ht_fin.f0).key_by('s')

ht_site = hl.read_table(SITES_TABLE)
ht_pheno = hl.read_table(PHENO_TABLE)
ht_impsex = hl.import_table(IMPUTED_SEX_TABLE, impute=True, key="s")

ht_discordSex = hl.import_table(DISCORDANT_SEX_LIST, impute=True, no_header=True)
ht_discordSex = ht_discordSex.annotate(s = ht_discordSex.f0).key_by('s')
ht_missingSex = hl.import_table(MISSING_SEX_LIST, impute=True, no_header=True)
ht_missingSex = ht_missingSex.annotate(s = ht_missingSex.f0).key_by('s')
ht_ibd = hl.import_table(POP_IBD_REMOVE_LIST, impute=True, no_header=True)
ht_ibd = ht_ibd.annotate(s = ht_ibd.f0).key_by('s')

ht_vqc1 = (hl.import_table(VARIANT_PREQC_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')
ht_vqc2 = (hl.import_table(POP_VARIANT_FINALQC_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')


print("Annotate with phenotype info and imputed sex")
mt = mt.annotate_cols(impsex = ht_impsex[mt.s], pheno_tmp=ht_pheno[mt.s])
mt = mt.annotate_cols(pheno = mt.pheno_tmp.pheno)
mt = mt.drop('pheno_tmp')

print("Annotate with variant annotations")
mt = mt.annotate_rows(site = ht_site[mt.row_key])
mt.describe()

print("")
print("Retain unrelated, sex checked POP samples")
mt = mt.filter_cols(hl.is_defined(ht_pop[mt.s]))
if POP == 'EUR': mt = mt.filter_cols(~hl.is_defined(ht_fin[mt.s]))
mt = mt.filter_cols(hl.is_defined(ht_discordSex[mt.s]), keep=False)
mt = mt.filter_cols(hl.is_defined(ht_missingSex[mt.s]), keep=False)
mt = mt.filter_cols(hl.is_defined(ht_ibd[mt.s]), keep=False)


print("")
print("Remove variants that fail initial and final QC")
mt = mt.filter_rows(hl.is_defined(ht_vqc1[mt.locus, mt.alleles]), keep=False)
mt = mt.filter_rows(hl.is_defined(ht_vqc2[mt.locus, mt.alleles]), keep=False)


print("")
print("Perform sample QC")
mt = hl.sample_qc(mt)


print("")
print("Perform variant QC and add an indicator for singleton")
mt = hl.variant_qc(mt)
mt = mt.annotate_rows(isSingleton = (mt.variant_qc.n_non_ref==1))


print("")
print("Add annotation-specific singleton count")
mt = mt.annotate_cols(
    n_nonref_singleton = (hl.agg.filter(mt.isSingleton, hl.agg.count_where(mt.GT.is_non_ref()))),
    n_syn_singleton = (hl.agg.filter(mt.isSingleton & (mt.site.consequence_category == "synonymous"), hl.agg.count_where(mt.GT.is_non_ref()))),
    n_ptv_singleton = (hl.agg.filter(mt.isSingleton & (mt.site.consequence_category == "pLoF"), hl.agg.count_where(mt.GT.is_non_ref()))),
    n_mis1_singleton = (hl.agg.filter(mt.isSingleton & (mt.site.consequence_category == "other_missense"), hl.agg.count_where(mt.GT.is_non_ref()))),
    n_mis3_singleton = (hl.agg.filter(mt.isSingleton & (mt.site.consequence_category == "damaging_missense"), hl.agg.count_where(mt.GT.is_non_ref()))),
    n_noncoding_singleton = (hl.agg.filter(mt.isSingleton & (mt.site.consequence_category == "non_coding"), hl.agg.count_where(mt.GT.is_non_ref())))
    )


print("")
print("Output sample QC metrics")
ht_sqc = mt.cols()
ht_sqc = ht_sqc.select(batch = ht_sqc.pheno.batch, cohort = ht_sqc.pheno.cohort,
                         epi25_site = ht_sqc.pheno.epi25_site,
                         case_control = ht_sqc.pheno.case_control,
                         epilepsy_type = ht_sqc.pheno.epilepsy_type,
                         capture_set = ht_sqc.pheno.bait_set,
                         contamination = ht_sqc.pheno.pct_contamination,
                         pct_chimeras = ht_sqc.pheno.pct_chimeras,
                         call_rate = ht_sqc.sample_qc.call_rate, 
                         n_het = ht_sqc.sample_qc.n_het,
                         dpMean = ht_sqc.sample_qc.dp_stats.mean, 
                         gqMean = ht_sqc.sample_qc.gq_stats.mean,
                         rTiTv = ht_sqc.sample_qc.r_ti_tv,
                         rHetHom = ht_sqc.sample_qc.r_het_hom_var,
                         rIndel = ht_sqc.sample_qc.r_insertion_deletion,
                         n_singleton = ht_sqc.sample_qc.n_singleton,
                         n_nonref_singleton = ht_sqc.n_nonref_singleton,
                         n_syn_singleton = ht_sqc.n_syn_singleton,
                         n_ptv_singleton = ht_sqc.n_ptv_singleton,
                         n_mis1_singleton = ht_sqc.n_mis1_singleton,
                         n_mis3_singleton = ht_sqc.n_mis3_singleton,
                         n_noncoding_singleton = ht_sqc.n_noncoding_singleton
                         )
ht_sqc.export(POP_SAMPLE_FINALQC_TABLE)

