import hail as hl
hl.init(log='/tmp/hail_write_qced_mt.log')

from pprint import pprint
# Input
MT = 'gs://epi25/wes/genetic_data/epi25_cc_wes_hg38_alqc_split_gqc.mt'

PHENO_TABLE = 'gs://epi25/wes/annotations/samples_annot.ht'
SITES_TABLE = 'gs://epi25/wes/annotations/sites_annot.ht'
HIGH_LD_REGIONS = 'gs://epi25/wes/annotations/high-LD-regions-hg38-GRCh38.txt'

IMPUTED_SEX_TABLE = 'gs://epi25/wes/qc_files/04_imputed_sex.tsv'
MISSING_SEX_LIST = 'gs://epi25/wes/qc_files/04_missing_impSex_sample.remove.list'
DISCORDANT_SEX_LIST = 'gs://epi25/wes/qc_files/04_discordant_sex_sample.remove.list'

POP = 'EUR'
POP_SAMPLE_LIST = 'gs://epi25/wes/qc_files/06_predicted_{0}_sample.list'.format(POP)
POP_IBD_REMOVE_LIST = 'gs://epi25/wes/qc_files/08_ibd_{0}_sample.remove.list'.format(POP.lower())
FIN_SAMPLE_LIST = 'gs://epi25/wes/qc_files/07_predicted_FIN_sample.list'

VARIANT_PREQC_LIST = 'gs://epi25/wes/qc_files/01_prefilter_variant.remove.list'
POP_VARIANT_FINALQC_LIST = 'gs://epi25/wes/qc_files/09_{0}_finalQC_variant.remove.list'.format(POP.lower())
POP_SAMPLE_FINALQC_LIST = 'gs://epi25/wes/qc_files/10_{0}_outlier_sample.remove.list'.format(POP.lower())
POP_NSYNAC1_OUTLIER_LIST = 'gs://epi25/wes/qc_files/10_{0}_nSYNsingleton_outlier_sample.list'.format(POP.lower())


# Output
MT_QCed = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_split_annot_qc_{0}.mt'.format(POP.lower())
SAMPLE_LIST_QCed = 'gs://epi25/wes/qc_files/12_{0}_finalQC_sample.keep.list'.format(POP.lower())



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

ht_sqc1 = hl.import_table(POP_SAMPLE_FINALQC_LIST, impute=True, no_header=True)
ht_sqc1 = ht_sqc1.annotate(s = ht_sqc1.f0).key_by('s')
ht_sqc2 = hl.import_table(POP_NSYNAC1_OUTLIER_LIST, impute=True, no_header=True)
ht_sqc2 = ht_sqc2.annotate(s = ht_sqc2.f0).key_by('s')


print("")
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
print("Remove samples failing final QC")
mt = mt.filter_cols(~hl.is_defined(ht_sqc1[mt.s]))
mt = mt.filter_cols(~hl.is_defined(ht_sqc2[mt.s]))


print("")
print("Remove variants that fail initial and final QC")
mt = mt.filter_rows(hl.is_defined(ht_vqc1[mt.locus, mt.alleles]), keep=False)
mt = mt.filter_rows(hl.is_defined(ht_vqc2[mt.locus, mt.alleles]), keep=False)


print("")
print("DROP EPI4k cohort!")
mt = mt.filter_cols(mt.pheno.cohort=="Epi4k", keep=False)


print("Perform sample and variant QC on the final data set")
mt = hl.sample_qc(mt)
mt = hl.variant_qc(mt)


print("Remove invariants/AC=0")
mt = mt.filter_rows((mt.variant_qc.AC[0] == 0) | (mt.variant_qc.AC[1] == 0), keep=False) 


print("Add annotations for later use")
# variant annotations
mt = mt.annotate_rows(
    ptv = mt.site.consequence_category == "pLoF",
    mis1 = mt.site.consequence_category == "other_missense",
    mis3 = mt.site.consequence_category == "damaging_missense",
    mis = (mt.site.consequence_category == "damaging_missense") | (mt.site.consequence_category == "other_missense"),
    syn = mt.site.consequence_category == "synonymous",
    mpc2 = mt.site.MPC >=2,
    lofteeHC = mt.site.vep.worst_csq_for_variant_canonical.lof == "HC",
    isSingleton = (mt.variant_qc.n_non_ref==1),
    nonRefSamples = hl.agg.filter((mt.variant_qc.AF[1] < 0.0005) & (mt.GT.is_non_ref()), hl.agg.collect(mt.s)),
    nonRefCohorts = hl.agg.filter((mt.variant_qc.AF[1] < 0.0005) & (mt.GT.is_non_ref()), hl.agg.collect(mt.pheno.cohort)),
    nonRefEpitype = hl.agg.filter((mt.variant_qc.AF[1] < 0.0005) & (mt.GT.is_non_ref()), hl.agg.collect(mt.pheno.epilepsy_type)),
    nonRefHetSamples = hl.agg.filter((mt.variant_qc.AF[1] < 0.0005) & (mt.GT.is_non_ref() & mt.GT.is_het() ), hl.agg.collect(mt.s)),
    nonRefHomSamples = hl.agg.filter((mt.variant_qc.AF[1] < 0.0005) & (mt.GT.is_non_ref() & mt.GT.is_hom_var() ), hl.agg.collect(mt.s)),
    )

print("Write out the qc'ed MT")
mt.write(MT_QCed, overwrite=True)

print("Write out the list of qc'ed samples to keep")
mt.cols().select().export(SAMPLE_LIST_QCed)

