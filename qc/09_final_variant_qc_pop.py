import hail as hl
hl.init(log='/tmp/epi25_final_vqc.log')

from pprint import pprint
# Input
MT = 'gs://epi25/wes/genetic_data/epi25_cc_wes_hg38_alqc_split_gqc.mt'

POP = 'EUR'
POP_SAMPLE_LIST = 'gs://epi25/wes/qc_files/06_predicted_{0}_sample.list'.format(POP)

PHENO_TABLE = 'gs://epi25/wes/annotations/samples_annot.ht'
IMPUTED_SEX_TABLE = 'gs://epi25/wes/qc_files/04_imputed_sex.tsv'
MISSING_SEX_LIST = 'gs://epi25/wes/qc_files/04_missing_impSex_sample.remove.list'
DISCORDANT_SEX_LIST = 'gs://epi25/wes/qc_files/04_discordant_sex_sample.remove.list'

POP_IBD_REMOVE_LIST = 'gs://epi25/wes/qc_files/08_ibd_{0}_sample.remove.list'.format(POP.lower())
FIN_SAMPLE_LIST = 'gs://epi25/wes/qc_files/07_predicted_FIN_sample.list'
VARIANT_PREQC_LIST = 'gs://epi25/wes/qc_files/01_prefilter_variant.remove.list'

# Output
POP_VARIANT_FINALQC_TABLE = 'gs://epi25/wes/qc_files/09_{0}_finalQC_variant.tsv'.format(POP.lower())
POP_VARIANT_FINALQC_LIST = 'gs://epi25/wes/qc_files/09_{0}_finalQC_variant.remove.list'.format(POP.lower())



print("")
print("Read in data")
mt = hl.read_matrix_table(MT)
ht_pop = hl.import_table(POP_SAMPLE_LIST, impute=True, no_header=True)
ht_pop = ht_pop.annotate(s = ht_pop.f0).key_by('s')
ht_fin = hl.import_table(FIN_SAMPLE_LIST, impute=True, no_header=True)
ht_fin = ht_fin.annotate(s = ht_fin.f0).key_by('s')

ht_pheno = hl.read_table(PHENO_TABLE)
ht_impsex = hl.import_table(IMPUTED_SEX_TABLE, impute=True, key="s")
ht_discordSex = hl.import_table(DISCORDANT_SEX_LIST, impute=True, no_header=True)
ht_discordSex = ht_discordSex.annotate(s = ht_discordSex.f0).key_by('s')
ht_missingSex = hl.import_table(MISSING_SEX_LIST, impute=True, no_header=True)
ht_missingSex = ht_missingSex.annotate(s = ht_missingSex.f0).key_by('s')
ht_ibd = hl.import_table(POP_IBD_REMOVE_LIST, impute=True, no_header=True)
ht_ibd = ht_ibd.annotate(s = ht_ibd.f0).key_by('s')

ht_vqc = (hl.import_table(VARIANT_PREQC_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')


print("")
print("Annotate with phenotype info and imputed sex")
mt = mt.annotate_cols(impsex = ht_impsex[mt.s], pheno_tmp=ht_pheno[mt.s])
mt = mt.annotate_cols(pheno = mt.pheno_tmp.pheno)
mt = mt.drop('pheno_tmp')


print("")
print("Retain unrelated, sex checked POP samples")
mt = mt.filter_cols(hl.is_defined(ht_pop[mt.s]))
if POP == "EUR": mt = mt.filter_cols(~hl.is_defined(ht_fin[mt.s]))
mt = mt.filter_cols(hl.is_defined(ht_discordSex[mt.s]), keep=False)
mt = mt.filter_cols(hl.is_defined(ht_missingSex[mt.s]), keep=False)
mt = mt.filter_cols(hl.is_defined(ht_ibd[mt.s]), keep=False)


print("")
print("Remove variants that fail pre-filter")
mt = mt.filter_rows(hl.is_defined(ht_vqc[mt.locus, mt.alleles]), keep=False)


print("")
print("Perform variant QC")
mt = hl.variant_qc(mt)

mt = mt.annotate_rows(variant_qc = mt.variant_qc.annotate(p_value_hwe = hl.case()
    .when(mt.locus.in_autosome(), mt.variant_qc.p_value_hwe)
    .default(hl.agg.filter(mt.impsex.is_female,
        hl.agg.hardy_weinberg_test(mt.GT)['p_value'])))
)


print("")
print("Perform variant QC among the cases and the controls separately")
mt_case = mt.filter_cols(mt.pheno.case_control == "case")
mt_ctrl = mt.filter_cols(mt.pheno.case_control == "control")
mt_case = hl.variant_qc(mt_case)
mt_ctrl = hl.variant_qc(mt_ctrl)


print("")
print("Annotate mt with case and control-specific QC metrics")
ht_case_vqc = mt_case.rows()
ht_ctrl_vqc = mt_ctrl.rows()

mt = mt.annotate_rows(case_variant_qc = ht_case_vqc[mt.locus, mt.alleles].variant_qc)
mt = mt.annotate_rows(ctrl_variant_qc = ht_ctrl_vqc[mt.locus, mt.alleles].variant_qc)


print("")
print('Output final variant QC metrics as a Table')
ht_vqc = mt.rows()
ht_vqc = ht_vqc.select(
    AC = ht_vqc.variant_qc.AC[1],
    AF = ht_vqc.variant_qc.AF[1],
    homozygote_count=ht_vqc.variant_qc.homozygote_count[1],
    call_rate = ht_vqc.variant_qc.call_rate,
    case_call_rate = ht_vqc.case_variant_qc.call_rate,
    ctrl_call_rate = ht_vqc.ctrl_variant_qc.call_rate,
    n_filtered = ht_vqc.variant_qc.n_filtered,
    n_het = ht_vqc.variant_qc.n_het,
    dpMean = ht_vqc.variant_qc.dp_stats.mean, 
    gqMean = ht_vqc.variant_qc.gq_stats.mean,
    pHWE = ht_vqc.variant_qc.p_value_hwe
    )
ht_vqc.export(POP_VARIANT_FINALQC_TABLE)


print("")
print("Filter to variants that fail these QC metrics for removal")
mt_rm = mt.filter_rows(((mt.variant_qc.AC[0] == 0) | (mt.variant_qc.AC[1] == 0)) | (mt.variant_qc.p_value_hwe < 1e-06) 
	| (mt.variant_qc.call_rate < 0.98)
	| (mt.case_variant_qc.call_rate < 0.98)
	| (mt.ctrl_variant_qc.call_rate < 0.98)
	| ((mt.case_variant_qc.call_rate - mt.ctrl_variant_qc.call_rate) > 0.02))


print("")
pprint("Write out the list of variants to remove")
mt_rm.rows().select().export(POP_VARIANT_FINALQC_LIST)

