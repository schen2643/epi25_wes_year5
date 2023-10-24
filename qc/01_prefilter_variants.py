import hail as hl
hl.init()

from pprint import pprint

# Input
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc_hardcalls.mt'
LCR_INTVL = 'gs://epi25/wes/annotations/LCRFromHengHg38.bed'
ICE_TARGET_INTVL = 'gs://epi25/wes/annotations/ice_coding_v1_targets.interval_list'
ICE_PADDED_INTVL = 'gs://epi25/wes/annotations/ice_coding_v1_padded_targets.interval_list'
ICE_PADDED_100_INTVL = 'gs://epi25/wes/annotations/ice_coding_v1_padded_targets_100.interval_list'
ICE_PADDED_150_INTVL = 'gs://epi25/wes/annotations/ice_coding_v1_padded_targets_150.interval_list'

TWIST_TARGET_INTVL = 'gs://epi25/wes/annotations/twist_jan18_grch38_targets.interval_list'
TWIST_PADDED_INTVL = 'gs://epi25/wes/annotations/twist_jan18_grch38_padded_targets.interval_list'
TWIST_PADDED_100_INTVL = 'gs://epi25/wes/annotations/twist_jan18_grch38_padded_targets_100.interval_list'
TWIST_PADDED_150_INTVL = 'gs://epi25/wes/annotations/twist_jan18_grch38_padded_targets_150.interval_list'


# Output
VARIANT_PREQC_LIST = 'gs://epi25/wes/qc_files/01_prefilter_variant.remove.list'



print('Read in data')
# Read in data
mt = hl.read_matrix_table(MT)

lcr = hl.import_bed(LCR_INTVL, reference_genome='GRCh38')
ice_target = hl.import_locus_intervals(ICE_TARGET_INTVL, reference_genome='GRCh38')
ice_padded = hl.import_locus_intervals(ICE_PADDED_INTVL, reference_genome='GRCh38')
ice_padded_100 = hl.import_locus_intervals(ICE_PADDED_100_INTVL, reference_genome='GRCh38')
ice_padded_150 = hl.import_locus_intervals(ICE_PADDED_150_INTVL, reference_genome='GRCh38')

twist_target = hl.import_locus_intervals(TWIST_TARGET_INTVL, reference_genome='GRCh38')
twist_padded = hl.import_locus_intervals(TWIST_PADDED_INTVL, reference_genome='GRCh38')
twist_padded_100 = hl.import_locus_intervals(TWIST_PADDED_100_INTVL, reference_genome='GRCh38')
twist_padded_150 = hl.import_locus_intervals(TWIST_PADDED_150_INTVL, reference_genome='GRCh38')



print('Annotate with VQSR pass/fail')
mt = mt.annotate_rows(failVQSR = (hl.len(mt.filters) != 0))

print('Annotate with LCR interval')
mt = mt.annotate_rows(inLCR = hl.is_defined(lcr[mt.locus]))

print('Annotate sites with platform')
mt = mt.annotate_rows(inICE_target = hl.is_defined(ice_target[mt.locus]),
                      inICE_padded = hl.is_defined(ice_padded[mt.locus]),
                      inICE_padded_100 = hl.is_defined(ice_padded_100[mt.locus]),
                      inICE_padded_150 = hl.is_defined(ice_padded_150[mt.locus]),
                      inTWIST_target = hl.is_defined(twist_target[mt.locus]),
                      inTWIST_padded = hl.is_defined(twist_padded[mt.locus]),
                      inTWIST_padded_100 = hl.is_defined(twist_padded_100[mt.locus]),
                      inTWIST_padded_150 = hl.is_defined(twist_padded_150[mt.locus]))



print('Annotate with variant AC')
mt = hl.variant_qc(mt, name='qc')


# Filter to variants that fail these QC metrics for removal
mt_remove = mt.filter_rows(mt.failVQSR | mt.inLCR | ~(mt.inICE_target & mt.inTWIST_target) | ((mt.qc.AC[0] == 0.0) | (mt.qc.AC[1] == 0.0)))



pprint('Write out the list of variants to remove')
mt_remove.rows().select().export(VARIANT_PREQC_LIST)

