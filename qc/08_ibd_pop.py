import hail as hl
hl.init()

from pprint import pprint

# Input
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc_hardcalls.mt'
POP = "EUR" # AFR, EAS, SAS, AMR, FIN
POP_SAMPLE_LIST = 'gs://epi25/wes/qc_files/06_predicted_{0}_sample.list'.format(POP)
PRUNED_AUTO_LIST = 'gs://epi25/wes/qc_files/03_pruned_autosomal_variant.keep.list'

# Output
POP_IBD_EST_TABLE = 'gs://epi25/wes/qc_files/08_ibd_{0}.tsv'.format(POP.lower())


print("")
print("Read in data")
mt = hl.read_matrix_table(MT)
ht_pop = hl.import_table(POP_SAMPLE_LIST, impute=True, no_header=True)
ht_pop = ht_pop.annotate(s = ht_pop.f0).key_by('s')
ht_pruned_auto = (hl.import_table(PRUNED_AUTO_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')


print("")
print("Retain single-pop samples and pruned variants")
mt_pop = mt.filter_cols(hl.is_defined(ht_pop[mt.s]))
mt_pop = mt_pop.filter_rows(hl.is_defined(ht_pruned_auto[mt_pop.row_key]))


print("")
print("Estimate IBD")
ibd_pop = hl.identity_by_descent(mt_pop, min=0.125)


print("")
print("Write out IBD metrics")
# columns i and j (keys here) are automatically included, so no need to export again
ibd_pop = ibd_pop.select(sample_i = ibd_pop.i, sample_j = ibd_pop.j, ibd_Z0 = ibd_pop.ibd.Z0,
                         ibd_Z1 = ibd_pop.ibd.Z1, ibd_Z2 = ibd_pop.ibd.Z2, PI_HAT = ibd_pop.ibd.PI_HAT,
                         ibs0 = ibd_pop.ibs0, ibs1 = ibd_pop.ibs1, ibs2 = ibd_pop.ibs2)
ibd_pop.export(POP_IBD_EST_TABLE)
