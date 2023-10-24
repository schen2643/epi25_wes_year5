import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Input
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc_hardcalls.mt'
VARIANT_PREQC_LIST = 'gs://epi25/wes/qc_files/01_prefilter_variant.remove.list'
SAMPLE_INITQC_LIST = 'gs://epi25/wes/qc_files/02_initQC_sample.remove.list'
HIGH_LD_REGIONS = 'gs://epi25/wes/annotations/high-LD-regions-hg38-GRCh38.txt'

# Output
PRUNED_chrX_LIST = 'gs://epi25/wes/qc_files/03_pruned_chrX_variant.keep.list'
PRUNED_AUTO_LIST = 'gs://epi25/wes/qc_files/03_pruned_autosomal_variant.keep.list'
#########


print("Read in data")
mt = hl.read_matrix_table(MT)
ht_vqc = (hl.import_table(VARIANT_PREQC_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')
ht_sqc = (hl.import_table(SAMPLE_INITQC_LIST, impute=True).key_by('s'))
ht_highLD = hl.import_bed(HIGH_LD_REGIONS, reference_genome='GRCh38')



print("")
print("Remove variants and samples that failed init QC")
mt = mt.filter_rows(hl.is_defined(ht_vqc[mt.locus, mt.alleles]), keep=False)
mt = mt.filter_cols(hl.is_defined(ht_sqc[mt.s]), keep=False)


print("")
print("Remove variants in high LD region")
mt = mt.filter_rows(hl.is_defined(ht_highLD[mt.locus]), keep=False)


print("") 
print("Keep in_x_nonpar or in_autosome_or_par")
mt = mt.filter_rows(mt.locus.in_x_nonpar() | mt.locus.in_autosome_or_par())


print("")
print("Perform variant QC")
mt = hl.variant_qc(mt)


print("")
print("Keep common (MAF>1%) variants for pruning")
mt_chrX = mt.filter_rows((mt['locus'].contig == 'chrX') & (mt.variant_qc.AF[1] >= 0.01) & (mt.variant_qc.AF[1] <= 0.99) & (mt.variant_qc.call_rate >= 0.98)).persist()
mt_auto = mt.filter_rows(mt.locus.in_autosome() & (mt.variant_qc.AF[1] >= 0.01) & (mt.variant_qc.AF[1] <= 0.99) & (mt.variant_qc.call_rate >= 0.98)).persist()


print("")
print("Perform LD pruning")
ht_pruned_chrX = hl.ld_prune(mt_chrX.GT, r2=0.2, bp_window_size=500000)
mt_pruned_chrX = mt_chrX.filter_rows(hl.is_defined(ht_pruned_chrX[mt_chrX.row_key]))

ht_pruned_auto = hl.ld_prune(mt_auto.GT, r2=0.005, bp_window_size=10000000)
mt_pruned_auto = mt_auto.filter_rows(hl.is_defined(ht_pruned_auto[mt_auto.row_key]))


print("")
pprint('Write out pruned variants')
mt_pruned_chrX.rows().select().export(PRUNED_chrX_LIST)
mt_pruned_auto.rows().select().export(PRUNED_AUTO_LIST)

