import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Input
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc_hardcalls.mt'
SAMPLE_INITQC_LIST = 'gs://epi25/wes/qc_files/02_initQC_sample.remove.list'
PRUNED_chrX_LIST = 'gs://epi25/wes/qc_files/03_pruned_chrX_variant.keep.list'

# Output
IMPUTED_SEX_TABLE = 'gs://epi25/wes/qc_files/04_imputed_sex.tsv'
Y_NCALLED = 'gs://epi25/wes/qc_files/04_ycalled.tsv'
#########


print("Read in data")
mt = hl.read_matrix_table(MT)
ht_sqc = (hl.import_table(SAMPLE_INITQC_LIST, impute=True).key_by('s'))
ht_pruned_chrX = (hl.import_table(PRUNED_chrX_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')

print("")
print("Retain samples that pass init QC")
mt = mt.filter_cols(hl.is_defined(ht_sqc[mt.s]), keep=False)

print("")
print("Keep chrX pruned-in SNPs")
mt = mt.filter_rows(hl.is_defined(ht_pruned_chrX[mt.row_key]))


print("")
print("Impute sex")
ht_imputed_sex = hl.impute_sex(mt.GT, include_par=False, female_threshold=0.5, male_threshold=0.8)


print("")
print("Write out sex imputation metrics")
ht_imputed_sex.export(IMPUTED_SEX_TABLE)

mt = hl.read_matrix_table(MT)
mt = mt.filter_cols(hl.is_defined(ht_sqc[mt.s]), keep=False)
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt)

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.sample_qc.n_called).export(Y_NCALLED)

