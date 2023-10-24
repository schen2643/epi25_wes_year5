import hail as hl
hl.init(log='/tmp/hail_pca_1kg.log')

from pprint import pprint
# Input
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc_hardcalls.mt'
MT_1KG = 'gs://hail-datasets-us/1000_Genomes/NYGC_30x/GRCh38/autosomes_phased.mt'

SAMPLE_INITQC_LIST = 'gs://epi25/wes/qc_files/02_initQC_sample.remove.list'
PRUNED_AUTO_LIST = 'gs://epi25/wes/qc_files/03_pruned_autosomal_variant.keep.list'

# Output
oneKG_POP_TABLE = 'gs://epi25/wes/misc-data/samples_1kg.ht'
PCA_1KG_SCORES_TABLE = 'gs://epi25/wes/qc_files/05_pca_w_1kg.tsv'



###
print("")
print("Read in epi25 data")
# Read in data
mt = hl.read_matrix_table(MT)
ht_sqc = (hl.import_table(SAMPLE_INITQC_LIST, impute=True).key_by('s'))
ht_pruned_auto = (hl.import_table(PRUNED_AUTO_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')


print("")
print("Retain samples that pass init QC")
mt = mt.filter_cols(hl.is_defined(ht_sqc[mt.s]), keep=False)


print("")
print("Load 1kG mt")
mt_1kg = hl.read_matrix_table(MT_1KG)
mt_1kg = hl.split_multi_hts(mt_1kg)


print("")
print("Empty the column fields in 1kG mt")
mt_1kg = mt_1kg.select_entries("GT")
mt_1kg.cols().write(output=oneKG_POP_TABLE, overwrite=True)
mt_1kg = mt_1kg.select_cols()


print("")
print("Join epi25 and 1kg mt")
mt_join = mt.union_cols(mt_1kg)


print("")
print("Keep pruned-in SNPs")
mt_join = mt_join.filter_rows(hl.is_defined(ht_pruned_auto[mt_join.row_key]))


print("")
print("Run PCA")
eigenvalues, scores, _ = hl.hwe_normalized_pca(mt_join.GT, k=20)


print("")
print("Export PCA with 1kG table")
scores.export(PCA_1KG_SCORES_TABLE)
