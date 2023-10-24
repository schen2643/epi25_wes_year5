import hail as hl
hl.init(log='/tmp/hail_pca_pop.log')

from pprint import pprint
# Input
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc_hardcalls.mt'
MT_1KG = 'gs://hail-datasets-us/1000_Genomes/NYGC_30x/GRCh38/autosomes_phased.mt'

POP = "EUR" # AFR, EAS, SAS, AMR, FIN
POP_SAMPLE_LIST = 'gs://epi25/wes/qc_files/06_predicted_{0}_sample.list'.format(POP)
PRUNED_AUTO_LIST = 'gs://epi25/wes/qc_files/03_pruned_autosomal_variant.keep.list'

# Output
POP_PCA_1KG_SCORES_TABLE = 'gs://epi25/wes/qc_files/07_pca_{0}_w_1kg.tsv'.format(POP.lower())


print("")
print("Read in data")
mt = hl.read_matrix_table(MT)
ht_pop = hl.import_table(POP_SAMPLE_LIST, impute=True, no_header=True)
ht_pop = ht_pop.annotate(s = ht_pop.f0).key_by('s')
ht_pruned_auto = (hl.import_table(PRUNED_AUTO_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')


print("")
print("Retain single-pop samples")
mt_pop = mt.filter_cols(hl.is_defined(ht_pop[mt.s]))


print("")
print("Load 1kG mt")
mt_1kg = hl.read_matrix_table(MT_1KG)
mt_1kg = hl.split_multi_hts(mt_1kg)
mt_1kg = mt_1kg.select_entries("GT")
mt_1kg.describe()
if POP == "FIN":
	mt_1kg_pop = mt_1kg.filter_cols(mt_1kg.Population == POP)
else:
	mt_1kg_pop = mt_1kg.filter_cols(mt_1kg.Superpopulation == POP)
mt_1kg_pop = mt_1kg_pop.select_cols()


print("")
print("Join epi25 and 1kg mt")
mt_pop = mt_pop.union_cols(mt_1kg_pop)


print("")
print("Keep pruned-in SNPs")
mt_pop = mt_pop.filter_rows(hl.is_defined(ht_pruned_auto[mt_pop.row_key]))


print("")
print("Run PCA")
eigenval_pop, score_pop, _ = hl.hwe_normalized_pca(mt_pop.GT, k=20)


print("")
print("Export PCA with 1kG table")
score_pop.export(POP_PCA_1KG_SCORES_TABLE)

