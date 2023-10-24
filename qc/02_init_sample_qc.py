import hail as hl
hl.init(log="/tmp/hail-fail.log")

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Input
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc.mt'
VARIANT_PREQC_LIST = 'gs://epi25/wes/qc_files/01_prefilter_variant.remove.list'
PHENO = 'gs://epi25/wes/annotations/samples_annot.ht'

# Output
SAMPLE_INITQC_TABLE = 'gs://epi25/wes/qc_files/02_initQC_sample.tsv'
SAMPLE_INITQC_LIST = 'gs://epi25/wes/qc_files/02_initQC_sample.remove.list'


print("Read in data")
mt = hl.read_matrix_table(MT)
ht_vqc = (hl.import_table(VARIANT_PREQC_LIST, types={'locus': hl.tlocus('GRCh38'), 'alleles': hl.tarray('str')})).key_by('locus','alleles')


print("Remove variants that fail prefilter")
mt = mt.filter_rows(hl.is_defined(ht_vqc[mt.locus, mt.alleles]), keep=False)


print("Annotate with phenotype file")
ht_pheno = hl.read_table(PHENO)
ht_pheno.describe()
mt = mt.annotate_cols(pheno_tmp = ht_pheno[mt.s])
mt = mt.annotate_cols(pheno = mt.pheno_tmp.pheno)
mt = mt.annotate_cols(isExclude = mt.pheno_tmp.isExclude)
mt = mt.drop('pheno_tmp')


print("Keep samples to be used in analysis")
mt = mt.filter_cols(~mt.isExclude)


print("Perform sample QC")
mt = hl.sample_qc(mt)
print("Summarize sample call rate")
pprint(mt.aggregate_cols(hl.agg.stats(mt.sample_qc.call_rate)))


print("Output sample QC metrics")
ht_sqc = mt.cols()
ht_sqc = ht_sqc.select('pheno', 'sample_qc').flatten()
ht_sqc.export(SAMPLE_INITQC_TABLE)


print("Filter to samples that fail each QC metric")
mt = mt.filter_cols((mt.sample_qc.call_rate < 0.9) | (mt.sample_qc.dp_stats.mean < 25) | (mt.sample_qc.gq_stats.mean < 57) | (mt.pheno.pct_contamination > 2.5) |
    (mt.pheno.pct_chimeras > 2) | (mt.sample_qc.n_singleton > 500))


print("")
print("Write out the list of samples to remove")
mt.cols().select().export(SAMPLE_INITQC_LIST)
