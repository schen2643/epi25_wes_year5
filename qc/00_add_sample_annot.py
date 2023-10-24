import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint


# Input
MT = 'gs://epi25/wes/genetic_data/epi25_cc_wes_hg38_alqc_split_gqc.mt'
PHENO = 'gs://epi25/wes/annotations/epi25_pheno.tsv'

# Output
HT_ANNOT = 'gs://epi25/wes/annotations/samples_annot.ht'



print('Read in data')
mt = hl.read_matrix_table(MT)


print('Add sample annotations')
ht_pheno = hl.import_table(PHENO, impute=True, key="sample")
mt = mt.annotate_cols(pheno = ht_pheno[mt.s])
# mt.describe()


# Annotate and remove samples that are not to be used in analysis
mt = mt.annotate_cols(isExclude = ( hl.is_missing(mt.pheno.case_control) | ((mt.pheno.case_control != "case") & (mt.pheno.case_control != "control"))))
# mt.describe()

print('n cases and controls')
pprint(mt.aggregate_cols(hl.agg.counter(mt.pheno.case_control)))

print('missing cases and controls')
pprint(mt.aggregate_cols(hl.agg.counter(hl.is_missing(mt.pheno.case_control))))

print('epilepsy types')
pprint(mt.aggregate_cols(hl.agg.counter(mt.pheno.epilepsy_type)))
print('n to be excluded')
pprint(mt.aggregate_cols(hl.agg.counter(mt.isExclude)))


print('Write the annotaiton table to file')
mt_cols = mt.cols()
mt_cols.write(HT_ANNOT, overwrite=True)

