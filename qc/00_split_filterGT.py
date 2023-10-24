import hail as hl
hl.init()

from pprint import pprint

MT_RAW = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38.mt'

# adding in allele length (al) fitler before split_multi, and save another mt
MT_SPLIT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc.mt'
MT_HARDCALLS = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc_hardcalls.mt'



print('')
print('Read in data')
mt = hl.read_matrix_table(MT_RAW)


print('Filter out alleles with length > 6')
mt = mt.filter_rows(mt.alleles.length() > 6, keep=False)


print('Split multiallelic sites')
mt = hl.split_multi_hts(mt, left_aligned=True)
# mt.write(MT_SPLIT, overwrite=True)


print('Filter genotypes')
filter_condition_ab = hl.is_defined(mt.GT) & (
                (mt.GT.is_hom_ref() & 
                        ((mt.GQ < 20) | (mt.DP < 10))) | 
                (mt.GT.is_het() & 
                    (((mt.AD[0]+mt.AD[1])/mt.DP < 0.8) | (mt.AD[1]/mt.DP < 0.2) | (mt.PL[0] < 20) | (mt.DP < 10))) |
                (mt.GT.is_hom_var() & 
                    ((mt.AD[1]/mt.DP < 0.8) | (mt.PL[0] < 20) | (mt.DP < 10))))


print('Summarize genotype QC')
mt = mt.filter_entries(filter_condition_ab, keep=False)


print('')
print('Wrtite split-multi and geno-filtered MT')
mt.write(MT_SPLIT, overwrite=True)
mt = mt.checkpoint(MT_SPLIT, overwrite=True)


print('')
print('Save a hardcalls copy')
mt = hl.read_matrix_table(MT_SPLIT)
mt.select_entries(mt.GT).naive_coalesce(800).write(MT_HARDCALLS, overwrite=True)

