import hail as hl
hl.init(log='/tmp/hail_write_qced_mt.log')

from pprint import pprint

# Input
MT = 'gs://epi25/wes/genetic_data/epi25_cc_wes_hg38_alqc_split_gqc.mt'

# combined sample list across pops
SAMPLE_LIST_QCed = 'gs://epi25/wes/qc_files/12_pops_postQC_sample.keep.list'
MT_EUR = 'gs://epi25/wes/genetic_data/epi25_cc_wes_split_annot_qc_eur.mt'
MT_AFR = 'gs://epi25/wes/genetic_data/epi25_cc_wes_split_annot_qc_afr.mt'
MT_EAS = 'gs://epi25/wes/genetic_data/epi25_cc_wes_split_annot_qc_eas.mt'
MT_SAS = 'gs://epi25/wes/genetic_data/epi25_cc_wes_split_annot_qc_sas.mt'
MT_FIN = 'gs://epi25/wes/genetic_data/epi25_cc_wes_split_annot_qc_fin.mt'
MT_AMR = 'gs://epi25/wes/genetic_data/epi25_cc_wes_split_annot_qc_amr.mt'

# Output
MT_QCed = 'gs://epi25/wes/genetic_data/epi25_cc_wes_split_annot_qc_joints.mt'



print("")
print("Read in data")
mt = hl.read_matrix_table(MT)

ht_sample = hl.import_table(SAMPLE_LIST, impute=True, no_header=True)
ht_sample = ht_sample.annotate(s = ht_sample.f0).key_by('s')

ht_vqc_eur = hl.read_matrix_table(MT_EUR).rows()
ht_vqc_afr = hl.read_matrix_table(MT_AFR).rows()
ht_vqc_eas = hl.read_matrix_table(MT_EAS).rows()
ht_vqc_sas = hl.read_matrix_table(MT_SAS).rows()
ht_vqc_fin = hl.read_matrix_table(MT_FIN).rows()
ht_vqc_amr = hl.read_matrix_table(MT_AMR).rows()

print("")
print("Keep variants that pass final QC")
mt = mt.filter_rows( 
    (hl.is_defined(ht_vqc_eur[mt.locus, mt.alleles])) |  
    (hl.is_defined(ht_vqc_afr[mt.locus, mt.alleles])) |
    (hl.is_defined(ht_vqc_eas[mt.locus, mt.alleles])) |
    (hl.is_defined(ht_vqc_sas[mt.locus, mt.alleles])) | 
    (hl.is_defined(ht_vqc_fin[mt.locus, mt.alleles])) |
    (hl.is_defined(ht_vqc_amr[mt.locus, mt.alleles]))
    )

print("")
print("Keep samples that pass final QC")
mt = mt.filter_cols(hl.is_defined(ht_sample[mt.s]))


print("Perform sample and variant QC on the final data set")
mt = hl.variant_qc(mt)


print("Write out the qc'ed MT")
mt.write(MT_QCed, overwrite=True)
