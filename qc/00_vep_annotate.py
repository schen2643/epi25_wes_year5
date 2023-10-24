import hail as hl
hl.init()


CONFIG = 'gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json'
MT = 'gs://epi25/wes/genetic_data/epi25_yr5_cc_wes_hg38_alqc_split_gqc.mt'
HT_VEP = 'gs://epi25/wes/genetic_data/vep_annot.ht'


ht = hl.read_matrix_table(MT).rows()
ht_vep = hl.vep(ht, CONFIG)
ht_vep.write(HT_VEP, overwrite=True)

