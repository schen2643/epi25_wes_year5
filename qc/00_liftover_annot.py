import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()


# Hail documentation on liftover: e.g., from GRCh37 to 38
# https://hail.is/docs/0.2/guides/genetics.html#liftover-variants-from-one-coordinate-system-to-another

# Prepare for liftover
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

def flip_base(base: str) -> str:
    """
    Returns the complement of a base
    :param str base: Base to be flipped
    :return: Complement of input base
    :rtype: str
    """
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('G', 'C')
            .when('C', 'G')
            .default(base))



######################
# [ MPC score ]
# Read in data
MPC_SCORE = 'gs://epi25/wes/misc_data/fordist_constraint_official_mpc_values_v2.tsv.bgz'
HT_MPC = 'gs://epi25/wes/annotations/fordist_constraint_official_mpc_values_v2_grch38.ht'

ht_mpc = hl.import_table(MPC_SCORE, impute=True)
ht_mpc = ht_mpc.annotate(
  locus = hl.locus(contig = ht_mpc.chrom, pos = ht_mpc.pos),
    alleles = [ht_mpc.ref, ht_mpc.alt]).key_by('locus','alleles').select('MPC')

# Liftover
ht_mpc = ht_mpc.annotate(new_locus=hl.liftover(ht_mpc.locus, 'GRCh38', include_strand=True))
ht_mpc = ht_mpc.filter(hl.is_defined(ht_mpc.new_locus))

ht_mpc = ht_mpc.annotate(
        new_alleles = hl.cond(ht_mpc.new_locus.is_negative_strand,
          [flip_base(ht_mpc.alleles[0]), flip_base(ht_mpc.alleles[1])], ht_mpc.alleles)
        )
ht_mpc = ht_mpc.key_by(locus=ht_mpc.new_locus.result, alleles=ht_mpc.new_alleles)

# Write the result to file
ht_mpc.write(HT_MPC, overwrite=True)



######################
# [ gnomAD ]: The entirety of the gnomad data
HT_GNOMAD_SITES_GRCh37 = 'gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht'
HT_GNOMAD_SITES_GRCh38 = 'gs://epi25/wes/annotations/gnomad.exomes.r2.1.1.nfe_nonpsych_sites_grch38.ht'

# Read in data
ht_gnomad = hl.read_table(HT_GNOMAD_SITES_GRCh37)
# pprint(ht_gnomad.freq_index_dict.collect())

# Extract variants that are present in the nfe non-neuro subset
non_neuro_index = hl.eval(ht_gnomad.freq_index_dict.get('non_neuro_nfe'))
ht_gnomad_nonpsych = ht_gnomad.filter(ht_gnomad.freq[non_neuro_index].AC > 0)

# Liftover
ht_gnomad_nonpsych = ht_gnomad_nonpsych.annotate(new_locus=hl.liftover(ht_gnomad_nonpsych.locus, 'GRCh38', include_strand=True))
ht_gnomad_nonpsych = ht_gnomad_nonpsych.filter(hl.is_defined(ht_gnomad_nonpsych.new_locus))

ht_gnomad_nonpsych = ht_gnomad_nonpsych.annotate(
        new_alleles = hl.cond(ht_gnomad_nonpsych.new_locus.is_negative_strand,
          [flip_base(ht_gnomad_nonpsych.alleles[0]), flip_base(ht_gnomad_nonpsych.alleles[1])], ht_gnomad_nonpsych.alleles)
        )

# Write out results
ht_gnomad_nonpsych = ht_gnomad_nonpsych.key_by(locus=ht_gnomad_nonpsych.new_locus.result, alleles=ht_gnomad_nonpsych.new_alleles)
ht_gnomad_nonpsych.write(HT_GNOMAD_SITES_GRCh38, overwrite=True)



######################
# [ DiscovEHR ]
VCF_DISCOVEHR = 'gs://epi25/wes/misc_data/DiscovEHR_GHS_Freeze_50.L3DP10.pVCF.frq.vcf'
HT_DISCOVEHR = 'gs://epi25/wes/annotations/DiscovEHR_GHS_Freeze_50.L3DP10.pVCF.frq_sites_grch38.ht'

# Read in data
mt_discR = hl.import_vcf(VCF_DISCOVEHR)
# mt_discR.describe()
ht_discR = mt_discR.rows()

# Liftover
ht_discR = ht_discR.annotate(new_locus=hl.liftover(ht_discR.locus, 'GRCh38', include_strand=True))
ht_discR = ht_discR.filter(hl.is_defined(ht_discR.new_locus))

ht_discR = ht_discR.annotate(
        new_alleles = hl.cond(ht_discR.new_locus.is_negative_strand,
                              [flip_base(ht_discR.alleles[0]), flip_base(ht_discR.alleles[1])], ht_discR.alleles)
        )
ht_discR = ht_discR.key_by(locus=ht_discR.new_locus.result, alleles=ht_discR.new_alleles)

# Write out results
ht_discR.write(HT_DISCOVEHR, overwrite=True)

