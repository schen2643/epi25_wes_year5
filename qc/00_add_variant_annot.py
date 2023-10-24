import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
from gnomad.utils.vep import process_consequences #see gnmoad python package documentation here: https://broadinstitute.github.io/gnomad_methods/examples/vep.html
hl.plot.output_notebook()


# Input
MT = 'gs://epi25/wes/genetic_data/epi25_cc_wes_hg38_alqc_split_gqc.mt'
HT_VEP = 'gs://epi25/wes/genetic_data/vep_annot.ht'
HT_MPC = 'gs://epi25/wes/annotations/fordist_constraint_official_mpc_values_v2_grch38.ht'
HT_GNOMAD_NONPSYCH = 'gs://epi25/wes/annotations/gnomad.exomes.r2.1.1.nfe_nonpsych_sites_grch38.ht'
HT_DISCOVEHR = 'gs://epi25/wes/annotations/DiscovEHR_GHS_Freeze_50.L3DP10.pVCF.frq_sites_grch38.ht'

# Output
HT_ANNOT = 'gs://epi25/wes/annotations/sites_annot.ht'



print('Read in data')
mt = hl.read_matrix_table(MT)

ht_vep = hl.read_table(HT_VEP)
ht_vep = ht_vep.select('vep')

ht_mpc = hl.read_table(HT_MPC)
ht_gnomad_nonpsych = hl.read_table(HT_GNOMAD_NONPSYCH)
ht_discov = hl.read_table(HT_DISCOVEHR)


print('Annotate with the vep information')
mt = mt.annotate_rows(vep_annot = ht_vep[mt.row_key])
mt = mt.annotate_rows(vep = mt.vep_annot.vep)
mt = mt.drop(mt.vep_annot)

mt = process_consequences(mt)
# this adds worst_consequence_term, worst_csq_for_variant, worst_csq_by_gene and other fields to ds.vep
# and returns an MT with better formatted consequences


# Case builder function
PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"]
MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant",
                 "splice_region_variant"]
SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]
OTHER_CSQS = ["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

def annotation_case_builder(worst_csq_for_variant_canonical_expr, lof_use_loftee: bool = True, mis_use_polyphen_and_sift: bool = False, mis_use_strict_def: bool = False, syn_use_strict_def: bool = False):
    case = hl.case(missing_false=True)
    if lof_use_loftee:
        case = (case
                .when(worst_csq_for_variant_canonical_expr.lof == 'HC', 'pLoF') #predicted loss-of-function
                .when(worst_csq_for_variant_canonical_expr.lof == 'LC', 'LC'))
    else:
        case = case.when(hl.set(PLOF_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'pLoF')
    if mis_use_polyphen_and_sift:
        case = (case
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence) &
                      (worst_csq_for_variant_canonical_expr.polyphen_prediction == "probably_damaging") &
                      (worst_csq_for_variant_canonical_expr.sift_prediction == "deleterious"), "damaging_missense")
                .when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), "other_missense"))
    else:
        if mis_use_strict_def:
            case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'missense_variant', 'missense')
        else:
            case = case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'missense')
    if syn_use_strict_def:
        case = case.when(worst_csq_for_variant_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
    else:
        case = case.when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'synonymous')
    case = case.when(hl.set(OTHER_CSQS).contains(worst_csq_for_variant_canonical_expr.most_severe_consequence), 'non_coding')
    return case.or_missing()

mt = mt.annotate_rows(consequence_category = annotation_case_builder(mt.vep.worst_csq_for_variant_canonical, False, True, False, False))
# pprint(mt.aggregate_rows(hl.agg.counter(mt.consequence_category)))


print('Add other annotations')
mt = mt.annotate_rows(mpc = ht_mpc[mt.row_key])
mt = mt.annotate_rows(MPC = mt.mpc.MPC)
mt = mt.drop(mt.mpc)
mt = mt.annotate_rows(inGnomAD_nonpsych = hl.is_defined(ht_gnomad_nonpsych[mt.row_key]))
mt = mt.annotate_rows(inDiscovEHR = hl.is_defined(ht_discov[mt.row_key].info))

# Annotate whether a site is a SNP, insertion, deletion, or none; and inframe insertion or deletion
mt = mt.annotate_rows(
    type = (hl.case()
    .when((hl.len(mt.alleles[0]) == 1) & (hl.len(mt.alleles[1]) == 1), "SNP")
    .when(hl.len(mt.alleles[0]) < hl.len(mt.alleles[1]), "Insertion")
    .when(hl.len(mt.alleles[0]) > hl.len(mt.alleles[1]), "Deletion")
    .or_missing()),
    infrIndel = ((mt.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_insertion") | (mt.vep.worst_csq_for_variant_canonical.most_severe_consequence == "inframe_deletion")),
    )


print('Write the annotaiton table to file')
mt_rows = mt.rows().repartition(64)
mt_rows.write(HT_ANNOT, overwrite=True)

