import sys
import python_modules.vcf as vcf

maxentscan = {
    "MES-NCSS_downstream_donor":1,
    "MES-NCSS_downstream_acceptor" :1,
    "MES-NCSS_upstream_acceptor":1,
    "MES-NCSS_upstream_donor":1,
    "MES-SWA_acceptor_alt":1,
    "MES-SWA_acceptor_diff":1,
    "MES-SWA_acceptor_ref":1,
    "MES-SWA_acceptor_ref_comp":1,
    "MES-SWA_donor_alt":1,
    "MES-SWA_donor_diff":1,
    "MES-SWA_donor_ref":1,
    "MES-SWA_donor_ref_comp":1,
    "MaxEntScan_alt":1,
    "MaxEntScan_diff":1,
    "MaxEntScan_ref": 1
}

clinmod = {
    "Pathogenic": "_5_",
    "Likely_pathogenic": "_4_",
    "Likely_benign": "_3_",
    "Benign": "_2_",
    "Uncertain_significance": "_0_",
    "not_provided": "_1_", 
    "drug_response": "_6_"
}

rank = {
    'transcript_ablation': 1,
    'initiator_codon_variant': 2,
    'frameshift_variant': 3,
    'stop_gained': 4,
    'start_lost': 5,
    'stop_lost': 6,
    'splice_acceptor_variant': 7,
    'splice_donor_variant': 8,
    'inframe_deletion': 9,
    'transcript_amplification': 10,
    'splice_region_variant': 11,
    'missense_variant': 12,
    'protein_altering_variant': 13,
    'inframe_insertion': 14,
    'incomplete_terminal_codon_variant': 15,
    'non_coding_transcript_exon_variant': 16,
    'synonymous_variant': 17,
    'mature_mirna_variant': 18,
    'non_coding_transcript_variant': 19,
    'regulatory_region_variant': 20,
    'upstream_gene_variant': 21,
    'regulatory_region_amplification': 22,
    'tfbs_amplification': 23,
    '5_prime_utr_variant': 24,
    'intron_variant': 25,
    '3_prime_utr_variant': 26,
    'feature_truncation': 27,
    'tf_binding_site_variant': 28,
    'start_retained_variant': 29,
    'stop_retained_variant': 30,
    'feature_elongation': 31,
    'regulatory_region_ablation': 32,
    'tfbs_ablation': 33,
    'coding_sequence_variant': 34,
    'downstream_gene_variant': 35,
    'nmd_transcript_variant': 36,
    'intergenic_variant': 37
}

# vep_csq?

def main():
    print(f"Hello world! First argument: {sys.argv[1]}")

    with open(sys.argv[1]) as VEP:
        for line in VEP:
            line = line.rstrip()
            # print(line)

            # Print and store Meta-info
            if (line.startswith("##")):

                [header_type, meta] = vcf.parse_metainfo(line)
                print(f"{header_type}-{meta}")

                # Convert in meta data
                # If INFO=CSQ?

            # Print and store header
            else if (line.startswith("#")):
                print("Header to output file here")
            # Print and store variant information
            # Add gnomadg
            # Add conservation scores
            else:
                pass

if __name__ == "__main__":
    main()
