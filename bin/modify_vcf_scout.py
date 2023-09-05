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

info_lines = [
    "##INFO\=<ID=GNOMADAF\,Number=1\,Type=Float,Description=\"Average AF GnomAD\">",
    "##INFO=<ID=GNOMADAF_MAX,Number=1,Type=Float,Description=\"Highest reported AF in gnomAD\">",
    "##INFO=<ID=GNOMADPOP_MAX,Number=1,Type=Float,Description=\"Population of highest AF\">",
    "##INFO=<ID=dbNSFP_GERP___RS,Number=1,Type=Float,Description=\"GERP score\">",
    "##INFO=<ID=dbNSFP_phyloP100way_vertebrate,Number=1,Type=Float,Description=\"phyloP100 score\">",
    "##INFO=<ID=dbNSFP_phastCons100way_vertebrate,Number=1,Type=Float,Description=\"phastcons score\">",
    "##INFO=<ID=CLNSIG_MOD,Number=.,Type=String,Description=\"Modified Variant Clinical Significance, for genmod score _0_ - Uncertain significance, _1_ - not provided, _2_ - Benign, _3_ - Likely benign, _4_ - Likely pathogenic, _5_ - Pathogenic, _6_ - drug response, _7_ - histocompatibility, _255_ - other\">",
    "##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description=\"Most severe genomic consequence.\">",
    "##INFO=<ID=CADD,Number=.,Type=String,Description=\"CADD phred score\">",
    "##INFO=<ID=nhomalt,Number=.,Type=Integer,Description=\"number of alt allele homozygous individuals in gnomad\">",
]

# vep_csq?

def debug(text, debug_info=''):
    if debug_info == '':
        print(f'DEBUG: {text}')
    else:
        print(f'DEBUG: {debug_info} {text}')


def main():

    if len(sys.argv) < 2:
        print(f'One argument is needed: ./modify_vcf_scout.py <vcf file>')
        sys.exit(1)

    headers = list()
    vcf_meta = dict()
    vcf_data = list()
    vep_csq = ''

    with open(sys.argv[1]) as VEP:
        for line in VEP:
            line = line.rstrip()
            # print(line)

            # Print and store Meta-info
            if (line.startswith("##")):

                debug(f'Working on line {line}')

                [header_type, meta] = vcf.parse_metainfo(line)
                if header_type is not None:
                    
                    vcf_meta[header_type]["ID"] = meta
                if line.startswith("##INFO=<ID=CSQ,Number/"):
                    vep_csq = line
                # print(f"{header_type}-{meta}")

                # Convert in meta data
                # If INFO=CSQ?

            # Print and store header
            elif (line.startswith("#")):
                print('\n'.join(info_lines))
                headers = line[1:].split("\t")

            # Print and store variant information
            # Add gnomadg
            # Add conservation scores
            else:
                print_variant_information(line, headers, vcf_meta, vep_csq)

def print_variant_information(line: str, header: list, vcf_meta: dict, vep_csq: str):
    parsed_variant_doobi = vcf.parse_variant(header, vcf_meta)

    add_info_field = list()
    variants = line.split("\t")
    info_field = [line.split(";"), variants[7]]

    # FIXME: What is this? Mitochondrial?
    if (parsed_variant_doobi["CHROM"].startswith("M")):
        field_names_str = vep_csq.replace("Consequence annotations from Ensembl VEP. Format: ", "")
        field_names = field_names_str.split("|")
        info_field_mt = ""
        trans_c = 0

        for trans in [var['INFO']['CSQ'] for var in parsed_variant_doobi]:
            csq_mt = list()
            for key in field_names:
                if maxentscan[key]:
                    csq_mt.append("")
                elif (key == 'Consequence'):
                    # FIXME: What is this?
                    tmps = parsed_variant_doobi['INFO']['CSQ'][trans_c][key]
                    csq_mt.append(tmps.join('&'))
                else:
                    # FIXME: What is this?
                    tmps = parsed_variant_doobi['INFO']['CSQ'][trans_c]
                    csq_mt.append(tmps)
            csq_trans = '|'.join(csq_mt)
            # FIXME: In perl a leading | was trimmed, is this needed?

            if trans_c == 0:
                info_field_mt = info_field_mt + csq_trans
            else:
                info_field_mt = info_field_mt + "," + csq_trans

            # Next transcript
            trans_c += 1
        
        tmpinfo = list()
        for info in tmpinfo:
            if info.find("CSQ") != -1:
                tmpinfo.append(f'CSQ={info_field_mt}')
            else:
                tmpinfo.append(info)
        info_field = tmpinfo
        info_field.append("GeneticModels=mt")

    print('\t'.join(variants[0:7]))
    print('\t')

    csq = parsed_variant_doobi['INFO']['CSQ'][0]

    # GNOMAD (?)
    # Overall
    my_max = csq['gnomADg_AF_popmax']
    if (my_max is not None):
        add_info_field.append(f'GNOMADAF_MAX={my_max}')
    
    # Population with max
    # max_pop = parsed_variant_doobi
    max_pop = csq['gnomADg_popmax']
    if max_pop is not None:
        add_info_field.append(f'GNOMADPOP_MAX={max_pop}')

    # GERP
    gerp = csq['GERP']
    add_info_field.append(f'dbNSFP_GERP___RS={gerp}')

    # PHASTCONS
    pC = csq['phastCons']
    if pC is not None:
        add_info_field.append(f'dbNSFP_phastCons10way_vertebrate={pC}')

    # PHYLOP
    pP = csq['phyloP100way']
    add_info_field.append(f'dbNSFP_phyloP100way_vertebrate={pP}')

    # CADD
    cadd = csq['CADD_PHRED']
    if cadd is not None:
        add_info_field.append(f'CADD={cadd}')

    # CLINSIG MODIFY
    modify_clinsig(add_info_field, parsed_variant_doobi)

    # MOST SEVERE CONSEQUENCE
    most_severe_consequence(add_info_field, parsed_variant_doobi)

    info_field.append(add_info_field)
    print(';'.join(info_field))
    print('\t')
    print('\t'.join(variants[8:]))

# FIXME: Return the result instead and mutate add_info_field from outside
def modify_clinsig(add_info_field, doobi):
    csM = doobi['INFO']['CLNSIG']
    mods = list()
    if csM is not None:
        csMs = csM.split(',')
        for entry in csMs:
            split_slash = entry.split('/')
            for entry_slash in split_slash:
                if clinmod[entry_slash]:
                    mods.append(clinmod[entry_slash])
                else:
                    # FIXME: What is the purpose?
                    mods.append('_255_')
    joined_mods = '|'.join(mods)
    add_info_field.push(f'CLNSIG_MOD={joined_mods}')

def most_severe_consequence(add_info_field, doobi):
    csq_ref = doobi['INFO']['CSQ']
    m_s_c = CSQ(csq_ref)
    most_severe = '.'
    if m_s_c is not None:
        m_s_c_low = [el.lower() for el in m_s_c]
        most_severe = sorted(m_s_c_low, key=lambda el: rank[el])[-1]
    add_info_field.push(f'most_severe_consequence={most_severe}')

# THESE ARE ALL COUNTED AS "other":

#      4 Affects
#     14 association
#      2 _association
#     92 Conflicting_interpretations_of_pathogenicity
#      3 other
#      4 _other
#      6 protective
#      2 _protective
#     11 _risk_factor
#     37 risk_factor

# FIXME: Document purpose
def CSQ(csq):
    all_csq = list()
    for b in csq:
        # If canon pick consensus consequence or if equal most severe
        tmp = b['Consequence']
        for conq in tmp:
            all_csq.append(conq)
    return all_csq

if __name__ == "__main__":
    main()
