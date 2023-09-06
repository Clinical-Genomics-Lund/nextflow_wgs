import sys
import pprint
import python_modules.vcf as vcf

pp = pprint.PrettyPrinter(depth=4)

maxentscan = {
    'MES-NCSS_downstream_donor': 1,
    'MES-NCSS_downstream_acceptor': 1,
    'MES-NCSS_upstream_acceptor': 1,
    'MES-NCSS_upstream_donor': 1,
    'MES-SWA_acceptor_alt': 1,
    'MES-SWA_acceptor_diff': 1,
    'MES-SWA_acceptor_ref': 1,
    'MES-SWA_acceptor_ref_comp': 1,
    'MES-SWA_donor_alt': 1,
    'MES-SWA_donor_diff': 1,
    'MES-SWA_donor_ref': 1,
    'MES-SWA_donor_ref_comp': 1,
    'MaxEntScan_alt': 1,
    'MaxEntScan_diff': 1,
    'MaxEntScan_ref': 1
}


clinmod = {
    'Pathogenic': '_5_',
    'Likely_pathogenic': '_4_',
    'Likely_benign': '_3_',
    'Benign': '_2_',
    'Uncertain_significance': '_0_',
    'not_provided': '_1_',
    'drug_response': '_6_'
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
    '##INFO=<ID=GNOMADAF,Number=1,Type=Float,Description="Average AF GnomAD">',
    '##INFO=<ID=GNOMADAF_MAX,Number=1,Type=Float,Description="Highest reported AF in gnomAD">',
    '##INFO=<ID=GNOMADPOP_MAX,Number=1,Type=Float,Description="Population of highest AF">',
    '##INFO=<ID=dbNSFP_GERP___RS,Number=1,Type=Float,Description="GERP score">',
    '##INFO=<ID=dbNSFP_phyloP100way_vertebrate,Number=1,Type=Float,Description="phyloP100 score">',
    '##INFO=<ID=dbNSFP_phastCons100way_vertebrate,Number=1,Type=Float,Description="phastcons score">',
    '##INFO=<ID=CLNSIG_MOD,Number=.,Type=String,Description="Modified Variant Clinical Significance, for genmod score _0_ - Uncertain significance, _1_ - not provided, _2_ - Benign, _3_ - Likely benign, _4_ - Likely pathogenic, _5_ - Pathogenic, _6_ - drug response, _7_ - histocompatibility, _255_ - other">',
    '##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description="Most severe genomic consequence.">',
    '##INFO=<ID=CADD,Number=.,Type=String,Description="CADD phred score">',
    '##INFO=<ID=nhomalt,Number=.,Type=Integer,Description="number of alt allele homozygous individuals in gnomad">',
]


def debug(text, debug_info=''):
    if debug_info == '':
        print(f'DEBUG: {text}', file=sys.stderr)
    else:
        print(f'DEBUG: {debug_info} {text}', file=sys.stderr)


def main():

    if len(sys.argv) < 2:
        print(f'One argument is needed: ./modify_vcf_scout.py <vcf file>')
        sys.exit(1)

    headers = list()
    vcf_meta = dict()
    vcf_data = list()
    vep_csq = ''

    hits = 0

    with open(sys.argv[1]) as VEP:
        for line in VEP:
            line = line.rstrip()

            # Print and store Meta-info
            if (line.startswith("##")):

                [header_type, meta] = vcf.parse_metainfo(line)
                if header_type is not None:
                    hits += 1
                    if header_type not in vcf_meta:
                        vcf_meta[header_type] = dict()
                    vcf_meta[header_type][meta['ID']] = meta

                if line.startswith('##INFO=<ID=CSQ,Number/'):
                    vep_csq = line

            # Print and store header
            elif (line.startswith("#")):
                print('\n'.join(info_lines))
                header_line = line
                headers = header_line[1:].split("\t")
                print(header_line)

            # Print and store variant information
            # Add gnomadg
            # Add conservation scores
            else:
                # print_variant_information(var_str, header, vcf_meta, vep_csq)
                # def print_variant_information(var_str: str, header: list, vcf_meta: dict, vep_csq: str):

                parsed_variant = vcf.parse_variant(line, headers, vcf_meta)

                variant_fields = line.split("\t")
                info_fields = variant_fields[7].split(';')
                # original_variant_fields = line.split(";")
                # info_field.append(variant_fields[7])

                # # FIXME: What is this? Mitochondrial?
                if (parsed_variant["CHROM"].startswith("M")):
                    mit_info = parse_mitochondrial_chr(
                        parsed_variant, vep_csq, info_fields)
                    info_fields.extend(mit_info)

                # print('\t'.join(variants[0:7]))
                # print('\t')

                csq = parsed_variant['INFO']['CSQ'][0]

                # debug(csq, 'CSQ')

                additional_info = gnomad_info(csq, parsed_variant)
                info_fields.extend(additional_info)

                output_fields = list()
                output_fields.extend(variant_fields[0:7])
                output_fields.append(';'.join(info_fields))
                output_fields.extend(variant_fields[8:])

                print('\t'.join(output_fields))

                # debug(info_field, 'info_field')

                # base_info = original_variant_fields[0].split('\t')

                # header = original_variant_fields
                # header_output = list()
                # header_output.extend(base_info[0:7])
                # info_field = [base_info[7]]
                # info_field.extend(original_variant_fields[1:])
                # info_field.extend(additional_info)
                # header_output.append(';'.join(info_field))
                # # header_output.extend(original_variant_fields[1:-1])
                # header_output.append(original_variant_fields[-1])

                # header_output.extend(info_field + ';'.join(additional_info))
                # header_output.extend(original_variant_fields[8:])

                # print('\t'.join(header_output))

                # info_fields.extend(additional_info)
                # output_variant_line = ';'.join(info_fields)
                # print(output_variant_line)
                # print('\t')
                # print('\t'.join(variants[8:]))


def parse_mitochondrial_chr(parsed_variant: dict, vep_csq: str, info_field: list[str]) -> list[str]:
    # FIXME: What is this? Mitochondrial?
    field_names_str = vep_csq.replace(
        "Consequence annotations from Ensembl VEP. Format: ", "")
    field_names = field_names_str.split("|")
    info_field_mt = ""
    trans_c = 0

    for var_csq in [var['INFO']['CSQ'] for var in parsed_variant]:
        csq_mt = list()
        for key in field_names:
            if maxentscan[key]:
                csq_mt.append("")
            elif (key == 'Consequence'):
                # FIXME: What is this?
                tmps = parsed_variant['INFO']['CSQ'][trans_c][key]
                csq_mt.append(tmps.join('&'))
            else:
                # FIXME: What is this?
                tmps = parsed_variant['INFO']['CSQ'][trans_c]
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
    for info in info_field:
        if info.find("CSQ") != -1:
            tmpinfo.append(f'CSQ={info_field_mt}')
        else:
            tmpinfo.append(info)
    # info_fields = tmpinfo
    tmpinfo.append("GeneticModels=mt")
    return tmpinfo


def gnomad_info(csq: dict, parsed_variant: dict) -> list[str]:
    # GNOMAD
    # Overall
    additional_info = list()
    my_max = csq.get('gnomADg_AF_popmax')
    if (my_max is not None):
        additional_info.append(f'GNOMADAF_MAX={my_max}')

    # Population with max
    # max_pop = parsed_variant_doobi
    max_pop = csq.get('gnomADg_popmax')
    if max_pop is not None:
        additional_info.append(f'GNOMADPOP_MAX={max_pop}')

    # GERP
    gerp = csq.get('GERP')
    if gerp is not None:
        additional_info.append(f'dbNSFP_GERP___RS={gerp}')

    # PHASTCONS
    pC = csq.get('phastCons')
    if pC is not None:
        additional_info.append(
            f'dbNSFP_phastCons10way_vertebrate={pC}')

    # PHYLOP
    pP = csq.get('phyloP100way')
    if pP is not None:
        additional_info.append(
            f'dbNSFP_phyloP100way_vertebrate={pP}')

    # CADD
    cadd = csq.get('CADD_PHRED')
    if cadd is not None:
        additional_info.append(f'CADD={cadd}')

    # CLINSIG MODIFY
    if 'CLNSIG' in parsed_variant['INFO']:
        clinsig = get_clinsig(parsed_variant)
        additional_info.append(f'CLNSIG_MOD={clinsig}')

    # MOST SEVERE CONSEQUENCE
    # if 'Consequence' in parsed_variant['INFO']:
    most_severe = most_severe_consequence(parsed_variant)
    additional_info.append(
        f'most_severe_consequence={most_severe}')

    return additional_info


def get_clinsig(parsed_variant: dict) -> str:
    csM = parsed_variant['INFO']['CLNSIG']
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
    return joined_mods
    # add_info_field.push(f'CLNSIG_MOD={joined_mods}')


def most_severe_consequence(parsed_variant: dict) -> str:
    csq_ref = parsed_variant['INFO']['CSQ']
    m_s_c = CSQ(csq_ref)
    most_severe = '.'
    if m_s_c is not None:
        m_s_c_low = [el.lower() for el in m_s_c]
        most_severe = sorted(m_s_c_low, key=lambda el: rank[el])[0]
    return most_severe
    # add_info_field.push(f'most_severe_consequence={most_severe}')


def CSQ(consequences: dict) -> list[str]:
    all_csq = list()
    for csq_dict in consequences:
        # If canon pick consensus consequence or if equal most severe
        consequence = csq_dict['Consequence']
        # for conq in consequence:
        all_csq.append(consequence)
    return all_csq


if __name__ == "__main__":
    main()
