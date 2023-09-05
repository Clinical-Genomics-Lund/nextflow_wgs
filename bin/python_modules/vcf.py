import re
import subprocess


def debug(text, debug_info=''):
    if debug_info == '':
        print(f'DEBUG: {text}')
    else:
        print(f'DEBUG: {debug_info} {text}')


def parse_vcf(file_name: str) -> tuple[dict, list, list]:

    header_lines = list()
    vcf_meta = dict()
    vcf_data = list()
    header = list()
    samples = list()

    # FIXME: Handle gzipped
    with open(file_name) as vcf_fh:

        for line in vcf_fh:
            line = line.rstrip()
            if line == '':
                continue

            if (line.startswith('#')):
                header_lines.append(line)

            if (line.startswith('##')):
                (my_type, meta) = parse_metainfo(line)
                if (my_type is not None):
                    vcf_meta[my_type]['ID'] = meta
            elif (line.startswith('##')):
                line_content = line[1:]
                header = line_content.split('\t')
            else:
                if header is None:
                    raise Exception(
                        "Malformed VCF: No column description header")
                variant = parse_variant(line, header, vcf_meta)
                if "CHROM" in variant:
                    vcf_data.append(variant)

    samples = header[9:]
    vcf_meta[header] = '\n'.join(header_lines)

    return (vcf_meta, vcf_data, samples)


# FIXME
def parse_metainfo(comment: str) -> tuple[str, dict] | tuple[None, None]:

    comment = re.sub('^##', '', comment)
    # debug(f'comment {comment}')
    fields = comment.split('=')
    line_type = fields[0]
    data = '='.join(fields[1:])
    # [line_type, data] = comment.split('=')
    valid_linetypes = ['FORMAT', 'INFO', 'SAMPLE', 'FILTER']

    # debug(f'{line_type}, {data}')

    if line_type in valid_linetypes:
        trimmed_data = remove_surrounding(data, '<', '>')
        pairs = keyval(trimmed_data, '=', ',')
        return (line_type, pairs)

    return (None, None)


def parse_variant(var_str: str, head: list[str], meta: dict[str, dict]) -> dict:
    var_data = var_str.split("\t")
    variants = dict()

    variants[var_str] = var_str

    # First seven fields
    for i in range(0, 7):
        variants[head[i]] = var_data[i]

    # Eight field, INFO
    variants["INFO"] = parse_info(var_data[7])

    # FIXME: Run with a sample containing CSQ field
    if ('CSQ' in variants["INFO"]):
        # debug(meta['INFO'], 'test')

        assert 'CSQ' in meta['INFO'], 'CSQ not found among: {meta}'.format(
            meta['INFO'])

        variants["INFO"]["CSQ"] = parse_VEP_CSQ(
            variants["INFO"]["CSQ"],
            meta["INFO"]["CSQ"]
        )

    for i in range(9, len(var_data)):
        if 'GT' not in variants:
            variants['GT'] = dict()
        variants["GT"][head[i]] = parse_genotype(var_data[8], var_data[i])

    return variants


def parse_genotype(format_str: str, data_str: str) -> dict[str, str]:
    format_arr = format_str.split(':')
    data = data_str.split(':')

    gt = dict(zip(format_arr, data))
    return gt


def parse_info(info_str: str) -> dict:
    info = keyval(info_str, '=', ';')
    return info


def parse_VEP_CSQ(CSQ_var: str, CSQ_meta: dict[str, str]) -> list[dict[str, str]]:
    field_names = CSQ_meta['Description']\
        .replace('Consequence annotations from Ensembl VEP. Format: ', '')\
        .split('|')
    transcripts = CSQ_var.split(',')
    data_transcripts = list()
    for transcript_CSQ in transcripts:
        values = transcript_CSQ.split('|')

        data = dict()
        for i in range(0, len(field_names)):
            if field_names[i] == 'Consequence':
                # FIXME: Something strange here
                conseq_array = values[i].split('&')
                # FIXME; Compare to the Perl logic
                data[field_names[i]] = conseq_array[0]
            else:
                data[field_names[i]] = values[i] if len(values) > i else ""
        data_transcripts.append(data)

    # debug(data_transcripts, 'data_transcripts')
    # sys.exit(1)

    return data_transcripts

# Remove character(s) defined in arg2 if first in string, and arg3 if last in string


def remove_surrounding(line: str, before: str, after: str) -> str:
    front_trimmed = re.sub(f"^{before}", "", line)
    back_trimmed = re.sub(f"{after}$", "", front_trimmed)
    return back_trimmed

# Parse string with key value pairs. Return dictionary.
# - Keys and values separated by 2nd argument
# - Pairs separated by 3rd argument
# - Handles commas in values if surrounded by double quotes


def keyval(my_str: str, keyval_sep: str, pair_sep: str) -> dict[str, str | int]:
    pair_strs = my_str.split(pair_sep)
    # debug(pair_strs, 'keyval1')

    pairs = dict()

    for pair in pair_strs:
        # If key-value separator exsts, save the value for the key
        if keyval_sep is not None:
            # debug(f'keyval {pair}')
            fields = pair.split(keyval_sep)
            key = fields[0]
            val = keyval_sep.join(fields[1:])
            pairs[key] = val
        # Otherwise treat the whole string as a flag and set it to one (true)
        else:
            pairs[pair] = 1
    return pairs


def excel_float(val: str) -> str:
    if val == ".":
        return '0'

    val = re.sub(".", ",", val)
    return val


def is_gzipped(file_name: str) -> bool:
    completed_process = subprocess.run(
        f"file -L {file_name}",
        shell=True,
        stdout=subprocess.PIPE,
        text=True
    )
    output = completed_process.stdout
    is_gzipped = output.find("gzip compressed") != -1
    return is_gzipped
