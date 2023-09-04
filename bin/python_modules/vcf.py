import re
import subprocess

def parse_vcf(file_name):
    # FIXME: Handle gzipped
    with open(file_name) as vcf_fh:

        header_lines = list()
        vcf_meta = dict()

        for line in vcf_fh:
            line = line.rstrip()
            if line == "":
                next
        
            if (line.startswith("#")):
                header_lines.push(line)
            
            if (line.startswith("##")):
                (my_type, meta) = parse_metainfo(line)
                if (type is not None):
                    vcf_meta[my_type]['ID'] = meta



# FIXME
def parse_metainfo(line):

    comment = re.sub('^##', '', line)

    [line_type, data] = comment.split('=')

    valid_linetypes = ['FORMAT', 'INFO', 'SAMPLE', 'FILTER']

    if line_type in valid_linetypes:
        trimmed_data = remove_surrounding(data, '<', '>')
        pairs = keyval(trimmed_data, '=')
        return (line_type, pairs)

    return (null, null)

def parse_variant():
    pass

def parse_genotype():
    pass

def parse_info():
    pass

def parse_VEP_CSQ(CSQ_var, CSQ_meta):
    CSQ_meta["Description"] = 

# Remove character(s) defined in arg2 if first in string, and arg3 if last in string
def remove_surrounding(line, before, after):
    front_trimmed = re.sub(f"^{before}", "", line)
    back_trimmed = re.sub(f"{after}$", "", line)
    return back_trimmed

# Parse string with key value pairs. Return dictionary.
# - Keys and values separated by 2nd argument
# - Pairs separated by 3rd argument
# - Handles commas in values if surrounded by double quotes
def keyval(str, keyval_sep, pair_sep):
    pair_str = str.split(pair_sep)
    assert(len(pair_str) == 2)

    pairs = dict()

    for pair in pair_str:
        if keyval_sep is not None:
            (key, val) = keyval_sep.split(" ")
            pairs[key] = val
        else:
            pairs[pair] = 1
    return pairs

def excel_float(val):
    if val == ".":
        return 0
    
    val = re.sub(".", ",", val)
    return val

def is_gzipped(file_name):
    completed_process = subprocess.run(
        f"file -L {file_name}",
        shell=True,
        stdout=subprocess.PIPE,
        text=True
    )
    output = completed_process.stdout
    is_gzipped = output.find("gzip compressed") != -1
    return is_gzipped

