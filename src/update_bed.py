#!/usr/bin/env python3

import argparse
import subprocess
import os
from pathlib import Path
import re


def main():
    args = parse_arguments()

    clinvar_vcf = args.new
    clinvar_vcf_old = args.old
    clinvardate = args.clinvardate
    out_dir = Path(args.out_dir)
    
    if not out_dir.exists():
        out_dir.mkdir(out_dir, parents=True, exist_ok=True)

    release = args.release
    # gtf_fp = f'Homo_sapiens.GRCh38.{release}.gtf.gz'
    out_base_fp = f'{out_dir}/exons_hg38_{release}.bed'
    
    small_pad = 5
    pad = 20
    print("Write initial base")
    write_base(args.ensembl, out_base_fp, release, args.skip_download, pad)

    final_bed_path = Path(f'exons_{release}.padded{pad}bp_clinvar-{args.clinvardate}padded{small_pad}bp.bed')
    if final_bed_path.exists() and final_bed_path.is_file():
        print(f'Removing file: {final_bed_path}')
        final_bed_path.unlink()

    if args.incl_bed is not None:
        print(f'Include {args.incl_bed}')
        for bed_fp in args.incl_bed:
            suffix = bed_fp
            append_to_bed(final_bed_path, bed_fp, suffix)
    else:
        print('No extra bed files to include')

    (clinvar_new, new_benign) = read_clinvar(clinvar_vcf)
    (clinvar_old, old_benign) = read_clinvar(clinvar_vcf_old)

    clinvar_final_bed_fp = f'clinvar_{clinvardate}.bed'
    (new_to_add, old_to_remove) = compare_clinvar(clinvar_new, clinvar_old, clinvar_final_bed_fp, out_dir)



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--old', required=True)
    parser.add_argument('--new', required=True)
    parser.add_argument('--release', required=True)
    parser.add_argument('--clinvardate', required=True)

    parser.add_argument('--out_dir', required=True)

    parser.add_argument('--ensembl')

    parser.add_argument('--incl_bed', nargs="*")
    parser.add_argument('--skip_download', action="store_true")
    args = parser.parse_args()
    return args


class Variant:
    def __init__(self, vcf_line):
        fields = vcf_line.split('\t')
        self.chr = fields[0]
        self.pos = int(fields[1])
        self.ref = fields[3]
        self.alt = fields[4]
        info_str = fields[7]
        self.reason = None

        self.info = {}
        for info_field in info_str.split(';'):
            (key, val) = info_field.split('=')
            self.info[key] = val
        
        self.clndn = 'Undefined'
        if self.info.get('CLNDN') is not None:
            self.clndn = self.info['CLNDN']
        
        self.key = 'FIXME'
        # clnacc = self.info['CLNACC'] if self.info.get('CLNACC') is not None else ""
        # self.fourth = '~'.join([self.key, self.reason, clnacc])

        

#     fourth = [key, variant.reason, variant.info['CLNACC']]



class BedEntry:
    def __init__(self, line: str):
        line = line.rstrip()
        fields = line.split('\t')
        self.chr = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.annot = None if len(fields) == 3 else fields[3]

    def __str__(self) -> str:
        out_fields = [self.chr, str(self.start), str(self.end)]
        if self.annot is not None:
            out_fields.append(self.annot)
        return '\t'.join(out_fields)


class GtfEntry:
    def __init__(self, line: str):
            line = line.rstrip()
            fields = line.split('\t')
            self.chr = fields[0]
            self.molecule = fields[2]
            self.start_pos = int(fields[3])
            self.end_pos = int(fields[4])
            self.annotation = fields[8]


def read_clinvar(vcf_fp: str) -> tuple[dict[str, Variant], dict[str, Variant]]:
    # FIXME: Use Viktor's module here?
    clinvar_variants = {}
    benign_clinvar_variants = {}

    with open(vcf_fp, 'r') as vcf_fh:
        for line in vcf_fh:
            if line.startswith('#'):
                continue
            variant = Variant(line)

            sig = ''
            confidence = ''
            haplo = ''

            info = variant.info

            if info.get('CLNSIG') is not None:
                sig = info['CLNSIG']
            if info.get('CLNSIGINCL') is not None:
                haplo = info['CLNSIGINCL']
            if info.get('CLNSIGCONF') is not None:
                confidence = info['CLNSIGCONF']
            
            # FIXME: Annotate why
            if sig is None and haplo is not None:
                continue

            keep = False
            reason = None

            if sig.find('Pathogenic') != -1 or sig.find('Likely_pathogenic') != -1:
                keep = True
                reason = sig
            elif sig.find('Conflicting_interpretations_of_pathogenicity'):
                if confidence.find('Pathogenic') != -1 or confidence.find('Likely_pathogenic'):
                    keep = True
                    reason = confidence
            # FIXME: What is 'athogenic'
            elif haplo.find('athogenic'):
                keep = True
                reason = haplo


            key = f'{variant.chr}:{variant.pos}_{variant.ref}_{variant.alt}'
            if keep:
                variant.reason = reason
                clinvar_variants[key] = variant
            else:
                variant.reason = reason
                benign_clinvar_variants[key] = variant
    return (clinvar_variants, benign_clinvar_variants)


def write_base(ensembl_fp: str|None, out_fp: str, release: str, skip_download: bool, padding: int):
    
    if not skip_download:
        gtf_request = f'https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz'
        tmp_fp = 'tmp.gtf.gz'
        download_gtf_cmd = ['wget', gtf_request, '-O', tmp_fp]
        subprocess.run(download_gtf_cmd, check=True)
        gunzip_cmd = ['gunzip', 'tmp.gtf.gz']
        subprocess.run(gunzip_cmd, check=True)
    else:
        # FIXME: Handle situation where not available
        pass
    
    with open(ensembl_fp, 'r') as gtf_fh, open(out_fp, 'w') as out_fh:
        keep = False
        for line in gtf_fh:
            if line.startswith('#'):
                continue
            gtf_entry = GtfEntry(line)

            # Filter out non chr chromosomes
            # Added "replace" to allow "chr" based reference
            chrom = re.sub("^chr", "", gtf_entry.chr)
            if len(chrom) > 2:
                continue
            if gtf_entry.molecule == 'transcript':
                keep = gtf_entry.annotation.find('transcript_biotype "protein_coding"') != -1
            if gtf_entry.molecule == 'exon' and keep:
                out_line = f'{gtf_entry.chr}\t{gtf_entry.start_pos - padding}\t{gtf_entry.end_pos + padding}'
                print(out_line, file=out_fh)


def append_to_bed(out_bed_fp: str, bed2add_fp: str, fourth_col: str):
    with open(bed2add_fp, 'r') as bed2add_fh, open(out_bed_fp, 'a') as out_fh:
        for line in bed2add_fh:
            bed_entry = BedEntry(line)
            if bed_entry.annot is None:
                bed_entry.annot = fourth_col
            print(bed_entry, file=out_fh)


def write_sorted_merged_output(bed_fp: str):
    tmp_bed_fp = 'tmp.sort.bed'
    sort_cmd = f'bedtools sort -i {bed_fp} > {tmp_bed_fp}'
    merge_cmd = f'bedtools merge -i {tmp_bed_fp} -c4 -o collapse > {bed_fp}'
    subprocess.call(sort_cmd, shell=True, check=True)
    subprocess.call(merge_cmd, shell=True, check=True)
    os.remove(tmp_bed_fp)


def compare_clinvar(new_clinvar: dict[str, Variant], old_clinvar: dict[str, Variant], final_bed_fp: str, out_dir: str) -> tuple[list[str], list[str]]:
    
    new_bed_fp = f'{out_dir}/clinvar_new.bed'
    old_bed_fp = f'{out_dir}/clinvar_old.bed'
    clinvar_in_common = set()
    padding = 5
    with open(new_bed_fp, 'w') as new_fh, open(old_bed_fp, 'w') as old_fh:
        for key, variant in new_clinvar.items():
            if old_clinvar.get(key) is not None:
                clinvar_in_common.add(key)
            # clinvar_info = get_clinvar_info(variant, key)
            clinvar_info = variant.fourth
            out_line = f'{variant.chr}\t{variant.pos - padding}\t{variant.pos + padding}\t{clinvar_info}'
            print(out_line, file=new_fh)

    # FIXME: This is a work in progress


def get_intersect(clinvar: str, bed: str) -> list[str]:
    not_in_bed_cmd = ['bedtools', 'intersect', '-a', clinvar, '-b', bed, '-v']
    # FIXME: How to gather and return results
    subprocess.call(not_in_bed_cmd)
    return []


# def get_clinvar_info(variant: Variant, key: str) -> str:
#     fourth = [key, variant.reason, variant.info['CLNACC']]
#     clndn = 'Undefined'
#     if variant.info.get('CLNDN') is not None:
#         clndn = variant.info.get('CLNDN')
#     fourth.append(clndn)
#     return '~'.join(fourth)

def log_clinvar_bed_and_info(clinvar: dict[str, Variant], new_or_old: dict[str, Variant], clinvar_log: str):
    with open(clinvar_log, 'a') as log_fh:
        for key, variant in clinvar.items():
            clinvarreason = variant.reason
            


if __name__ == "__main__":
    main()


