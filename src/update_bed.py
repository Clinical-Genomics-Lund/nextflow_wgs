#!/usr/bin/env python3

import argparse
import subprocess
import os
from pathlib import Path
import re


def main(clinvar_new: Path, clinvar_old: Path, clinvardate: str, out_dir: Path, release: str, ensembl: str, skip_download: bool, incl_bed: list[str]):

    clinvar_vcf = clinvar_new
    clinvar_vcf_old = clinvar_old
    clinvardate = clinvardate
    out_dir = Path(out_dir)
    
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    release = release
    # gtf_fp = f'Homo_sapiens.GRCh38.{release}.gtf.gz'
    out_base_fp = f'{out_dir}/exons_hg38_{release}.bed'
    
    small_pad = 5
    pad = 20
    print("Write initial base")
    write_base(ensembl, out_base_fp, release, skip_download, pad)

    final_bed_path = Path(f'exons_{release}.padded{pad}bp_clinvar-{clinvardate}padded{small_pad}bp.bed')
    if final_bed_path.exists() and final_bed_path.is_file():
        print(f'Removing file: {final_bed_path}')
        final_bed_path.unlink()

    if incl_bed is not None:
        print(f'Include {incl_bed}')
        for bed_fp in incl_bed:
            suffix = bed_fp
            append_to_bed(final_bed_path, bed_fp, suffix)
    else:
        print('No extra bed files to include')

    (clinvar_new, new_benign) = read_clinvar(clinvar_vcf)
    (clinvar_old, old_benign) = read_clinvar(clinvar_vcf_old)

    clinvar_all = {**clinvar_new, **clinvar_old}

    clinvar_final_bed_fp = f'clinvar_{clinvardate}.bed'
    (new_to_add, old_to_remove) = compare_clinvar(clinvar_new, clinvar_old, clinvar_final_bed_fp, out_dir, small_pad)

    clinvar_log_path = Path('{out_dir}/clinvar_{clinvardate}.log')
    clinvar_final_bed_path = Path('{out_dir}/clinvar_{clinvardate}.bed')
    clinvar_log_path.unlink()
    clinvar_final_bed_path.unlink()

    log_changes(clinvar_log_path, clinvar_all, clinvar_new, clinvar_old)

    # log_clinvar_bed_and_info(new_to_add, 'new', clinvar_log_path, new_benign)
    # log_clinvar_bed_and_info(old_to_remove, 'old', clinvar_log_path, new_benign)
    # log_clinvar_bed_and_info(new_to_add, 'bed', clinvar_final_bed_fp, new_benign)

    append_to_bed(final_bed_path, clinvar_final_bed_fp, '.')

    sort_merge_output(final_bed_path)

    # sort_merge_output()



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

    @staticmethod
    def from_line(vcf_line) -> 'Variant':

        vcf_line = vcf_line.rstrip()
        fields = vcf_line.split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        info_str = fields[7]

        info = {}
        for info_field in info_str.split(';'):
            (key, val) = info_field.split('=')
            info[key] = val

        return Variant(chrom, pos, ref, alt, info)

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, info: dict[str, str] = {}):

        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

        self.reason = None
        
        self.clndn = 'Undefined'
        if self.info.get('CLNDN') is not None:
            self.clndn = self.info['CLNDN']
        
        self.key = 'FIXME'
        # clnacc = self.info['CLNACC'] if self.info.get('CLNACC') is not None else ""
        # self.fourth = '~'.join([self.key, self.reason, clnacc])
    
    def get_bed(self, padding: int) -> list[str, int, int, str]:
        return f'{self.chrom}\t{self.pos - padding}\t{self.pos + padding}\t<INFO PLACEHOLDER>'

    def __str__(self):
        return '\t'.join(self.fields)
    
    def get_bed_annot(self):
        reason = self.reason
        clnacc = self.info['CLNACC']
        clndn = self.info['CLNDN'] if self.info.get('CLNDN') is not None else 'Undefined'
        return f'{self.chrom}:{self.pos}_{self.ref}_{self.alt}~{reason}~{clnacc}~{clndn}'

    def get_clnsig(self) -> str:
        clnsig = self.info.get('CLNSIG')
        if clnsig != None:
            return clnsig
        else:
            return 'MISSING'
        

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
    # FIXME: Use Viktor's VCF module here?
    clinvar_variants = {}
    benign_clinvar_variants = {}

    with open(vcf_fp, 'r') as vcf_fh:
        for line in vcf_fh:
            if line.startswith('#'):
                continue
            variant = Variant.from_line(line)

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

            # print("New entry")
            # FIXME: Look over this part
            if 'Pathogenic' in sig or 'Likely_pathogenic' in sig:
                # print('First hit')
                keep = True
                reason = sig
            elif sig in 'Conflicting_interpretations_of_pathogenicity':
                if 'Pathogenic' in confidence or 'Likely_pathogenic' in confidence:
                    # print('Second hit')
                    keep = True
                    reason = confidence
            # FIXME: What is 'athogenic'
            elif 'athogenic' in haplo:
                # print('Third hit')
                keep = True
                reason = haplo

            key = f'{variant.chrom}:{variant.pos}_{variant.ref}_{variant.alt}'
            if keep:
                variant.reason = reason
                clinvar_variants[key] = variant
            else:
                variant.reason = reason
                benign_clinvar_variants[key] = variant
    return (clinvar_variants, benign_clinvar_variants)


def write_base(ensembl_fp: str, out_fp: str, release: str, skip_download: bool, padding: int):
    
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


def append_to_bed(out_bed_fp: str, bed2add_fp: str, fourth_col_default: str):
    with open(bed2add_fp, 'r') as bed2add_fh, open(out_bed_fp, 'a') as out_fh:
        for line in bed2add_fh:
            bed_entry = BedEntry(line)
            if bed_entry.annot is None:
                bed_entry.annot = fourth_col_default
            print(bed_entry, file=out_fh)


def sort_merge_output(bed_fp: str):
    tmp_bed_fp = 'tmp.sort.bed'
    sort_cmd = f'bedtools sort -i {bed_fp} > {tmp_bed_fp}'
    merge_cmd = f'bedtools merge -i {tmp_bed_fp} -c4 -o collapse > {bed_fp}'
    subprocess.call(sort_cmd, shell=True, check=True)
    subprocess.call(merge_cmd, shell=True, check=True)
    os.remove(tmp_bed_fp)


def compare_clinvar(new_clinvar: dict[str, Variant], old_clinvar: dict[str, Variant], final_bed_fp: str, padding: int, tmp_dir: str = '/tmp') -> tuple[list[str], list[str]]:
    
    new_bed_path = Path(f'{tmp_dir}/clinvar_new.bed')
    old_bed_path = Path(f'{tmp_dir}/clinvar_old.bed')
    clinvar_in_common = set()
    with open(str(new_bed_path), 'w') as new_fh:
        for key, variant in new_clinvar.items():
            if old_clinvar.get(key) is not None:
                clinvar_in_common.add(key)

            # FIXME: Look over this
            clinvar_info = variant.get_bed_annot()
            out_line = f'{variant.chrom}\t{variant.pos - padding}\t{variant.pos + padding}\t{clinvar_info}'
            print(out_line, file=new_fh)

    clinvar_new_added = set()
    for new_key in new_clinvar:
        if new_key not in old_clinvar:
            clinvar_new_added.add(new_key)

    with open(old_bed_path, 'w') as old_fh:
        clinvar_old_removed = set()
        for old_key in old_clinvar:
            if old_key not in new_clinvar:
                old_var = old_clinvar[old_key]
                clinvar_old_removed.add(old_key)
                # FIXME: Write this to old BED
                print(old_var.get_bed(padding), file=old_fh)

    new_to_add = get_bed_intersect(new_bed_path, final_bed_fp)
    old_to_remove = get_bed_intersect(old_bed_path, final_bed_fp)
    new_bed_path.unlink()
    old_bed_path.unlink()

    print(f'Clinvar in common between versions: {len(clinvar_in_common)}')
    print(f'Added new (unique targets): {len(clinvar_new_added)} ({len(new_to_add)})')
    print(f'Removed old (unique targets): {len(clinvar_old_removed)} ({len(old_to_remove)})')

    return (new_to_add, old_to_remove)


def get_bed_intersect(clinvar: str, bed: str) -> list[str]:
    not_in_bed_cmd = ['bedtools', 'intersect', '-a', clinvar, '-b', bed, '-v']
    # FIXME: How to gather and return results
    result = subprocess.run(not_in_bed_cmd, capture_output=True, text=True, check=True)
    intersected_regions = result.stdout.strip().splitlines()
    return intersected_regions


def log_changes(log_fp: str, clinvar_all: dict[str, Variant], added_keys: set[str], removed_keys: set[str], new_benign: dict[str, Variant]):
    '''
    Log variants added and removed between clinvar versions
    Check if part of latest ClinVar with benign status if removed
    '''
    with open(log_fp, 'w') as log_fh:
        for clin_key, clin_var in clinvar_all.items():
            if clin_key in added_keys:
                clnsig = 'MISSING'
                if clin_key in new_benign:
                    clnsig = new_benign[clin_key].get_clnsig()
                print(f'REMOVED: {clin_var.pos} reason[2] reason[1] => {clnsig}', file=log_fh)
                
            elif clin_key in removed_keys:
                print(f'ADDED: {clin_var.pos} reason[2] reason[1]', file=log_fh)


def append_clinvar_to_bed(bed_fp: str, clinvar_new: dict[str, Variant]):
    with open(bed_fp, 'w') as bed_fh:
        for clin_var in clinvar_new.values():
            print(f'{clin_var.chrom}\t{clin_var.pos}\t{clin_var.pos}\tCLINVAR-{clin_var.get_clnsig()}', file=bed_fh)


if __name__ == "__main__":
    args = parse_arguments()
    main(
        clinvar_new=Path(args.new),
        clinvar_old=Path(args.old),
        clinvardate=args.clinvardate,
        out_dir=Path(args.out_dir),
        release=args.release,
        ensembl=args.ensembl,
        skip_download=args.skip_download,
        incl_bed=args.incl_bed
    )


    # parser.add_argument('--old', required=True)
    # parser.add_argument('--new', required=True)
    # parser.add_argument('--release', required=True)
    # parser.add_argument('--clinvardate', required=True)

    # parser.add_argument('--out_dir', required=True)

    # parser.add_argument('--ensembl')

    # parser.add_argument('--incl_bed', nargs="*")
    # parser.add_argument('--skip_download', action="store_true")
# def main(new: str, old: str, clinvardate: str, out_dir: str, release: str, ensembl: str, skip_download: bool, incl_bed: list[str]):
