#!/usr/bin/env python3

import argparse
import subprocess
import os


def main():
    args = parse_arguments()

    clinvar_vcf = args.new
    clinvar_vcf_old = args.old

    release = args.build
    gtf = f'Homo_sapiens.GRCh38.{release}.gtf.gz'
    gtf_request = f'http://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/gft'
    base_bed = f'exons_hg38_{release}.bed'

    get_base(gtf_request, gtf, base_bed, release, args.skip_download)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--old', required=True)
    parser.add_argument('--new', required=True)
    parser.add_argument('--build', required=True)
    parser.add_argument('--clinvardate', required=True)

    parser.add_argument('--incl_bed', nargs="*")
    parser.add_argument('--skip_download', action="store_true")
    args = parser.parse_args()
    return args


class Variant:
    def __init__(self, vcf_line):
        fields = vcf_line.split('\t')
        self.chr = fields[0]
        self.pos = int(fields[1])
        info_str = fields[7]
        reason = None

        info = {}
        for info_field in info_str.split(';'):
            (key, val) = info_field.split('=')
            info[key] = val


def read_clinvar(vcf_fp: str) -> tuple[dict[str, Variant], dict[str, Variant]]:
    # FIXME: Use Viktor's module here?
    clinvar_variants = {}
    benign_clinvar_variants = {}

    with open(vcf_fp, 'r') as vcf_fh:
        for line in vcf_fh:
            if line.startswith('#'):
                continue
            variant = Variant(line)

            sig = None
            confidence = None
            haplo = None

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


def get_base(gtf_req: str, gtf: str, out_fp: str, release: str, skip_download: bool):
    
    padding = 20

    if not skip_download:
        tmp_fp = 'tmp.gtf.gz'
        download_gtf_cmd = ['wget', gtf_req, '-O', tmp_fp]
        subprocess.run(download_gtf_cmd, check=True)
        gunzip_cmd = ['gunzip', 'tmp.gtf.gz']
        subprocess.run(gunzip_cmd, check=True)
    
    with open(gtf, 'r') as gtf_fh, open(out_fp, 'w') as out_fh:
        keep = False
        for line in gtf_fh:
            if line.startswith('#'):
                continue
            gtf_entry = GtfEntry(line)

            # FIXME: Guess this filters away non regular chromosomes?
            if len(chr) > 2:
                continue
            if gtf_entry.molecule == 'transcript':
                keep = gtf_entry.annotation.find('transcript_biotype "protein_coding"') != -1
            if gtf_entry.molecule == 'exon' and keep:
                out_line = f'{chr}\t{gtf_entry.start_pos - padding}\t{gtf_entry.end_pos + padding}\n'
                print(out_line, file=out_fh)


def write_to_bed(out_bed_fp: str, bed2add_fp: str, fourth_col: str):
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


def compare_clinvar(new_clinvar: dict[str, str], old_clinvar: dict[str, str], final_bed_fp: str) -> tuple[list[str], list[str]]:
    
    new_bed_fp = 'clinvar_new.bed'
    old_bed_fp = 'clinvar_old.bed'
    clinvar_in_common = {}
    padding = 5
    with open(new_bed_fp, 'w') as new_fh, open(old_bed_fp, 'w') as old_fh:
        out_line = [

        ]


def get_intersect(clinvar: str, bed: str) -> list[str]:
    not_in_bed_cmd = ['bedtools', 'intersect', '-a', clinvar, '-b', bed, '-v']
    # FIXME: How to gather and return results
    subprocess.call(not_in_bed_cmd)
    return []


def get_clinvar_info(info: dict[str,str], var: str) -> str:
    pass

def log_clinvar_bed_and_info(clinvar: str, new_or_old: str, clinvar_log: str):
    pass

if __name__ == "__main__":
    main()


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