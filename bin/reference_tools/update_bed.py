#!/usr/bin/env python3

import argparse
import subprocess
import os
from pathlib import Path
import re

# import pdb

"""
CLNSIG:      Clinical significance of variant
CLNSIGINCL:  Clinical significance for a haplotype or genotype that includes
                this variant. Reported as pairs (VariantID:ClinSig).
                It can for instance be benign by itself, but when found together with a separate 
                variant on the same location, or different location, it could be Pathogenic.
CLNSIGCONF:  Reviewer certainty
CLNDN:       ClinVar preferred disease name
CLNACC:      ClinVar accession (deprecated? https://github.com/Clinical-Genomics/scout/issues/695)
"""

CLNSIG = "CLNSIG"
CLNACC = "CLNACC"
CLNDN = "CLNDN"
CLNSIGINCL = "CLNSIGINCL"
CLNSIGCONF = "CLNSIGCONF"

VARIANT_PAD = 5
EXON_PAD = 20

# Testing steps
# 1. Are the base exons files identical?
# 2.


def main(
    clinvar_vcf_path: Path,
    clinvar_vcf_old_path: Path,
    clinvardate: str,
    out_dir: Path,
    release: str,
    skip_download: bool,
    incl_bed: list[str],
):
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    final_bed_path = ensure_new_empty(
        f"{out_dir}/exons_{release}padded{EXON_PAD}bp_clinvar-{clinvardate}padded{VARIANT_PAD}bp.bed"
    )

    ensembl_bed_path = Path(f"{out_dir}/exons_hg38_{release}.bed")

    print("Write initial base (ENSEMBL exons)")
    write_ensembl(ensembl_bed_path, release, skip_download, EXON_PAD)

    if len(incl_bed) > 0:
        for bed_fp in incl_bed:
            print(f"Additional included BED file: {bed_fp}")
            suffix = bed_fp
            append_to_bed(final_bed_path, Path(bed_fp), suffix)
    else:
        print("No extra BED files to include")

    append_to_bed(final_bed_path, ensembl_bed_path, f"EXONS-{release}")

    (clinvar_new, new_benign) = read_clinvar(clinvar_vcf_path)
    (clinvar_old, _old_benign) = read_clinvar(clinvar_vcf_old_path)

    print(f"Number new: {len(clinvar_new)}")
    print(f"Number old: {len(clinvar_old)}")
    print(f"Number new benign: {len(new_benign)}")
    print(f"Number old benign: {len(_old_benign)}")

    clinvar_all: dict[str, Variant] = {**clinvar_new, **clinvar_old}

    (new_bed_path, new_clinvar_bed_rows, old_clinvar_bed_rows) = compare_clinvar(
        clinvar_new, clinvar_old, final_bed_path, VARIANT_PAD, out_dir
    )

    clinvar_log_path = ensure_new_empty(f"{out_dir}/clinvar_{clinvardate}.log")
    log_changes(
        clinvar_log_path,
        clinvar_all,
        new_clinvar_bed_rows,
        old_clinvar_bed_rows,
        new_benign,
    )

    append_to_bed(final_bed_path, new_bed_path, ".")

    sort_merge_output(str(final_bed_path))


def ensure_new_empty(filepath: str) -> Path:
    path = Path(filepath)
    if path.exists():
        print(f"Removing old file: {filepath}")
        path.unlink()
    path.write_text("")
    return path


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", required=True)
    parser.add_argument("--new", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--clinvardate", required=True)

    parser.add_argument("--out_dir", required=True)

    parser.add_argument("--incl_bed", nargs="*")
    parser.add_argument("--skip_download", action="store_true")
    args = parser.parse_args()
    return args


class Variant:

    @staticmethod
    def from_line(vcf_line: str) -> "Variant":

        vcf_line = vcf_line.rstrip()
        fields = vcf_line.split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        info_str = fields[7]

        info: dict[str, str] = {}
        for info_field in info_str.split(";"):
            (key, val) = info_field.split("=")
            info[key] = val

        return Variant(chrom, pos, ref, alt, info)

    def __init__(
        self, chrom: str, pos: int, ref: str, alt: str, info: dict[str, str] = {}
    ):

        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

        self.reason: str | None = None

        self.clndn = "Undefined"
        if self.info.get("CLNDN") is not None:
            self.clndn = self.info["CLNDN"]

        self.key = "FIXME"

    def get_key(self) -> str:
        return f"{self.chrom}:{self.pos}_{self.ref}_{self.alt}"

    def get_padded_bed(self, padding: int, annot: str | None = None) -> str:
        out_fields = [self.chrom, str(self.pos - padding), str(self.pos + padding)]
        if annot is not None:
            out_fields.append(annot)
        return "\t".join(out_fields)

    def get_tmp_bed_str_fixme(self, padding: int) -> str:
        return f"{self.chrom}\t{self.pos - padding}\t{self.pos + padding}\t{self.get_tmp_bed_annot_fixme()}"

    def get_tmp_bed_annot_fixme(self):
        reason = self.reason
        clnacc = self.info[CLNACC] if self.info.get(CLNACC) is not None else "Undefined"
        clndn = self.info[CLNDN] if self.info.get(CLNDN) is not None else "Undefined"
        return (
            f"{self.chrom}:{self.pos}_{self.ref}_{self.alt}~{reason}~{clnacc}~{clndn}"
        )

    def get_clnsig(self) -> str:
        clnsig = self.info.get(CLNSIG)
        if clnsig != None:
            return clnsig
        else:
            return "MISSING"


def read_clinvar(vcf_path: Path) -> tuple[dict[str, Variant], dict[str, Variant]]:
    """
    FIXME: Steps

    -
    - Mitochondria variants are not included
    """

    clinvar_variants: dict[str, Variant] = {}
    benign_clinvar_variants: dict[str, Variant] = {}

    with vcf_path.open("r") as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            variant = Variant.from_line(line)

            sig = ""
            confidence = ""
            haplo = ""

            info = variant.info

            if info.get(CLNSIG) is not None:
                sig = info[CLNSIG]
            if info.get(CLNSIGINCL) is not None:
                haplo = info[CLNSIGINCL]
            if info.get(CLNSIGCONF) is not None:
                confidence = info[CLNSIGCONF]

            # FIXME: Annotate why
            # Skipping entries where no significance is assigned while
            # the haplotype is
            if sig == "" and haplo != "":
                continue

            keep = False
            reason = None

            if "Pathogenic" in sig or "Likely_pathogenic" in sig:
                keep = True
                reason = sig
            elif "Conflicting_interpretations_of_pathogenicity" in sig:
                if "Pathogenic" in confidence or "Likely_pathogenic" in confidence:
                    keep = True
                    reason = confidence
            elif "athogenic" in haplo:
                keep = True
                reason = haplo

            # FIXME: Multiple variants have the same key, so the last will be used
            # Is this something we should think about?
            key = f"{variant.chrom}:{variant.pos}_{variant.ref}_{variant.alt}"
            if keep:
                # FIXME: Unsure if relevant to keep these for benign
                # This is how it was done in the pl script
                if variant.chrom == "MT":
                    continue
                variant.reason = reason
                clinvar_variants[key] = variant
            else:
                variant.reason = reason
                benign_clinvar_variants[key] = variant

    return (clinvar_variants, benign_clinvar_variants)


def write_ensembl(out_path: Path, release: str, skip_download: bool, exon_padding: int):

    ensembl_tmp_gz_fp = Path("tmp.gtf.gz")
    ensembl_tmp_fp = Path("tmp.gtf")

    if not skip_download:
        gtf_request = f"https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz"
        download_gtf_cmd = ["wget", gtf_request, "-O", ensembl_tmp_gz_fp]
        subprocess.run(download_gtf_cmd, check=True)
        gunzip_cmd = ["gunzip", ensembl_tmp_gz_fp]
        subprocess.run(gunzip_cmd, check=True)
    else:
        if not Path(ensembl_tmp_fp).exists():
            raise ValueError(
                f"To skip download, the ENSEMBL GTF needs to be present at the location: {ensembl_tmp_fp}"
            )

    with ensembl_tmp_fp.open("r") as gtf_fh, out_path.open("w") as out_fh:
        keep = False
        for line in gtf_fh:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom = fields[0]
            molecule = fields[2]
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            annotation = fields[8]

            # Filter out non chr chromosomes
            # Added "replace" to allow "chr" based reference
            chrom = re.sub("^chr", "", chrom)
            if len(chrom) > 2:
                continue
            if molecule == "transcript":
                keep = annotation.find('transcript_biotype "protein_coding"') != -1
            if molecule == "exon" and keep:
                out_line = (
                    f"{chrom}\t{start_pos - exon_padding}\t{end_pos + exon_padding}"
                )
                print(out_line, file=out_fh)


def append_to_bed(out_bed_fp: Path, bed2add_fp: Path, bed_annot_default: str):
    with bed2add_fp.open("r") as bed2add_fh, out_bed_fp.open("a") as out_fh:
        for line in bed2add_fh:
            line = line.rstrip()
            fields = line.split("\t")
            if len(fields) > 3:
                print(line, file=out_fh)
            else:
                out_fields = fields + [bed_annot_default]
                print("\t".join(out_fields), file=out_fh)


def sort_merge_output(bed_fp: str):
    tmp_bed_fp = "tmp.sort.bed"
    sort_cmd = f"bedtools sort -i {bed_fp} > {tmp_bed_fp}"
    # Annotation column is concatenated together
    # "distinct" means the same value won't be reused multiple times
    merge_cmd = f"bedtools merge -i {tmp_bed_fp} -c 4 -o collapse > {bed_fp}"
    # merge_cmd = f"bedtools merge -i {tmp_bed_fp} -c 4 -o distinct > {bed_fp}"
    subprocess.call(sort_cmd, shell=True)
    subprocess.call(merge_cmd, shell=True)
    os.remove(tmp_bed_fp)


def compare_clinvar(
    new_clinvar: dict[str, Variant],
    old_clinvar: dict[str, Variant],
    final_bed_path: Path,
    variant_padding: int,
    tmp_dir: Path,
) -> tuple[Path, list[str], list[str]]:
    """
    1. Generate bed files for ClinVar variants (pos +/- padding)
    2. Add info about pathogenicity as the fourth column
    Returns lists of chrom:pos_ref_alt keys
    """

    new_clinvar_keys = set(new_clinvar.keys())
    old_clinvar_keys = set(old_clinvar.keys())

    def write_tmp_bed(
        path: Path, variant_dict: dict[str, Variant], filter_keys: set[str]
    ):
        subset_dict = {
            key: val
            for (key, val) in variant_dict.items()
            if key in filter_keys or len(filter_keys) == 0
        }
        bed_text = "\n".join(
            [
                variant.get_padded_bed(variant_padding, annot=var_key)
                for (var_key, variant) in subset_dict.items()
            ]
        )
        path.write_text(bed_text)

    new_clinvar_tmp_bed = tmp_dir / "clinvar_new_python.bed"
    write_tmp_bed(new_clinvar_tmp_bed, new_clinvar, set())
    old_clinvar_tmp_bed = tmp_dir / "clinvar_old_python.bed"
    clinvar_old_removed = old_clinvar_keys.difference(new_clinvar_keys)
    write_tmp_bed(old_clinvar_tmp_bed, old_clinvar, clinvar_old_removed)

    clinvar_in_common = new_clinvar_keys.intersection(old_clinvar_keys)
    clinvar_new_added = new_clinvar_keys.difference(old_clinvar_keys)

    new_bed_rows = get_bed_intersect(str(new_clinvar_tmp_bed), str(final_bed_path))
    old_bed_rows = get_bed_intersect(str(old_clinvar_tmp_bed), str(final_bed_path))
    # new_clinvar_tmp_bed.unlink()
    # old_clinvar_tmp_bed.unlink()

    # new_bed_keys = [row.split("\t")[] for row in new_bed_rows]

    print(f"Clinvar in common between versions: {len(clinvar_in_common)}")
    print(f"Added new (unique targets): {len(clinvar_new_added)} ({len(new_bed_rows)})")
    print(
        f"Removed old (unique targets): {len(clinvar_old_removed)} ({len(old_bed_rows)})"
    )

    return (new_clinvar_tmp_bed, new_bed_rows, old_bed_rows)


def get_bed_intersect(left_bed: str, right_bed: str) -> list[str]:
    # not_in_bed_cmd = ["bedtools", "intersect", "-h"]
    not_in_bed_cmd = [
        "bedtools",
        "intersect",
        "-a",
        left_bed,
        "-b",
        right_bed,
        "-v",
    ]
    print(f"Running command: {' '.join(not_in_bed_cmd)}")
    result = subprocess.run(not_in_bed_cmd, capture_output=True, text=True, check=True)
    intersected_regions = result.stdout.strip().splitlines()
    # FIXME: This differs from the Perl script, and gives subtly differences in unique targets
    # In this case 1270 vs 1274 (i.e. four are duplicated)
    return intersected_regions
    # return set(intersected_regions)


def log_changes(
    log_path: Path,
    clinvar_all: dict[str, Variant],
    added_keys: list[str],
    removed_keys: list[str],
    new_benign: dict[str, Variant],
):
    """
    Log variants added and removed between clinvar versions
    Check if part of latest ClinVar with benign status if removed
    """
    with log_path.open("w") as log_fh:
        for clin_key, clin_var in clinvar_all.items():
            if clin_key in added_keys:
                clnsig = "MISSING"
                if clin_key in new_benign:
                    clnsig = new_benign[clin_key].get_clnsig()
                print(
                    f"REMOVED: {clin_var.pos} reason[2] reason[1] => {clnsig}",
                    file=log_fh,
                )

            elif clin_key in removed_keys:
                print(f"ADDED: {clin_var.pos} reason[2] reason[1]", file=log_fh)


# def append_clinvar_to_bed(bed_fp: str, clinvar_new: dict[str, Variant]):
#     with open(bed_fp, "w") as bed_fh:
#         for clin_var in clinvar_new.values():
#             print(
#                 f"{clin_var.chrom}\t{clin_var.pos}\t{clin_var.pos}\tCLINVAR-{clin_var.get_clnsig()}",
#                 file=bed_fh,
#             )


# def write_clinvar_bed(
#     clinvar_out_path: Path, bed_rows: list[str], variant_padding: int
# ):
#     with clinvar_out_path.open("a") as out_fh:
#         for variant in clinvar_new:
#             print(variant.get_tmp_bed_str_fixme(variant_padding), file=out_fh)


if __name__ == "__main__":
    args = parse_arguments()
    incl_bed: list[str] | None = args.incl_bed
    main(
        clinvar_vcf_path=Path(args.new),
        clinvar_vcf_old_path=Path(args.old),
        clinvardate=args.clinvardate,
        out_dir=Path(args.out_dir),
        release=args.release,
        skip_download=args.skip_download,
        incl_bed=incl_bed if incl_bed is not None else [],
    )
