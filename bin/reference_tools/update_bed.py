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


def main(
    clinvar_vcf_path: Path,
    clinvar_vcf_old_path: Path,
    clinvardate: str,
    out_dir: Path,
    release: str,
    skip_download: bool,
    incl_bed: list[str],
):

    out_dir = Path(out_dir)

    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    ensembl_bed_path = Path(f"{out_dir}/exons_hg38_{release}.bed")

    variant_pad = 5
    exon_pad = 20
    print("Write initial base")
    write_ensembl(ensembl_bed_path, release, skip_download, exon_pad)

    final_bed_path = ensure_new_path(
        f"{out_dir}/exons_{release}padded{exon_pad}bp_clinvar-{clinvardate}padded{variant_pad}bp.bed"
    )
    final_bed_path.write_text("")

    append_to_bed(final_bed_path, ensembl_bed_path, f"EXONS-{release}")

    if len(incl_bed) > 0:
        for bed_fp in incl_bed:
            print(f"Additional included BED file: {bed_fp}")
            suffix = bed_fp
            append_to_bed(final_bed_path, Path(bed_fp), suffix)
    else:
        print("No extra BED files to include")

    (clinvar_new, new_benign) = read_clinvar(clinvar_vcf_path)
    (clinvar_old, _old_benign) = read_clinvar(clinvar_vcf_old_path)

    print(f"Number new: {len(clinvar_new)}")
    print(f"Number old: {len(clinvar_old)}")
    print(f"Number new benign: {len(new_benign)}")
    print(f"Number old benign: {len(_old_benign)}")

    clinvar_all: dict[str, Variant] = {**clinvar_new, **clinvar_old}

    (new_to_add, old_to_remove) = compare_clinvar(
        clinvar_new, clinvar_old, final_bed_path, variant_pad, out_dir
    )

    clinvar_log_path = ensure_new_path(f"{out_dir}/clinvar_{clinvardate}.log")
    log_changes(
        str(clinvar_log_path), clinvar_all, new_to_add, old_to_remove, new_benign
    )

    append_new_clinvar(final_bed_path, clinvar_new, variant_pad)

    sort_merge_output(str(final_bed_path))


def ensure_new_path(filepath: str) -> Path:
    path = Path(filepath)
    if path.exists():
        print(f"Removing old file: {filepath}")
        path.unlink()
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

    def get_bed_str(self, padding: int) -> str:
        return f"{self.chrom}\t{self.pos - padding}\t{self.pos + padding}\t{self.get_bed_annot()}"

    def get_bed_annot(self):
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
    merge_cmd = f"bedtools merge -i {tmp_bed_fp} -c 4 -o distinct > {bed_fp}"
    subprocess.call(sort_cmd, shell=True)
    subprocess.call(merge_cmd, shell=True)
    os.remove(tmp_bed_fp)


def compare_clinvar(
    new_clinvar: dict[str, Variant],
    old_clinvar: dict[str, Variant],
    final_bed_path: Path,
    padding: int,
    tmp_dir: Path,
) -> tuple[set[str], set[str]]:

    new_bed_path = tmp_dir / "clinvar_new_python.bed"
    old_bed_path = tmp_dir / "clinvar_old_python.bed"
    clinvar_in_common: set[str] = set()
    with open(str(new_bed_path), "w") as new_fh:
        for key, variant in new_clinvar.items():
            if old_clinvar.get(key) is not None:
                clinvar_in_common.add(key)

            # FIXME: Look over this
            clinvar_info = variant.get_bed_annot()
            out_line = f"{variant.chrom}\t{variant.pos - padding}\t{variant.pos + padding}\t{clinvar_info}"
            print(out_line, file=new_fh)

    clinvar_new_added: set[str] = set()
    for new_key in new_clinvar:
        if new_key not in old_clinvar:
            clinvar_new_added.add(new_key)

    with open(old_bed_path, "w") as old_fh:
        clinvar_old_removed: set[str] = set()
        for old_key in old_clinvar:
            if old_key not in new_clinvar:
                old_var = old_clinvar[old_key]
                clinvar_old_removed.add(old_key)
                # FIXME: Write this to old BED
                print(old_var.get_bed_str(padding), file=old_fh)

    new_to_add = get_bed_intersect(str(new_bed_path), str(final_bed_path))
    old_to_remove = get_bed_intersect(str(old_bed_path), str(final_bed_path))
    new_bed_path.unlink()
    old_bed_path.unlink()

    print(f"Clinvar in common between versions: {len(clinvar_in_common)}")
    print(f"Added new (unique targets): {len(clinvar_new_added)} ({len(new_to_add)})")
    print(
        f"Removed old (unique targets): {len(clinvar_old_removed)} ({len(old_to_remove)})"
    )

    return (new_to_add, old_to_remove)


def get_bed_intersect(left_bed: str, right_bed: str) -> set[str]:
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
    return set(intersected_regions)


def log_changes(
    log_fp: str,
    clinvar_all: dict[str, Variant],
    added_keys: set[str],
    removed_keys: set[str],
    new_benign: dict[str, Variant],
):
    """
    Log variants added and removed between clinvar versions
    Check if part of latest ClinVar with benign status if removed
    """
    with open(log_fp, "w") as log_fh:
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


def append_new_clinvar(
    final_bed_path: Path, clinvar_new: dict[str, Variant], variant_padding: int
):
    with final_bed_path.open("a") as out_fh:
        for variant in clinvar_new.values():
            print(variant.get_bed_str(variant_padding), file=out_fh)


if __name__ == "__main__":
    args = parse_arguments()
    main(
        clinvar_vcf_path=Path(args.new),
        clinvar_vcf_old_path=Path(args.old),
        clinvardate=args.clinvardate,
        out_dir=Path(args.out_dir),
        release=args.release,
        skip_download=args.skip_download,
        incl_bed=args.incl_bed,
    )
