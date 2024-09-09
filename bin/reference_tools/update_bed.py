#!/usr/bin/env python3

"""
This script calculates the so-called "intersect bed", which contains the intervals relevant for annotation
and scoring in the wgs pipeline.

It is based on ENSEMBL exons which are downloaded and padded
It is extended with ClinVar variants, also padded
It compares with an older ClinVar version and provides a summary of what has been added / removed
Further custom bed files can be provided as well with

The final output is a sorted and merged bed file with regions based on the steps above
"""

import argparse
import subprocess
from pathlib import Path
import re
import gzip
import requests
import logging

LOG = logging.getLogger(__name__)

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
UNDEFINED = "Undefined"

VARIANT_PAD = 5
EXON_PAD = 20


def get_gtf_url(release: str) -> str:
    return f"https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz"


def main(
    clinvar_vcf_path: Path,
    clinvar_vcf_old_path: Path,
    clinvardate: str,
    out_dir: Path,
    release: str,
    skip_download: bool,
    incl_bed: list[str],
    keep_tmp: bool,
):
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    final_bed_path = ensure_new_empty(
        f"{out_dir}/exons_{release}padded{EXON_PAD}bp_clinvar-{clinvardate}padded{VARIANT_PAD}bp.bed"
    )

    ensembl_bed_path = Path(f"{out_dir}/exons_hg38_{release}.bed")

    LOG.info("Write initial base (ENSEMBL exons)")
    write_ensembl_bed(out_dir, ensembl_bed_path, release, skip_download, EXON_PAD)

    if len(incl_bed) > 0:
        for bed_fp in incl_bed:
            LOG.info(f"Additional included BED file: {bed_fp}")
            suffix = bed_fp
            append_to_bed(final_bed_path, Path(bed_fp), suffix)
    else:
        LOG.info("No extra BED files to include")

    append_to_bed(final_bed_path, ensembl_bed_path, f"EXONS-{release}")

    (clinvar_new, new_benign) = read_clinvar(clinvar_vcf_path)
    (clinvar_old, _old_benign) = read_clinvar(clinvar_vcf_old_path)

    (new_to_add_clinvar_bed_rows, old_to_remove_clinvar_bed_rows) = compare_clinvar(
        clinvar_new, clinvar_old, final_bed_path, VARIANT_PAD, out_dir, keep_tmp
    )

    new_to_remove_keys = [bed_row[4] for bed_row in new_to_add_clinvar_bed_rows]
    old_to_remove_keys = [bed_row[4] for bed_row in old_to_remove_clinvar_bed_rows]

    clinvar_log_path = ensure_new_empty(f"{out_dir}/clinvar_{clinvardate}.log")
    log_changes(
        clinvar_log_path,
        clinvar_old,
        clinvar_new,
        set(new_to_remove_keys),
        set(old_to_remove_keys),
        new_benign,
    )

    new_clinvar_to_add_path = Path(f"{out_dir}/new_clinvar_to_add.bed")
    new_clinvar_to_add_path.write_text(
        "\n".join(["\t".join(row[0:3]) for row in new_to_add_clinvar_bed_rows]) + "\n"
    )

    append_to_bed(final_bed_path, new_clinvar_to_add_path, ".")
    sort_merge_output(out_dir, final_bed_path, keep_tmp)


def ensure_new_empty(filepath: str) -> Path:
    path = Path(filepath)
    if path.exists():
        path.unlink()
    path.touch()
    return path


class ClinVarVariant:

    @staticmethod
    def from_line(vcf_line: str) -> "ClinVarVariant":

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

        return ClinVarVariant(chrom, pos, ref, alt, info)

    def __init__(
        self, chrom: str, pos: int, ref: str, alt: str, info: dict[str, str] = {}
    ):

        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

        self.reason = UNDEFINED

        self.CLNDN = self.info[CLNDN] if self.info.get(CLNDN) is not None else UNDEFINED
        self.CLNSIG = (
            self.info[CLNSIG] if self.info.get(CLNSIG) is not None else UNDEFINED
        )

    @property
    def reason(self):
        return self._reason

    @reason.setter
    def reason(self, value: str):
        self._reason = value

    def generate_key(self) -> str:
        return f"{self.chrom}:{self.pos}_{self.ref}_{self.alt}"

    def get_padded_bed(self, padding: int) -> str:
        out_fields = [
            self.chrom,
            str(self.pos - padding),
            str(self.pos + padding),
            f"CLINVAR-{self.reason}",
        ]
        return "\t".join(out_fields)


def read_clinvar(
    vcf_path: Path,
) -> tuple[dict[str, ClinVarVariant], dict[str, ClinVarVariant]]:
    """
    1. Given a ClinVar VCF file
    2. Retrieve ClinVar significance (CLNSIG), significance when co-occuring (CLNSIGINCL)
       and reviewer confidence (CLNSIGCONF)
    """

    clinvar_variants: dict[str, ClinVarVariant] = {}
    benign_clinvar_variants: dict[str, ClinVarVariant] = {}

    with vcf_path.open("r") as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            variant = ClinVarVariant.from_line(line)

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
            # FIXME: Investigate
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

                if reason is not None:
                    variant.reason = reason
                clinvar_variants[key] = variant
            else:
                if reason is not None:
                    variant.reason = reason
                benign_clinvar_variants[key] = variant

    return (clinvar_variants, benign_clinvar_variants)


def write_ensembl_bed(
    out_dir: Path, out_bed: Path, release: str, skip_download: bool, exon_padding: int
):
    """Print padded exons for protein coding transcripts"""

    ensembl_tmp_path = out_dir / "tmp.gtf.gz"

    if not skip_download:
        gtf_url = get_gtf_url(release)
        requests.get(gtf_url, stream=True)

        # FIXME: Cleanup
        # download_gtf_cmd = ["wget", gtf_request, "-O", ensembl_tmp_fp]
        # subprocess.run(download_gtf_cmd, check=True)
        # gunzip_cmd = ["gunzip", ensembl_tmp_gz_fp]
        # subprocess.run(gunzip_cmd, check=True)
    else:
        if not Path(ensembl_tmp_path).exists():
            raise ValueError(
                f"To skip download, the ENSEMBL GTF needs to be present at the location: {ensembl_tmp_path}"
            )

    def is_valid_chromosome(chr: str) -> bool:
        chr_stripped_chr = re.sub("^chr", "", chr)
        return len(chr_stripped_chr) > 2

    with gzip.open(str(ensembl_tmp_path), "rt") as gtf_fh, out_bed.open("w") as out_fh:
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
            if is_valid_chromosome(chrom):
                continue
            if molecule == "transcript":
                keep = annotation.find('transcript_biotype "protein_coding"') != -1
            if molecule == "exon" and keep:
                out_line = (
                    f"{chrom}\t{start_pos - exon_padding}\t{end_pos + exon_padding}"
                )
                print(out_line, file=out_fh)


def append_to_bed(out_bed_fp: Path, bed_to_add_path: Path, bed_annot_default: str):
    with bed_to_add_path.open("r") as bed2add_fh, out_bed_fp.open("a") as out_fh:
        for line in bed2add_fh:
            line = line.rstrip()
            fields = line.split("\t")
            if len(fields) > 3:
                print(line, file=out_fh)
            else:
                out_fields = fields + [bed_annot_default]
                print("\t".join(out_fields), file=out_fh)


def sort_merge_output(out_dir: Path, bed_fp: Path, keep_tmp: bool):
    tmp_bed_path = Path(out_dir / "tmp.sort.bed")
    sort_cmd = f"bedtools sort -i {str(bed_fp)} > {str(tmp_bed_path)}"
    # Annotation column is concatenated together
    # "distinct" means the same value won't be reused multiple times
    merge_cmd = (
        f"bedtools merge -i {str(tmp_bed_path)} -c 4 -o distinct > {str(bed_fp)}"
    )
    subprocess.call(sort_cmd, shell=True)
    subprocess.call(merge_cmd, shell=True)
    if not keep_tmp:
        tmp_bed_path.unlink()


def compare_clinvar(
    new_clinvar: dict[str, ClinVarVariant],
    old_clinvar: dict[str, ClinVarVariant],
    final_bed_path: Path,
    variant_padding: int,
    out_dir: Path,
    keep_tmp: bool,
) -> tuple[list[list[str]], list[list[str]]]:
    """
    1. Check what variants are new / kept / removed among keys
    2. Generate bed files for padded ClinVar variants to run "bedtools intersect"
    3. Keep info about pathogenicity as the fourth column and variant key (chr:pos_ref_alt) in fifth
    4. Return bed field lists
    """

    new_clinvar_keys = set(new_clinvar.keys())
    old_clinvar_keys = set(old_clinvar.keys())

    def write_tmp_bed(
        path: Path, variant_dict: dict[str, ClinVarVariant], filter_keys: set[str]
    ):
        subset_dict = {
            key: val
            for (key, val) in variant_dict.items()
            if key in filter_keys or len(filter_keys) == 0
        }
        bed_text = (
            "\n".join(
                [
                    variant.get_padded_bed(variant_padding) + f"\t{var_key}"
                    for (var_key, variant) in subset_dict.items()
                ]
            )
            + "\n"
        )
        path.write_text(bed_text)

    new_clinvar_tmp_bed = out_dir / "clinvar_new_python.bed"
    write_tmp_bed(new_clinvar_tmp_bed, new_clinvar, set())
    old_clinvar_tmp_bed = out_dir / "clinvar_old_python.bed"
    clinvar_old_removed = old_clinvar_keys.difference(new_clinvar_keys)
    write_tmp_bed(old_clinvar_tmp_bed, old_clinvar, clinvar_old_removed)

    clinvar_in_common = new_clinvar_keys.intersection(old_clinvar_keys)
    clinvar_new_added = new_clinvar_keys.difference(old_clinvar_keys)

    new_to_add_bed_rows = get_bed_intersect(
        str(new_clinvar_tmp_bed), str(final_bed_path)
    )
    old_to_remove_bed_rows = get_bed_intersect(
        str(old_clinvar_tmp_bed), str(final_bed_path)
    )
    if not keep_tmp:
        new_clinvar_tmp_bed.unlink()
        old_clinvar_tmp_bed.unlink()

    LOG.info(f"Clinvar in common between versions: {len(clinvar_in_common)}")
    LOG.info(
        f"Added new (unique targets): {len(clinvar_new_added)} ({len(new_to_add_bed_rows)})"
    )
    LOG.info(
        f"Removed old (unique targets): {len(clinvar_old_removed)} ({len(old_to_remove_bed_rows)})"
    )

    return (new_to_add_bed_rows, old_to_remove_bed_rows)


def get_bed_intersect(left_bed: str, right_bed: str) -> list[list[str]]:
    not_in_bed_cmd = [
        "bedtools",
        "intersect",
        "-a",
        left_bed,
        "-b",
        right_bed,
        "-v",
    ]
    result = subprocess.run(not_in_bed_cmd, capture_output=True, text=True, check=True)
    intersected_regions = result.stdout.strip().splitlines()
    # FIXME: This differs from the Perl script, and gives subtly differences in unique targets
    # In this case 1270 vs 1274 (i.e. four are duplicated)
    return [row.split("\t") for row in intersected_regions]
    # return set(intersected_regions)


def log_changes(
    log_path: Path,
    clinvar_old: dict[str, ClinVarVariant],
    clinvar_new: dict[str, ClinVarVariant],
    new_clinvar_to_add_keys: set[str],
    old_clinvar_to_remove_keys: set[str],
    new_benign: dict[str, ClinVarVariant],
):
    """
    Log variants added and removed between clinvar versions
    Check if part of latest ClinVar with benign status if removed
    """
    with log_path.open("w") as log_fh:

        for new_key in new_clinvar_to_add_keys:

            new_variant = clinvar_new.get(new_key)
            if new_variant is None:
                raise ValueError(f"New variant should be present for key {new_key}")

            print(
                f"ADDED: {new_variant.chrom}:{new_variant.pos}:{new_variant.CLNDN}:{new_variant.CLNSIG}",
                file=log_fh,
            )

        for old_key in old_clinvar_to_remove_keys:
            reason_for_removal = "MISSING"
            if old_key in new_benign:
                reason_for_removal = new_benign[old_key].CLNSIG
            old_variant = clinvar_old.get(old_key)
            if old_variant is None:
                raise ValueError(f"New variant should be present for key {old_key}")
            print(
                f"REMOVED: {old_variant.chrom}:{old_variant.pos}:{old_variant.CLNDN}:{old_variant.CLNSIG} => {reason_for_removal}",
                file=log_fh,
            )


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", required=True)
    parser.add_argument("--new", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--clinvardate", required=True)

    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--incl_bed", nargs="*")
    parser.add_argument("--skip_download", action="store_true")
    parser.add_argument(
        "--keep_tmp", action="store_true", help="Don't remove tmp files for debugging"
    )
    args = parser.parse_args()
    return args


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
        keep_tmp=args.keep_tmp,
    )
