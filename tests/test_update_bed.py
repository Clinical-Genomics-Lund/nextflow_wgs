from unittest.mock import patch
import sys
from pathlib import Path

# sys.path.insert(0, str(Path(__file__).resolve().parent / 'bin'))

from src.update_bed import write_base, append_to_bed, compare_clinvar, read_clinvar, Variant

mock_gtf_content = '''\
# Example GTF file
1\tENSEMBL\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gene1"; transcript_biotype "protein_coding";
1\tENSEMBL\texon\t1100\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
2\tENSEMBL\ttranscript\t3000\t4000\t.\t+\t.\tgene_id "gene2"; transcript_biotype "non_coding";
2\tENSEMBL\texon\t3100\t3200\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2";
3\tENSEMBL\ttranscript\t4000\t5000\t.\t+\t.\tgene_id "gene3"; transcript_biotype "protein_coding";
3\tENSEMBL\texon\t4400\t4500\t.\t+\t.\tgene_id "gene3"; transcript_id "transcript3";
'''

mock_chr_gtf_content = '''\
# Example GTF file
chr1\tENSEMBL\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gene1"; transcript_biotype "protein_coding";
chr1\tENSEMBL\texon\t1100\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr2\tENSEMBL\ttranscript\t3000\t4000\t.\t+\t.\tgene_id "gene2"; transcript_biotype "non_coding";
chr2\tENSEMBL\texon\t3100\t3200\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2";
chr3\tENSEMBL\ttranscript\t4000\t5000\t.\t+\t.\tgene_id "gene3"; transcript_biotype "protein_coding";
chr3\tENSEMBL\texon\t4400\t4500\t.\t+\t.\tgene_id "gene3"; transcript_id "transcript3";
'''



def get_file_string(rows: list[list[str]]) -> str:
    """
    Convert 
    [['chr1', '1', '2'], ['chr2', '2', '3']]
    to
    'chr1\t1\t2\nchr2\t2\t3'
    """
    row_strs = []
    for row in rows:
        row_str = '\t'.join(row)
        row_strs.append(row_str)
        
    return '\n'.join(row_strs) + '\n'


def test_write_base(tmp_path: Path):

    release = "108"
    skip_download = True
    padding = 20
    out_path = tmp_path / 'output.bed'

    # Non-chr based reference
    ensembl_nonchr_path = tmp_path / 'test_nonchr.gtf'
    ensembl_nonchr_path.write_text(mock_gtf_content)
    write_base(str(ensembl_nonchr_path), str(out_path), release, skip_download, padding)
    output_content = out_path.read_text()
    row_fields = [
        ['1', '1080', '1220'],
        ['3', '4380', '4520']
    ]
    expected_output = get_file_string(row_fields)
    assert output_content == expected_output, f'Output was: {output_content}'

    # Chr based reference
    ensembl_chr_path = tmp_path / 'test_chr.gtf'
    ensembl_chr_path.write_text(mock_chr_gtf_content)
    write_base(str(ensembl_chr_path), str(out_path), release, skip_download, padding)
    output_content = out_path.read_text()
    row_fields = [
        ['chr1', '1080', '1220'],
        ['chr3', '4380', '4520']
    ]
    expected_output = get_file_string(row_fields)
    assert output_content == expected_output, f'Output was: {output_content}'


def test_append_to_bed(tmp_path: Path):
    
    base_bed_content = '''\
chr1\t1\t5\told_annot
chr1\t10\t15\told_annot2
'''

    mock_bed_content = '''\
chr1\t1000\t2000
chr1\t3000\t4000\tgene1
chr2\t5000\t6000
'''

    expected_output = '''\
chr1\t1\t5\told_annot
chr1\t10\t15\told_annot2
chr1\t1000\t2000\tdefault_annot
chr1\t3000\t4000\tgene1
chr2\t5000\t6000\tdefault_annot
'''

    out_bed_path = tmp_path / 'out_bed.bed'
    bed2add_path = tmp_path / 'bed2add.bed'

    out_bed_path.write_text(base_bed_content)
    bed2add_path.write_text(mock_bed_content)

    fourth_col_default = 'default_annot'
    append_to_bed(str(out_bed_path), str(bed2add_path), fourth_col_default)
    result = out_bed_path.read_text()

    assert result == expected_output, f'Output was: {result}'


# FIXME: Investigate the CLNSIG, CLNSIGCONF, CLNSIGINCL combinations in real ClinVar data
def test_read_clinvar(tmp_path: Path):

    mock_vcf_content = """\
##fileformat=VCFv4.2
##source=ClinVar
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t12345\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNSIGCONF=Conflicting_interpretations_of_pathogenicity;CLNSIGINCL=athogenic
chr1\t54321\t.\tC\tT\t.\t.\tCLNSIG=Benign
chr2\t12345\t.\tG\tA\t.\t.\tCLNSIG=Likely_pathogenic
"""

    clinvar_vcf_path = tmp_path / 'tmp.vcf'
    clinvar_vcf_path.write_text(mock_vcf_content)

    clinvar_variants, benign_clinvar_variants = read_clinvar(str(clinvar_vcf_path))

    expected_clinvar_reasons = {
        'chr1:12345_A_G': 'Pathogenic',
        'chr2:12345_G_A': 'Likely_pathogenic'
    }
    expected_benign_clinvar_reasons = {
        'chr1:54321_C_T': None
    }

    for key, variant in clinvar_variants.items():
        assert key in expected_clinvar_reasons
        assert expected_clinvar_reasons[key] == variant.reason

    for key, variant in benign_clinvar_variants.items():
        assert key in expected_benign_clinvar_reasons
        assert expected_benign_clinvar_reasons[key] == variant.reason

    # print(clinvar_variants)
    # print(benign_clinvar_variants)



# def test_compare_clinvar(tmp_path: Path):

#     new_clinvar = {}
#     old_clinvar = {}
#     final_bed_fp = ''
#     out_dir = ''

#     compare_clinvar(
        
#     )