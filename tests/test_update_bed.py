from unittest.mock import patch
import sys
from pathlib import Path

# sys.path.insert(0, str(Path(__file__).resolve().parent / 'bin'))

from src.update_bed import write_base

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


def test_write_base(tmp_path: str):

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

# def test_full(tmp_path: str):
#     print("Test")
#     script_path = os.path.abspath('bin/reference_tools/update_bed.py')
#     args = [
#         script_path,
#         '--old', 'data/testdata/clinvar38_20231230.vcf',
#         '--new', 'data/testdata/clinvar38_20240624.vcf',
#         '--release', '108',
#         '--clinvardate', '2024-01-01',
#         '--ensembl', os.path.abspath('Homo_sapiens.GRCh38.108.gtf'),
#         '--skip_download',
#         '--out_dir', 'testout'
#     ]

#     result = subprocess.run(['python3'] + args, cwd=tmp_path, capture_output=True, text=True)
#     print(result)
#     # assert result.returncode == 0

#     # expected_output_files = ['output1.bed', 'output2.bed']
#     # for file_name in expected_output_files:
#     #     assert os.path.exists(tmp_path / file_name)
