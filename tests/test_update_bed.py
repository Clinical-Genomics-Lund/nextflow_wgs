import pytest
import os
import subprocess


def test_full(tmp_path: str):
    print("Test")
    script_path = os.path.abspath('bin/reference_tools/update_bed.py')
    args = [
        script_path,
        '--old', 'data/testdata/clinvar38_20231230.vcf',
        '--new', 'data/testdata/clinvar38_20240624.vcf',
        '--release', '108',
        '--clinvardate', '2024-01-01',
        '--ensembl', os.path.abspath('Homo_sapiens.GRCh38.108.gtf'),
        '--skip_download',
        '--out_dir', 'testout'
    ]

    result = subprocess.run(['python3'] + args, cwd=tmp_path, capture_output=True, text=True)
    print(result)
    # assert result.returncode == 0

    # expected_output_files = ['output1.bed', 'output2.bed']
    # for file_name in expected_output_files:
    #     assert os.path.exists(tmp_path / file_name)
