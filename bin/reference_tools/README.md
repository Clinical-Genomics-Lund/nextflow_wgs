Legacy, how to run the Perl script:

```
perl bin/reference_tools/update_bed.pl \
    --old data/testdata/clinvar38_20231230.vcf \
    --new data/testdata/clinvar38_20240624.vcf \
    --build 108 \
    --clinvardate 20240624 \
    --skip_download \
    --incl_bed bin/reference_tools/agilient_hg38_nochr_noalt_1-3.bed
```

How to run the Python script:

```
python3 bin/reference_tools/update_bed.py \
    --old data/testdata/clinvar38_20231230.vcf \
    --new data/testdata/clinvar38_20240624.vcf \
    --release 108 \
    --clinvardate 20240624 \
    --out_dir out/4 \
    --skip_download \
    --incl_bed bin/reference_tools/agilient_hg38_nochr_noalt_1-3.bed \
    --keep_tmp
```


