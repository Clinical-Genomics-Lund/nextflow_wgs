# create GATK ref

Build intervals from reference, 1000bp intervals, -XL exclude areas, centromeres, PAR and such

`gatk PreprocessIntervals -R $FASTA_FILE --padding 0 -imr OVERLAPPING_ONLY -O grch38.preprocessed.blacklisted.interval_list -XL hg38-blacklist.v2.bed(GENCODE)`

Collect read counts from COHORT SAMPLES, need to do for GC-content and extreme counts annotation later. Did this for around 200 samples

`gatk CollectReadCounts -R $FASTA_FILE -imr OVERLAPPING_ONLY --format TSV -O sampleN.tsv -I sampleN.bam -L grch38.preprocessed.blacklisted.interval_list`

Annotate reference with GC-content

`gatk AnnotateIntervals -L grch38.preprocessed.blacklisted.interval_list -R $FASTA_FILE -imr OVERLAPPING_ONLY -O grch38.annotated.tsv`

Filter out GC rich outliers and exreme counts from cohort

`gatk FilterIntervals -L grch38.preprocessed.blacklisted.interval_list --annotated-intervals grch38.annotated.tsv -imr OVERLAPPING_ONLY -O grch38.preprocessed.blacklisted.gcfiltered.interval_list -I sample1.tsv -I sample2.tsv ... -I sampleN.tsv`

scatter intervals, Divided genome into 7. Might need to go even more depending on performance of cluster

`gatk IntervalListTools --INPUT grch38.preprocessed.blacklisted.gcfiltered.interval_list --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_CONTENT 15000 --OUTPUT scatter`

Priors file needed in next step. Basically taking care of X and Y at different expected ploidys

```
CONTIG_NAME     PLOIDY_PRIOR_0  PLOIDY_PRIOR_1  PLOIDY_PRIOR_2  PLOIDY_PRIOR_3
1       0.01    0.01    0.97    0.01
2       0.01    0.01    0.97    0.01
3       0.01    0.01    0.97    0.01
4       0.01    0.01    0.97    0.01
5       0.01    0.01    0.97    0.01
6       0.01    0.01    0.97    0.01
7       0.01    0.01    0.97    0.01
8       0.01    0.01    0.97    0.01
9       0.01    0.01    0.97    0.01
10      0.01    0.01    0.97    0.01
11      0.01    0.01    0.97    0.01
12      0.01    0.01    0.97    0.01
13      0.01    0.01    0.97    0.01
14      0.01    0.01    0.97    0.01
15      0.01    0.01    0.97    0.01
16      0.01    0.01    0.97    0.01
17      0.01    0.01    0.97    0.01
18      0.01    0.01    0.97    0.01
19      0.01    0.01    0.97    0.01
20      0.01    0.01    0.97    0.01
21      0.01    0.01    0.97    0.01
22      0.01    0.01    0.97    0.01
X       0.01    0.49    0.49    0.01
Y       0.50    0.50    0.00    0.00
```

call ploidy, needed for CNVcaller in next step, COHORT MODE, did this for 30 samples

```
gatk DetermineGermlineContigPloidy \
    -L grch38.preprocessed.blacklisted.gcfiltered.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --contig-ploidy-priors priors \
    --output . \
    --output-prefix ploidy \
    -I sample1.tsv ... -I sampleN.tsv
```

call GermlineCNVCaller per scatter, COHORT MODE. Did this for 30 samples, for 7 scatters


```
for i in $( ls scatter_fewer); 
    do gatk --java-options "-Djava.io.tmpdir=/local/tmp_gatk" \
        GermlineCNVCaller --run-mode COHORT \
        -L scatter_7/${i}/scattered.interval_list \
        --contig-ploidy-calls ploidy-calls \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output cohort30_fewer \
        --output-prefix cohort_${i} \
        -I sample1.tsv -I sample2.tsv ... -I sampleN.tsv;
        done
```

## Call per case

```
 ## collect read counts ##
 export THEANO_FLAGS="base_compiledir=."
 export MKL_NUM_THREADS=${task.cpus}
 export OMP_NUM_THREADS=${task.cpus}
 gatk --java-options "-Xmx20g" CollectReadCounts \\
     -L $params.gatk_intervals \\
     -R $params.genome_file \\
     -imr OVERLAPPING_ONLY \\
     -I $bam \\
     --format TSV -O ${id}.tsv

 ## call ploidy  per case ##
 export THEANO_FLAGS="base_compiledir=."
 export MKL_NUM_THREADS=${task.cpus}
 export OMP_NUM_THREADS=${task.cpus}
 gatk --java-options "-Xmx20g" DetermineGermlineContigPloidy \\
     --model $params.ploidymodel \\
     -I $tsv \\
     -O ploidy/ \\
     --output-prefix $group
 tar -cvf ploidy.tar ploidy/

 ## CNVcaller PER SCATTER ##
 export THEANO_FLAGS="base_compiledir=."
 export HOME=/local/scratch
 export MKL_NUM_THREADS=${task.cpus}
 export OMP_NUM_THREADS=${task.cpus}
 tar -xvf ploidy.tar
 mkdir ${group}_${i}
 gatk --java-options "-Xmx25g" GermlineCNVCaller \\
     --run-mode CASE \\
     -I $tsv \\
     --contig-ploidy-calls ploidy/${group}-calls/ \\
     --model ${refpart} \\
     --output ${group}_${i}/ \\
     --output-prefix ${group}_${i}
 tar -cvf ${group}_${i}.tar ${group}_${i}/

 ## post process ##
 THEANO_FLAGS="base_compiledir=/fs1/resources/theano"
 for model in !{tar}; do
 tar -xvf $model
 done
 tar -xvf !{ploidy}
 export MKL_NUM_THREADS=!{task.cpus}
 export OMP_NUM_THREADS=!{task.cpus}
 gatk --java-options "-Xmx25g" PostprocessGermlineCNVCalls \
     --allosomal-contig X --allosomal-contig Y \
     --contig-ploidy-calls ploidy/!{group}-calls/ \
     --sample-index 0 \\
     --output-genotyped-intervals genotyped-intervals-!{group}-vs-cohort30.vcf.gz \
     --output-genotyped-segments genotyped-segments-!{group}-vs-cohort30.vcf.gz \
     --output-denoised-copy-ratios denoised-!{group}-vs-cohort30.vcf.gz \
     --sequence-dictionary !{params.GENOMEDICT} \
     --calls-shard-path !{caseshards} \
     --model-shard-path !{modelshards}
```
