These parameters are specified in the `configs/nextflow.hopper.config`. Some of these are specified as shared for all profiles, and some have profile specific settings.

Should we keep the profile specific ones out from this?

Further hard-coded annotations within the script?

Organize better. For instance, the set of references required to run VEP:

- VEP cache
- VEP-fasta
- gnomAD (exomes, genomes, mitochondria) VCFs
- phyloP
- phastcons
- CADD
- maxentscan plugin.

##### Reference genome

Related params:

* `genome_file`: Path to the `.fna` file (`<fasta>_nochr.fna` below)
* `GENOMEDICT`: Path to the `.dict` file
* `rCRS_fasta`: A mitochondria FASTA reference sequence as found here: https://www.ncbi.nlm.nih.gov/nuccore/251831106. Standard sequence to report mitochondria. It is used to normalize positions of mitochondria indels. FIXME: Verify

The pipeline is running with the hg38 reference genome, but with no `chr` prefix. It can be prepared the following way:

Retrieve the genome (FIXME: Where to get the masked genome?)

Uncompress.

```
gunzip <fasta>.fna.gz
```

Remove the `chr` prefix from each FASTA entry.

```
sed 's/^>chr/>/' <fasta>.fna > <fasta>_nochr.fna
```

Generate a `bwa` index.

```
bwa index <fasta>_nochr.fna
```

Create the sequence dictionary.

```
java -jar picard.jar CreateSequenceDictionary \
    REFERENCE=<fasta>_nochr.fna \
    OUTPUT=<fasta>_nochr.dict
```

Generate the samtools fasta index.

```
samtools faidx <fasta>_nochr.fna
```
##### Alignment and base quality score recalibration

FIXME: Ask about the purpose of the sharding. I guess it is to split the alignment to multiple processes (which could allow using multiple nodes perhaps). But let's verify.

The alignment can optionally be sharded, aligning subsets of the FASTQ reads in different Nextflow processes. 

* `bwa_shards`: The number of shards to split the reads into.
* `shardbwa`:  A boolean specifying whether to do alignment in sharded mode. 
* `KNOWN`: Set of gold standard indels to use in Sentieon's base quality score recalibration (BQSR)
##### VEP

https://www.ensembl.org/info/docs/tools/vep/index.html

VEP is a commonly used software to determine the effect of the variants (Variant Effect Predictor). It can be run as a web service, or offline using a cache. Some of the predicted information includes what genes and transcripts are affected, whether the variants are located in certain features (i.e. in coding sequence, in intron and so on), the consequence on the protein sequence, and changes in the protein sequence.

The results are added to the `CSQ` field in the `INFO` column of the provided VCF. One CSQ entry is added per matching transcript feature, meaning that while for other `INFO` annotations there will be just one value, for the VEP there can be a set of consequence values for instances.

* `VEP_CACHE`: Used for offline run and contains information about transcript locations, gene identifiers, regulatory regions, prediction scores for tools such as SIFT and PolyPhem. (https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html). This cache can be downloaded from the ENSEMBL web page.
* `VEP_FASTA`: Used to look up reference sequence (FIXME: Why do we use it? Why not the same FASTA as the actual reference?)
* `SYNONYMS`: Used to allow VEP to recognize mitochondria `M` chromosomes using the `MT` string (which is currently used within CMD)

**CADD** 

* `CADD` parameter

CADD is a tool for scoring deleteriousness of SNVs and indels in human (https://cadd.gs.washington.edu/). The specified file contains precalculated scores for specific variants at specified locations. It is only valid for SNVs, indels are handled separately.

First lines of the file:

```
#Chrom  Pos     Ref     Alt     RawScore        PHRED
1       10001   T       A       0.702541        8.478                                   
1       10001   T       C       0.750954        8.921                                   
1       10001   T       G       0.719549        8.634                                   
1       10002   A       C       0.713993        8.583
```

This file can be downloaded from the following page: https://cadd.gs.washington.edu/download . It will need to be indexed before use. Verify the md5sum of the retrieved file against the md5sums provided at the same page.

**Gnomad**

Gnomad is a database with allele frequencies for different variants, i.e. how commonly encountered they are in a general population. Variants causing rare diseases are expected to be rare.

The download page for Gnomad is available here: https://gnomad.broadinstitute.org/downloads

For this pipeline we don't need all the Gnomad annotations. Before usage, the VEP information and most of the INFO fields are removed. Retained are `nhomalt`, `AF`, `popomax` and `AF_popmax` (FIXME: Double check).

```
1       12198   rs62635282      G       C       9876.24 AC0     nhomalt=0               
1       12237   rs1324090652    G       A       81.96   AC0     nhomalt=0               
1       12259   rs1330604035    G       C       37.42   AC0     AF=0;nhomalt=0          
1       12266   rs1442951560    G       A       2721.48 AC0     nhomalt=0               
1       12272   rs1281272113    G       A       2707.42 AC0     AF=0;nhomalt=0
1       13499   rs1461668262    C       A       7144.56 PASS    AF=2.19603e-05;nhomalt=0;popmax=eas;AF_popmax=0.00027819
```

* `GNOMAD_GENOMES`: Gnomad frequencies for the full genome (FIXME)
* `GNOMAD_EXOMES`: Gnomad frequencies reduced to only show those residing within exons (FIXME)
* `GNOMAD_MT`: Gnomad frequencies for mitochondria (FIXME)

**MaxEntScan**

MaxEntScan is used to identify splice site predictions.

* `MAXENTSCAN`: A direct path to the MaxEntScan scripts (FIXME: Cannot this be specified just by )

**dbNSFP**

Used to retrieve data for missense variants, with functional predictions.

* `DBNSFP` (FIXME): Path to gzipped tab delimited file containing the dbNSFP annotations.

**Custom annotations**

Custom annotations are directly loaded from tab-delimited files.

* `PHYLOP`: Measures evolutionary conservation at individual alignment sites. This parameter expects a tab-delimited file with pre-calculated phyloP scores.
* `PHASTCONS`: Conservation scores. This parameter expects a tab-delimited file with pre-calculated phyloP scores.
##### VCF anno

VCFanno is a tool used for rapidly load annotations into the `INFO` field of a VCF file. The fields can be directly retrieved from tab-delimited input files, or preprocessed using lua-script.

* `LUA`: FIXME: Why do we need this one? Should it be a param, or maybe better stored within the repo?
* `vcfanno`: Config file pointing to which files to retrieve annotation from, which fields that should be extracted and how these should be inserted in the target VCF. This additional annotations are different depending on which profile is used for processing.

Example `vcf_anno` as used by CMD. Replace `<base_dir>` with the base folder for each file. Here the following information is retrieved:

* Frequency information from `loqusdb`
* Clinical significance information and reviewer quality from `ClinVar`
* Ensemble gene annotations

```
[[annotation]]
file="<base_dir>/loqusdb_latest.vcf.gz"
fields=["Frq"]
names=["loqusdb_freq"]
ops=["concat"]

[[annotation]]
file="<base_dir>/clinvar38_latest.vcf.gz"
fields=["CLNSIG","CLNREVSTAT","CLNACC"]
names=["CLNSIG","CLNREVSTAT","CLNACC"]
ops=["first","first","first"]

[[annotation]]
file="<base_dir>/ensembl_genes_38.sort.padded200bp.noNeg.bed.gz"
columns=[4]
names=["Annotation"]
ops=["concat"]
```

##### Copy number variants using GATK

* `COV_INTERVAL_LIST`: List of pre-calculated intervals to use in `gatk CollectReadCounts`, which will calculate how many reads have mapped within each of these bins. These bins can be calculated using the `GATK PreprocessIntervals` command (https://gatk.broadinstitute.org/hc/en-us/articles/360037427371-PreprocessIntervals-BETA-)

##### Coverage profile using GATK

* `GATK_PON_MALE` and `GATK_PON_FEMALE`: Gender specific panels of normals used in the `DenoiseReadCounts` (FIXME: How are these prepared?)

##### Short tandem repeats (STR)

These are called using `expansionhunter` and annotated using `stranger`.

* `expansionhunter_catalog`: JSON file outlining variants to look for and their genomic locations. This catalog can be created manually, or downloaded from the ExpansionHunter repository at: https://github.com/Illumina/ExpansionHunter/tree/master/variant_catalog/grch38
* `expansionhunter_catalog_gav`
##### Rank models

* `rank_model_s`
* `rank_model`
* `svrank_model_s`
* `svrank_model`
##### Onco-specific parameters

* `gene_regions`: Bed-file with genes of interest specifically for onco (FIXME)
##### Other params

* `FASTGNOMAD_REF`: Path to .dat file containing only highly frequent locations. Used for filtering out common variants. (FIXME: But how is that dat-file generated?)
* `GENS_GNOMAD`: Tab-separated Gnomad frequencies specifically for Gens
* `loqusdb`: Path-name to the `loqusdb` from which to retrieve local artefact frequencies
* `panelsdef`: Location of panels used in Scout (FIXME: What exactly is this?)
* `scoutbed`: Bedfile used in coverage calculations (FIXME: What is this?)
* `svdb`:  Location of artefact data bases used to identify artefacts within the data, which is used in prescoring (FIXME)
