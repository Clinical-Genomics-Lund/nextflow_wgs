### Annotation files

To run the pipeline, you need to setup a number of annotation files. Default values are specified in the `configs/nextflow.hopper.config`. 

<!-- Organize better. For instance, the set of references required to run VEP: -->

<!-- - VEP cache
- VEP-fasta
- gnomAD (exomes, genomes, mitochondria) VCFs
- phyloP
- phastcons
- CADD
- maxentscan plugin. -->

| Parameter        | Function             | Format  | Description                                                                                                                                                                                                                                             |
| ---------------- | -------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `genome_file`    | Reference genome     | `fasta` | The main reference                                                                                                                                                                                                                                      |
| `GENOMEDICT`     | Reference genome     | `dict`  | Dict file for the FASTA file                                                                                                                                                                                                                            |
| `rCRS_fasta`     | Reference genome     | `fasta` | Mitochondrial FASTA reference sequence found [here](https://www.ncbi.nlm.nih.gov/nuccore/251831106)                                                                                                                                                     |
| `bwa_shards`     | Alignment            |         | Number of shards to split the reads into prior to alignment                                                                                                                                                                                             |
| `shardbwa`       | Alignment            |         | Boolean specifying whether to do alignment in sharded mode                                                                                                                                                                                              |
| `KNOWN`          | Alignment            |         | Gold-standard indels used in Sentieon's base quality score recalibration (BQSR)                                                                                                                                                                         |
| `VEP_CACHE`      | Annotation (VEP)     |         | VEP files for offline run. Instructions on how to setup a cache can be found [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).                                                                                           |
| `VEP_FASTA`      | Annotation (VEP)     | `fasta` | Reference sequence used for optimization within the VEP cache                                                                                                                                                                                           |
| `SYNONYMS`       | Annotation (VEP)     | `tsv`   | Chromosome synonyms used by VEP (for instance to recognize `M` as `MT`)                                                                                                                                                                                 |
| `CADD`           | Annotation (VEP)     | `tsv`   | Precalculated CADD indices. ([download page](https://cadd.gs.washington.edu/download)) [CADD](https://cadd.gs.washington.edu) is used to score deleteriousness of SNVs and indels in human.                                                             |
| `GNOMAD_GENOMES` | Annotation (VEP)     | `tsv`   | Allele frequencies for variants ([download page](https://gnomad.broadinstitute.org/downloads))                                                                                                                                                          |
| `GNOMAD_EXOMES`  | Annotation (VEP)     | `tsv`   | FIXME                                                                                                                                                                                                                                                   |
| `GNOMAD_MT`      | Annotation (VEP)     | `tsv`   | FIXME                                                                                                                                                                                                                                                   |
| `MAXENTSCAN`     | Annotation (VEP)     |         | A direct path to the MaxEntScan scripts (FIXME: Cannot this be specified just by )                                                                                                                                                                      |
| `dbNSFP`         | Annotation (VEP)     | `tsv`   | (FIXME): Path to gzipped tab delimited file containing the dbNSFP annotations.                                                                                                                                                                          |
| `PHYLOP`         | Annotation (VEP)     | `tsv`   | Measures evolutionary conservation at individual alignment sites. This parameter expects a tab-delimited file with pre-calculated phyloP scores.                                                                                                        |
| `PHASTCONS`      | Annotation (VEP)     | `tsv`   | Conservation scores. This parameter expects a tab-delimited file with pre-calculated phyloP scores.                                                                                                                                                     |
| `vcfanno`        | Annotation (VCFAnno) | `toml`  | Config file pointing to which files to retrieve annotation from, which fields that should be extracted and how these should be inserted in the target VCF. This additional annotations are different depending on which profile is used for processing. |
| `LUA`            | Annotation (VCFAnno) | `lua`   | Custom script for VCFAnno annotations                                                                                                                                                                                                                   |

### Preparing annotation files

Some pointers on how the reference files can be prepared.

#### Preparing reference and indexes

The pipeline is running with the hg38 reference genome, but with no `chr` prefix. 

The `chr` prefix from each FASTA entry.

```
sed 's/^>chr/>/' <fasta>.fna > <fasta>_nochr.fna
```

Generate a bwa index. This requires that you have [bwa](https://github.com/lh3/bwa) installed.

```
bwa index <fasta>_nochr.fna
```

Create the sequence dictionary. This requires [Picard Tools](https://broadinstitute.github.io/picard/).

```
java -jar picard.jar CreateSequenceDictionary \
    REFERENCE=<fasta>_nochr.fna \
    OUTPUT=<fasta>_nochr.dict
```

Generate the samtools fasta index using [samtools](http://www.htslib.org/).

```
samtools faidx <fasta>_nochr.fna
```

##### VCF anno

<!-- VCFanno is a tool used for rapidly load annotations into the `INFO` field of a VCF file. The fields can be directly retrieved from tab-delimited input files, or preprocessed using lua-script.

* `LUA`: FIXME: Why do we need this one? Should it be a param, or maybe better stored within the repo?
* `vcfanno`: Config file pointing to which files to retrieve annotation from, which fields that should be extracted and how these should be inserted in the target VCF. This additional annotations are different depending on which profile is used for processing. -->

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
