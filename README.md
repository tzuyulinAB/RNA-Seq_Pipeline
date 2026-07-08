# Metatranscriptomics Nextflow Pipeline

This DSL2 Nextflow workflow follows the reference workflow in `Snakemake_ref/Snakefile`.

## Workflow

1. Validate the sample sheet.
2. Report command availability.
3. Prepare the selected Trimmomatic adapter FASTA.
4. Prepare a SortMeRNA rRNA reference FASTA.
5. Trim paired-end RNA reads with Trimmomatic.
6. Run FastQC on paired trimmed reads.
7. Remove rRNA reads with SortMeRNA.
8. Optionally map rRNA-removed reads to dereplicated MAGs with BBMap.

## Inputs

Edit `config/samples.tsv`:

```tsv
sample_id	condition	read1	read2
CAAI_1	AD	raw_data/sample_R1.fastq.gz	raw_data/sample_R2.fastq.gz
```

Relative FASTQ paths are resolved from the directory where you launch Nextflow.

## Run

Dry-run:

```bash
nextflow run main.nf -preview
```

Run with Docker:

```bash
nextflow run main.nf -profile docker
```

Include optional BBMap expression mapping against dereplicated genomes:

```bash
nextflow run main.nf -profile docker \
  --drep_genomes_dir /path/to/metagenomics/results/mags/drep/dereplicated_genomes
```

Use existing resource files instead of preparing/downloading them:

```bash
nextflow run main.nf -profile docker \
  --adapter_fasta /path/to/NexteraPE-PE.fa \
  --sortmerna_ref /path/to/smr_v4.3_default_db.fasta
```

## Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--samples` | `config/samples.tsv` | TSV with `sample_id`, `condition`, `read1`, and `read2`. |
| `--outdir` | `results` | Main output directory. |
| `--adapter_selected` | `NexteraPE-PE.fa` | Trimmomatic adapter file to prepare. |
| `--adapter_fasta` | `null` | Existing adapter FASTA to use directly. |
| `--sortmerna_ref` | `resources/sortmerna/smr_v4.3_default_db.fasta` | Existing or target SortMeRNA reference FASTA. |
| `--drep_genomes_dir` | empty | Optional directory containing dereplicated MAG `.fa` files for BBMap. |

## Outputs

| Step | Outputs |
| --- | --- |
| Validation | `config/validation.ok`, `logs/config/validate_config.log` |
| Dependency check | `reports/dependency_check.tsv`, `logs/config/check_dependencies.log` |
| Trimming | `results/rna/trim/*_paired.fq.gz`, `results/rna/trim/*_unpaired.fq.gz` |
| FastQC | `results/rna/fastqc_trimmed/*_fastqc.html` |
| rRNA removal | `results/rna/sortmerna/*_rRNArm_fwd.fq.gz`, `results/rna/sortmerna/*_rRNArm_rev.fq.gz` |
| BBMap expression | `results/expression/bbmap/<sample>/<genome>/scafstats.tsv`, `covstats.tsv` |

## Platform Note

The workflow is configured for Docker only. Docker Desktop must be installed and running before launching the pipeline.
