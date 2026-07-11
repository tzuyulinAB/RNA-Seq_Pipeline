# Metatranscriptomics Nextflow Pipeline

## Workflow

1. Validate the sample sheet.
2. Report command availability.
3. Prepare the selected Trimmomatic adapter FASTA.
4. Prepare a SortMeRNA rRNA reference FASTA.
5. Build or reuse a SortMeRNA index.
6. Trim paired-end RNA reads with Trimmomatic.
7. Run FastQC on paired trimmed reads.
8. Remove rRNA reads with SortMeRNA.
9. Optionally map rRNA-removed reads to dereplicated MAGs with BBMap.

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
  --ref_dir /path/to/metagenomics/results/mags/drep/dereplicated_genomes
```

Use existing resource files instead of preparing/downloading them:

```bash
nextflow run main.nf -profile docker \
  --adapter_fasta /path/to/NexteraPE-PE.fa \
  --sortmerna_ref /path/to/smr_v4.3_default_db.fasta
```

Use the smaller SortMeRNA fast database for a test run:

```bash
nextflow run main.nf -profile docker \
  --sortmerna_ref resources/sortmerna/smr_v4.3_fast_db.fasta
```

When the SortMeRNA database archive is downloaded, all included database FASTA files are extracted into `resources/sortmerna/`. The pipeline uses the file selected by `--sortmerna_ref`; if not assigned, it uses `resources/sortmerna/smr_v4.3_default_db.fasta`.

Reuse an existing SortMeRNA index directory:

```bash
nextflow run main.nf -profile docker \
  --sortmerna_ref /path/to/smr_v4.3_default_db.fasta \
  --sortmerna_index_dir /path/to/sortmerna_index
```

## Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--samples` | `config/samples.tsv` | TSV with `sample_id`, `condition`, `read1`, and `read2`. |
| `--outdir` | `results` | Main output directory. |
| `--adapter_selected` | `NexteraPE-PE.fa` | Trimmomatic adapter file to prepare. |
| `--adapter_fasta` | `null` | Existing adapter FASTA to use directly. |
| `--sortmerna_ref` | `resources/sortmerna/smr_v4.3_default_db.fasta` | Existing or target SortMeRNA reference FASTA. |
| `--sortmerna_index_dir` | `null` | Existing SortMeRNA work/index directory to reuse. If omitted, the workflow reuses `resources/sortmerna/index/sortmerna_index` when present, otherwise builds one once. |
| `--ref_dir` | empty | Optional directory containing reference `.fa`, `.fna`, or `.fasta` files for BBMap expression mapping. |

## Outputs

| Step | Outputs |
| --- | --- |
| Validation | `config/validation.ok`, `logs/config/validate_config.log` |
| Dependency check | `reports/dependency_check.tsv`, `logs/config/check_dependencies.log` |
| SortMeRNA index | `resources/sortmerna/index/sortmerna_index`, `logs/resources/sortmerna_index.log` |
| Trimming | `results/rna/trim/*_paired.fq.gz`, `results/rna/trim/*_unpaired.fq.gz` |
| FastQC | `results/rna/fastqc_trimmed/*_fastqc.html` |
| rRNA removal | `results/rna/sortmerna/*_rRNArm_fwd.fq.gz`, `results/rna/sortmerna/*_rRNArm_rev.fq.gz` |
| BBMap expression | `results/expression/bbmap/<sample>/<genome>/scafstats.tsv`, `covstats.tsv` |
| Assigned read table | `results/expression/assignedReads_feature_table.tsv` |
| Gene length table | `results/reference/gene_lengths.tsv` |

## Platform Note

The workflow is configured for Docker only. Docker Desktop must be installed and running before launching the pipeline.
