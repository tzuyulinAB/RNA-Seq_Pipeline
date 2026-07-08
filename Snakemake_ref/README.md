# Metatranscriptomics Snakemake Workflow

This workflow contains the RNA-specific pieces moved out of the metagenomics workflow:

1. RNA paired-end trimming with Trimmomatic.
2. Trimmed-read FastQC.
3. rRNA removal with SortMeRNA.
4. Optional BBMap expression mapping against dereplicated MAGs from the DNA workflow.

## Configure

Edit `config/samples.tsv`:

```tsv
sample_id	condition	read1	read2
CAAI_1	AD	/absolute/path/sample_R1.fastq.gz	/absolute/path/sample_R2.fastq.gz
```

Edit `config/config.yaml` if needed:

- `adapters.selected`: Trimmomatic adapter FASTA to use, for example `NexteraPE-PE.fa`.
- `databases.sortmerna_ref`: defaults to `resources/sortmerna/smr_v4.3_default_db.fasta`.
- `databases.drep_genomes_dir`: optional path to DNA workflow dereplicated genomes, such as `/path/to/metagenomics/results/mags/drep/dereplicated_genomes`.
- `resources.tmpdir`: optional scratch directory. Leave empty to use `results/tmp`.

The workflow will download the selected Trimmomatic adapters and the SortMeRNA v4.3.4 default database if they are not already present.

## Run

Dry-run:

```bash
./scripts/run_snakemake_local.sh -n --cores 8 --config samples=config/samples_test.tsv
```

Run:

```bash
./scripts/run_snakemake_local.sh --cores 8 --config samples=config/samples_test.tsv
```

To include BBMap expression mapping, set `databases.drep_genomes_dir`:

```bash
./scripts/run_snakemake_local.sh --cores 12 \
  --config samples=config/samples.tsv \
           databases.drep_genomes_dir=/path/to/metagenomics/results/mags/drep/dereplicated_genomes
```

## macOS Note

Native Apple Silicon Conda did not have a Bioconda `sortmerna` package during the local test. The RNA workflow is therefore expected to run best on Linux, AWS EC2/AWS Batch, HPC, or a compatible container environment.
