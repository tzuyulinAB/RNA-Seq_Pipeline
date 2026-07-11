#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Metatranscriptomics RNA workflow translated from Snakemake_ref/Snakefile.
 */

def readSamples(String samplesheet) {
    def path = file(samplesheet)
    if (!path.exists()) {
        throw new IllegalArgumentException("Samplesheet not found: ${samplesheet}")
    }

    def records = []
    def header = null

    path.eachLine { line ->
        def trimmed = line.trim()
        if (!trimmed || trimmed.startsWith('#')) {
            return
        }

        def fields = line.split('\t', -1).collect { it.trim() }
        if (header == null) {
            header = fields
            return
        }

        def row = [:]
        header.eachWithIndex { key, idx ->
            row[key] = idx < fields.size() ? fields[idx] : ''
        }
        if (row.values().any { it }) {
            records << row
        }
    }

    def required = ['sample_id', 'read1', 'read2']
    def missing = required.findAll { !(header ?: []).contains(it) }
    if (missing) {
        throw new IllegalArgumentException("Samplesheet is missing required column(s): ${missing.join(', ')}")
    }

    return records
}

def existingPathChannel(String path, boolean checkIfExists = true) {
    Channel.fromPath(path, checkIfExists: checkIfExists)
}

process VALIDATE_CONFIG {
    tag 'samples'
    label 'python'

    publishDir '.', mode: 'copy', overwrite: true

    input:
    path samplesheet

    output:
    path 'config/validation.ok', emit: ok
    path 'logs/config/validate_config.log', emit: log

    script:
    """
    mkdir -p config logs/config
    python3 "${projectDir}/bin/validate_config.py" --samples ${samplesheet} --base-dir "${workflow.launchDir}" > logs/config/validate_config.log 2>&1 || {
      cat logs/config/validate_config.log >&2
      exit 1
    }
    touch config/validation.ok
    """
}

process CHECK_DEPENDENCIES {
    tag 'tools'
    label 'python'

    publishDir '.', mode: 'copy', overwrite: true

    input:
    val include_bbmap

    output:
    path 'reports/dependency_check.tsv', emit: report
    path 'logs/config/check_dependencies.log', emit: log

    script:
    def bbmap_flag = include_bbmap ? '--include-bbmap' : ''
    """
    mkdir -p reports logs/config
    python3 "${projectDir}/bin/check_dependencies.py" --output reports/dependency_check.tsv ${bbmap_flag} > logs/config/check_dependencies.log 2>&1 || {
      cat logs/config/check_dependencies.log >&2
      exit 1
    }
    """
}

process PREPARE_TRIMMOMATIC_ADAPTER {
    tag params.adapter_selected
    label 'download'

    publishDir params.adapter_dir, mode: 'copy', overwrite: true, pattern: params.adapter_selected
    publishDir 'logs/resources', mode: 'copy', overwrite: true, pattern: 'trimmomatic_adapters.log'

    output:
    path params.adapter_selected, emit: adapter
    path 'trimmomatic_adapters.log', emit: log

    script:
    """
    python3 "${projectDir}/bin/prepare_trimmomatic_adapters.py" \
      --adapter-dir . \
      --selected ${params.adapter_selected} \
      > trimmomatic_adapters.log 2>&1 || {
      cat trimmomatic_adapters.log >&2
      exit 1
    }
    """
}

process DOWNLOAD_SORTMERNA_DATABASE {
    tag 'sortmerna_db'
    label 'download'

    publishDir params.sortmerna_dir, mode: 'copy', overwrite: true, pattern: 'smr_v4.3_*.fasta'
    publishDir 'logs/resources', mode: 'copy', overwrite: true, pattern: 'sortmerna_database.log'

    output:
    path 'selected_sortmerna_ref.fasta', emit: ref
    path 'smr_v4.3_*.fasta', emit: databases
    path 'sortmerna_database.log', emit: log

    script:
    """
    python3 "${projectDir}/bin/download_sortmerna_database.py" \
      --outdir . \
      --output selected_sortmerna_ref.fasta \
      --selected ${params.sortmerna_ref.toString().tokenize('/').last()} \
      > sortmerna_database.log 2>&1 || {
      cat sortmerna_database.log >&2
      exit 1
    }
    """
}

process PREPARE_SORTMERNA_INDEX {
    tag 'sortmerna_index'
    label 'sortmerna'
    cpus { params.threads.sortmerna as int }

    publishDir params.sortmerna_dir, mode: 'copy', overwrite: true, pattern: 'index/**'
    publishDir 'logs/resources', mode: 'copy', overwrite: true, pattern: 'sortmerna_index.log'

    input:
    path ref

    output:
    path 'index', emit: index
    path 'sortmerna_index.log', emit: log

    script:
    """
    cat > index_probe_R1.fq <<'EOF'
@index_probe/1
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

    cat > index_probe_R2.fq <<'EOF'
@index_probe/2
TGCATGCATGCATGCATGCATGCATGCATGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

    sortmerna \
      --workdir index \
      --ref ${ref} \
      --reads index_probe_R1.fq \
      --reads index_probe_R2.fq \
      --aligned index_probe_aligned \
      --other index_probe_other \
      --fastx \
      --threads ${task.cpus} \
      --out2 \
      --sout \
      > sortmerna_index.log 2>&1 || {
      cat sortmerna_index.log >&2
      exit 1
    }

    test -d index/idx
    rm -rf index/kvdb index/readb
    mkdir -p index/kvdb index/readb
    """
}

process TRIM_RNA {
    tag sample_id
    label 'trimmomatic'
    cpus { params.threads.trimmomatic as int }
    memory { params.memory.trimmomatic }

    publishDir "${params.outdir}/rna/trim", mode: 'copy', overwrite: true, pattern: '*.fq.gz'
    publishDir 'logs/trimmomatic/rna', mode: 'copy', overwrite: true, pattern: '*.log'

    input:
    tuple val(sample_id), val(condition), path(read1), path(read2), path(adapter)

    output:
    tuple val(sample_id), val(condition),
          path("${sample_id}_F_paired.fq.gz"),
          path("${sample_id}_R_paired.fq.gz"),
          path("${sample_id}_F_unpaired.fq.gz"),
          path("${sample_id}_R_unpaired.fq.gz"),
          emit: trimmed
    path "${sample_id}.log", emit: log

    script:
    """
    export JAVA_TOOL_OPTIONS="${params.trimmomatic_java_opts}"
    trimmomatic PE -threads ${task.cpus} \
      ${read1} ${read2} \
      ${sample_id}_F_paired.fq.gz ${sample_id}_F_unpaired.fq.gz \
      ${sample_id}_R_paired.fq.gz ${sample_id}_R_unpaired.fq.gz \
      ILLUMINACLIP:${adapter}:${params.trimmomatic_illuminaclip} \
      ${params.trimmomatic_options} \
      > ${sample_id}.log 2>&1
    """
}

process FASTQC_RNA_TRIMMED {
    tag sample_id
    label 'fastqc'
    cpus { params.threads.fastqc as int }

    publishDir "${params.outdir}/rna/fastqc_trimmed", mode: 'copy', overwrite: true, pattern: '*_fastqc.*'
    publishDir 'logs/fastqc/rna_trimmed', mode: 'copy', overwrite: true, pattern: '*.log'

    input:
    tuple val(sample_id), val(condition), path(fwd), path(rev), path(fwd_unpaired), path(rev_unpaired)

    output:
    tuple val(sample_id), path("${sample_id}_F_paired_fastqc.html"), path("${sample_id}_R_paired_fastqc.html"), emit: html
    path "${sample_id}.log", emit: log

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${fwd} ${rev} > ${sample_id}.log 2>&1
    """
}

process SORTMERNA {
    tag sample_id
    label 'sortmerna'
    cpus { params.threads.sortmerna as int }

    publishDir "${params.outdir}/rna/sortmerna", mode: 'copy', overwrite: true, pattern: '*_rRNArm_*.fq.gz'
    publishDir 'logs/sortmerna', mode: 'copy', overwrite: true, pattern: '*.log'

    input:
    tuple val(sample_id), val(condition), path(fwd), path(rev), path(fwd_unpaired), path(rev_unpaired), path(ref), path(index_dir)

    output:
    tuple val(sample_id), val(condition),
          path("${sample_id}_rRNArm_fwd.fq.gz"),
          path("${sample_id}_rRNArm_rev.fq.gz"),
          emit: rrna_removed
    path "${sample_id}.log", emit: log

    script:
    """
    mkdir -p sortmerna_work_${sample_id}
    cp -aL ${index_dir}/idx sortmerna_work_${sample_id}/idx
    mkdir -p sortmerna_work_${sample_id}/kvdb sortmerna_work_${sample_id}/readb

    sortmerna \
      --workdir sortmerna_work_${sample_id} \
      --ref ${ref} \
      --reads ${fwd} \
      --reads ${rev} \
      --aligned ${sample_id} \
      --other ${sample_id}_rRNArm \
      --fastx \
      --threads ${task.cpus} \
      --paired_in \
      --out2 \
      > ${sample_id}.log 2>&1 || {
      cat ${sample_id}.log >&2
      exit 1
    }

    for fastq in *.fq; do
      if [ -e "\${fastq}" ]; then
        gzip "\${fastq}"
      fi
    done

    if [ -s "${sample_id}_rRNArm_paired_fwd.fq.gz" ] && [ ! -s "${sample_id}_rRNArm_fwd.fq.gz" ]; then
      mv ${sample_id}_rRNArm_paired_fwd.fq.gz ${sample_id}_rRNArm_fwd.fq.gz
    fi

    if [ -s "${sample_id}_rRNArm_paired_rev.fq.gz" ] && [ ! -s "${sample_id}_rRNArm_rev.fq.gz" ]; then
      mv ${sample_id}_rRNArm_paired_rev.fq.gz ${sample_id}_rRNArm_rev.fq.gz
    fi

    test -s ${sample_id}_rRNArm_fwd.fq.gz || {
      echo "Missing expected SortMeRNA output: ${sample_id}_rRNArm_fwd.fq.gz" >&2
      ls -lh >&2
      cat ${sample_id}.log >&2
      exit 1
    }
    test -s ${sample_id}_rRNArm_rev.fq.gz || {
      echo "Missing expected SortMeRNA output: ${sample_id}_rRNArm_rev.fq.gz" >&2
      ls -lh >&2
      cat ${sample_id}.log >&2
      exit 1
    }
    """
}

process BBMAP_EXPRESSION {
    tag sample_id
    label 'bbmap'
    cpus { params.threads.bbmap as int }

    publishDir "${params.outdir}/expression/bbmap", mode: 'copy', overwrite: true, pattern: '*/**'
    publishDir 'logs/bbmap', mode: 'copy', overwrite: true, pattern: '*.log'

    input:
    tuple val(sample_id), val(condition), path(fwd), path(rev), path(genomes)

    output:
    path "${sample_id}/.done", emit: done
    path "${sample_id}/*/scafstats.tsv", emit: scafstats
    path "${sample_id}/*/covstats.tsv", emit: covstats
    tuple val(sample_id), path("${sample_id}/${sample_id}_assignedReads.tsv"), emit: assigned_reads
    path "${sample_id}.log", emit: log

    script:
    def sample_prefix = sample_id.toString().tokenize('_')[0]
    """
    mkdir -p ${sample_id}
    printf "Sample\\tORF_id\\tassignedReads\\n" > ${sample_id}/${sample_id}_assignedReads.tsv
    sample_prefix="${sample_prefix}"

    for genome in ${genomes}; do
      name=\$(basename "\${genome}")
      name="\${name%.*}"
      if [[ "\${name}" == "\${sample_prefix}_"* ]]; then
        out_name="${sample_id}_\${name#\${sample_prefix}_}"
      else
        out_name="${sample_id}_\${name}"
      fi
      mkdir -p ${sample_id}/"\${out_name}"
      bbmap.sh threads=${task.cpus} nodisk=t ref="\${genome}" \
        in=${fwd} in2=${rev} \
        scafstats=${sample_id}/"\${out_name}"/scafstats.tsv \
        covstats=${sample_id}/"\${out_name}"/covstats.tsv

      awk 'BEGIN{FS=OFS="\\t"} NR==1{print; next} {sub(/[[:space:]]*#.*/, "", \$1); print}' \
        ${sample_id}/"\${out_name}"/scafstats.tsv > ${sample_id}/"\${out_name}"/scafstats.clean.tsv
      mv ${sample_id}/"\${out_name}"/scafstats.clean.tsv ${sample_id}/"\${out_name}"/scafstats.tsv

      awk -v sample="${sample_id}" 'BEGIN{FS=OFS="\\t"} NR==1{for(i=1;i<=NF;i++) if(\$i=="assignedReads") assigned=i; next} assigned && \$1 !~ /^#/ {print sample, \$1, \$assigned}' \
        ${sample_id}/"\${out_name}"/scafstats.tsv >> ${sample_id}/${sample_id}_assignedReads.tsv
    done > ${sample_id}.log 2>&1 || {
      cat ${sample_id}.log >&2
      exit 1
    }

    touch ${sample_id}/.done
    """
}

process MERGE_ASSIGNED_READS {
    tag 'assigned_reads'
    label 'python'

    publishDir "${params.outdir}/expression", mode: 'copy', overwrite: true, pattern: 'assignedReads_feature_table.tsv'

    input:
    path assigned_tables

    output:
    path 'assignedReads_feature_table.tsv', emit: table

    script:
    """
    python3 "${projectDir}/bin/merge_assigned_reads.py" \
      --output assignedReads_feature_table.tsv \
      ${assigned_tables}
    """
}

process BUILD_GENE_LENGTH_TABLE {
    tag 'gene_lengths'
    label 'python'

    publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true, pattern: 'gene_lengths.tsv'

    input:
    path genomes

    output:
    path 'gene_lengths.tsv', emit: table

    script:
    """
    python3 "${projectDir}/bin/build_gene_length_table.py" \
      --output gene_lengths.tsv \
      ${genomes}
    """
}

workflow {
    samplesheet_ch = existingPathChannel(params.samples)
    VALIDATE_CONFIG(samplesheet_ch)
    ref_dir = params.ref_dir ?: params.drep_genomes_dir
    CHECK_DEPENDENCIES(ref_dir ? true : false)

    sample_rows = readSamples(params.samples)
    samples_ch = Channel
        .fromList(sample_rows)
        .map { row -> tuple(row.sample_id, row.condition ?: '', file(row.read1, checkIfExists: true), file(row.read2, checkIfExists: true)) }

    if (params.adapter_fasta) {
        adapter_ch = existingPathChannel(params.adapter_fasta)
    } else {
        PREPARE_TRIMMOMATIC_ADAPTER()
        adapter_ch = PREPARE_TRIMMOMATIC_ADAPTER.out.adapter
    }

    if (params.sortmerna_ref && file(params.sortmerna_ref).exists()) {
        sortmerna_ref_ch = existingPathChannel(params.sortmerna_ref)
    } else {
        DOWNLOAD_SORTMERNA_DATABASE()
        sortmerna_ref_ch = DOWNLOAD_SORTMERNA_DATABASE.out.ref
    }

    if (params.sortmerna_index_dir) {
        sortmerna_index_ch = existingPathChannel(params.sortmerna_index_dir)
    } else if (file("${params.sortmerna_dir}/index/idx").exists()) {
        sortmerna_index_ch = existingPathChannel("${params.sortmerna_dir}/index")
    } else if (file("${params.sortmerna_dir}/index/sortmerna_index/idx").exists()) {
        sortmerna_index_ch = existingPathChannel("${params.sortmerna_dir}/index/sortmerna_index")
    } else {
        PREPARE_SORTMERNA_INDEX(sortmerna_ref_ch)
        sortmerna_index_ch = PREPARE_SORTMERNA_INDEX.out.index
    }

    trim_input_ch = samples_ch.combine(adapter_ch)
    TRIM_RNA(trim_input_ch)

    FASTQC_RNA_TRIMMED(TRIM_RNA.out.trimmed)

    sortmerna_input_ch = TRIM_RNA.out.trimmed.combine(sortmerna_ref_ch).combine(sortmerna_index_ch)
    SORTMERNA(sortmerna_input_ch)

    if (ref_dir) {
        genomes_ch = Channel
            .fromPath("${ref_dir}/*.{fa,fna,fasta}", checkIfExists: true)
            .collect()

        bbmap_input_ch = SORTMERNA.out.rrna_removed
            .combine(genomes_ch)
            .map { values ->
                def genome_files = values[4..-1]
                tuple(values[0], values[1], values[2], values[3], genome_files)
            }
        BBMAP_EXPRESSION(bbmap_input_ch)
        BUILD_GENE_LENGTH_TABLE(genomes_ch)

        assigned_tables_ch = BBMAP_EXPRESSION.out.assigned_reads
            .map { sample_id, table -> table }
            .collect()
        MERGE_ASSIGNED_READS(assigned_tables_ch)
    }
}
