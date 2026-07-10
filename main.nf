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

    publishDir "${params.sortmerna_dir}/index", mode: 'copy', overwrite: true, pattern: 'sortmerna_index/**'
    publishDir 'logs/resources', mode: 'copy', overwrite: true, pattern: 'sortmerna_index.log'

    input:
    path ref

    output:
    path 'sortmerna_index', emit: index
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
      --workdir sortmerna_index \
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

    test -d sortmerna_index/idx
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
    cp -a ${index_dir} sortmerna_work_${sample_id}

    sortmerna \
      --workdir sortmerna_work_${sample_id} \
      --ref ${ref} \
      --reads ${fwd} \
      --reads ${rev} \
      --aligned ${sample_id} \
      --other ${sample_id}_rRNArm \
      --fastx \
      --threads ${task.cpus} \
      --out2 \
      --sout \
      > ${sample_id}.log 2>&1 || {
      cat ${sample_id}.log >&2
      exit 1
    }

    for suffix in fwd rev; do
      if [ -s "${sample_id}_rRNArm_\${suffix}.fq" ] && [ ! -s "${sample_id}_rRNArm_\${suffix}.fq.gz" ]; then
        gzip "${sample_id}_rRNArm_\${suffix}.fq"
      fi
    done

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
    path "${sample_id}/**", emit: expression
    path "${sample_id}.log", emit: log

    script:
    """
    mkdir -p ${sample_id}

    for genome in ${genomes}; do
      name=\$(basename "\${genome}")
      name="\${name%.*}"
      mkdir -p ${sample_id}/"\${name}"
      bbmap.sh threads=${task.cpus} nodisk=t ref="\${genome}" \
        in=${fwd} in2=${rev} \
        scafstats=${sample_id}/"\${name}"/scafstats.tsv \
        covstats=${sample_id}/"\${name}"/covstats.tsv
    done > ${sample_id}.log 2>&1

    touch ${sample_id}/.done
    """
}

workflow {
    samplesheet_ch = existingPathChannel(params.samples)
    VALIDATE_CONFIG(samplesheet_ch)
    CHECK_DEPENDENCIES(params.drep_genomes_dir ? true : false)

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
    } else {
        PREPARE_SORTMERNA_INDEX(sortmerna_ref_ch)
        sortmerna_index_ch = PREPARE_SORTMERNA_INDEX.out.index
    }

    trim_input_ch = samples_ch.combine(adapter_ch)
    TRIM_RNA(trim_input_ch)

    FASTQC_RNA_TRIMMED(TRIM_RNA.out.trimmed)

    sortmerna_input_ch = TRIM_RNA.out.trimmed.combine(sortmerna_ref_ch).combine(sortmerna_index_ch)
    SORTMERNA(sortmerna_input_ch)

    if (params.drep_genomes_dir) {
        genomes_ch = Channel
            .fromPath("${params.drep_genomes_dir}/*.{fa,fna,fasta}", checkIfExists: true)
            .collect()

        bbmap_input_ch = SORTMERNA.out.rrna_removed.combine(genomes_ch)
        BBMAP_EXPRESSION(bbmap_input_ch)
    }
}
