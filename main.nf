// main.nf
// Nextflow Pipeline for Phage Isolate Assembly & Annotation
// Workflow: QC -> Flye -> Medaka -> CheckV/ViralVerify -> Pharokka -> Phold -> Phynteny -> Report
// Usage:
//   Real run:  nextflow run main.nf -c nextflow.config --inputFile long_reads.fastq 
//   Test run:  nextflow run main.nf -c nextflow.config --inputFile long_reads.fastq -stub

// --- Parameters ---
params.inputFile  = null
params.outdir     = './results'

// ============================================================
// 1. QC & Trimming: fastplong
// ============================================================
process QC_TRIM {
    publishDir "${params.outdir}/01_qc", mode: 'copy'

    input:
    path reads

    output:
    path "*_trimmed.fastq.gz", emit: trimmed_reads
    path "fastplong.html",  emit: qc_report
    path "fastplong.json",  emit: qc_report_json

    script:
    """
    fastplong \
        -i ${reads} \
        -o ${reads.baseName}_trimmed.fastq.gz \
        --n_base_limit 500 \
        --thread ${task.cpus}
    """

    stub:
    """
    echo "Dummy trimmed content" | gzip > ${reads.baseName}_trimmed.fastq.gz
    echo "<html><body>Stub QC Report</body></html>" > fastplong_report.html
    echo '{"summary": "stub"}' > fastplong_report.json
    """
}

// ============================================================
// 2. Assembly: Flye (Standard Mode)
// ============================================================
process ASSEMBLY {
    publishDir "${params.outdir}/02_assembly", mode: 'copy'

    input:
    path trimmed_reads

    output:
    path "flye_out/assembly.fasta",       emit: assembly
    path "flye_out/assembly_info.txt",    emit: assembly_info

    script:
    """
    flye \
        --nano-corr ${trimmed_reads} \
        --out-dir flye_out \
        --meta \
        --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p flye_out
    printf '>contig_1 length=50000\\nATCGATCGATCG\\n' > flye_out/assembly.fasta
    echo "assembly_info: stubbed" > flye_out/assembly_info.txt
    """
}

// ============================================================
// 3. Polishing: Medaka
// ============================================================
process POLISH {
    publishDir "${params.outdir}/02_assembly", mode: 'copy'

    input:
    path assembly
    path trimmed_reads

    output:
    path "medaka_out/consensus.fasta", emit: polished_assembly

    script:
    """
    medaka_consensus \
        -i ${trimmed_reads} \
        -d ${assembly} \
        -o medaka_out \
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p medaka_out
    printf '>contig_1_polished length=50000\\nATCGATCGATCG\\n' > medaka_out/consensus.fasta
    """
}

// ============================================================
// 4. Assembly QC: CheckV
// ============================================================
process CHECKV {
    publishDir "${params.outdir}/03_qc/checkv", mode: 'copy'

    input:
    path polished_assembly

    output:
    path "checkv_out/*", emit: checkv_results

    script:
    """
    checkv end_to_end ${polished_assembly} checkv_out -t ${task.cpus} -d ~/micromamba/envs/phage-pipeline/db/checkv-db-v1.5
    """

    stub:
    """
    mkdir -p checkv_out
    printf 'contig_id,completeness,contamination\\nviral_contig_1,95.5,1.2\\n' > checkv_out/quality_summary.tsv
    """
}

// ============================================================
// 5. Taxonomic Verification: ViralVerify
// ============================================================
process VIRAL_VERIFY {
    publishDir "${params.outdir}/03_qc/viralverify", mode: 'copy'

    input:
    path polished_assembly

    output:
    path "viralverify_report.txt", emit: taxonomy_report

    script:
    """
    viralverify -f ${polished_assembly} -o viralverify_report.txt --hmm ~/micromamba/envs/phage-pipeline/db/nbc_hmms.hmm
    """

    stub:
    """
    printf 'Contig\\tPredicted_Label\\tFamily\\tGenus\\tviral_score\\n' > viralverify_report.txt
    printf 'contig_1\\tVirus\\tSiphoviridae\\tUnknown\\t0.98\\n' >> viralverify_report.txt
    """
}

// ============================================================
// 6a. Phage Annotation: Pharokka
// ============================================================
process PHAROKKA {
    publishDir "${params.outdir}/04_annotation/", mode: 'copy'

    input:
    path polished_assembly

    output:
    path "pharokka_out/pharokka.gff",         emit: pharokka_gff
    path "pharokka_out/pharokka.gbk",         emit: pharokka_gbk
    path "pharokka_out",                      emit: pharokka_dir

    script:
    """
    pharokka.py \
        -i ${polished_assembly} \
        -o pharokka_out \
        -t ${task.cpus} \
        -d ~/micromamba/envs/phage-pipeline/db/ \
        -f
    """

    stub:
    """
    mkdir -p pharokka_out
    printf '##gff-version 3\\nviral_contig_1\\tpharokka\\tCDS\\t1\\t500\\t.\\t+\\t0\\tID=CDS_1;phrog=phrog_1;function=tail fiber\\n' \
        > pharokka_out/pharokka.gff
    echo "LOCUS viral_contig_1" > pharokka_out/pharokka.gbk
    """
}

// ============================================================
// 6b. Structure-based Annotation: Phold
// ============================================================
process PHOLD {
    publishDir "${params.outdir}/04_annotation/phold", mode: 'copy'

    input:
    path gbk_file

    output:
    path "phold_out/phold.gff",         emit: phold_gff
    path "phold_out/phold.gbk",         emit: phold_gbk
    path "phold_out",                   emit: phold_dir

    script:
    """
    phold run \
        -i ${gbk_file} \
        -o phold_out \
        -t ${task.cpus} \
        -f
    """

    stub:
    """
    mkdir -p phold_out
    printf '##gff-version 3\\nviral_contig_1\\tphold\\tCDS\\t1\\t500\\t.\\t+\\t0\\tID=CDS_1;function=tail fiber protein (3D)\\n' \
        > phold_out/phold.gff
    cp ${gbk_file} phold_out/phold.gbk
    """
}

// // ============================================================
// // 6c. Synteny-based Context: Phynteny
// // ============================================================
// process PHYNTENY {
//     publishDir "${params.outdir}/04_annotation/phynteny", mode: 'copy'

//     input:
//     path gbk_file

//     output:
//     path "phynteny_out/phynteny.gbk",         emit: phynteny_gbk
//     path "phynteny_out/phynteny_summary.tsv", emit: phynteny_summary

//     script:
//     """
//     mkdir -p phynteny_input
//     cp ${gbk_file} phynteny_input/
    
//     phynteny \
//         phynteny_input \
//         phynteny_out
//     """

//     stub:
//     """
//     mkdir -p phynteny_out
//     echo "LOCUS viral_contig_1_phynteny" > phynteny_out/phynteny.gbk
//     printf 'contig\\tsynteny_group\\trelated_phages\\ncontig_1\\tSG_042\\t12\\n' \
//         > phynteny_out/phynteny_summary.tsv
//     """
// }

// ============================================================
// 6. Report: MultiQC
// ============================================================
process REPORT {
    tag "multiqc"
    publishDir "${params.outdir}/05_report", mode: 'copy'

    input:
    path qc_report
    path qc_report_json
    path checkv_results
    path taxonomy_report  
    // path phynteny_summary

    output:
    path "multiqc_out/multiqc_report.html", emit: final_report

    script:
    """
    multiqc . -o multiqc_out
    """

    stub:
    """
    mkdir -p multiqc_out
    cat <<'EOF' > multiqc_out/multiqc_report.html
    <html>
    <body>
        <h1>Stub MultiQC Report</h1>
        <p>All steps completed successfully (stub mode).</p>
        <p>Includes: QC, CheckV, ViralVerify, Pharokka, Phold, Phynteny</p>
    </body>
    </html>
    EOF
    """
}

// ============================================================
// Workflow
// ============================================================
workflow {

    if (!params.inputFile) {
        error "Please provide an input file: --inputFile <path_to_fastq>"
    }

    reads_ch = Channel.fromPath(params.inputFile, checkIfExists: true)

    // 1. QC
    qc_out       = QC_TRIM(reads_ch)

    // 2. Assembly
    asm_out      = ASSEMBLY(qc_out.trimmed_reads)

    // 3. Polishing
    polish_out   = POLISH(asm_out.assembly, qc_out.trimmed_reads)

    // 4. QC & Verification (Parallel)
    checkv_out   = CHECKV(polish_out.polished_assembly)
    verify_out   = VIRAL_VERIFY(polish_out.polished_assembly)

    // 5. Annotation Chain
    pharokka_out = PHAROKKA(polish_out.polished_assembly)
    phold_out    = PHOLD(pharokka_out.pharokka_gbk)
    // phynteny_out = PHYNTENY(phold_out.phold_gbk)

    // 6. Report
    report_out   = REPORT(
        qc_out.qc_report,
        qc_out.qc_report_json,
        checkv_out.checkv_results,
        verify_out.taxonomy_report
    )

    report_out.final_report.view { path ->
        """
        =========================================
        Pipeline completed!
        Report: ${path}
        =========================================
        """
    }
}