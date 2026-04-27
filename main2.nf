// main2.nf
// Nextflow Pipeline for Phage Isolate Assembly & Annotation
// Workflow: QC -> Flye -> Medaka -> CheckV/ViralVerify -> Pharokka -> Report
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
    echo "<html><body>Stub QC Report</body></html>" > fastplong.html
    echo '{"summary": "stub"}' > fastplong.json
    """
}

// ============================================================
// 2. Assembly: Flye 
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
    publishDir "${params.outdir}/03_qc", mode: 'copy'

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
    publishDir "${params.outdir}/03_qc", mode: 'copy'

    input:
    path polished_assembly

    output:
    path "viralverify_out", emit: viralverify_dir

    script:
    """
    viralverify -f ${polished_assembly} -o viralverify_out --hmm ~/micromamba/envs/phage-pipeline/db/nbc_hmms.hmm
    """

    stub:
    """
    mkdir -p viralverify_out
    printf 'Contig\\tPredicted_Label\\tFamily\\tGenus\\tviral_score\\n' > viralverify_out/viralverify_report.txt
    printf 'contig_1\\tVirus\\tSiphoviridae\\tUnknown\\t0.98\\n' >> viralverify_out/viralverify_report.txt
    """
}

// ============================================================
// 6. Phage Annotation: Pharokka
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
// 7. Report Generation (HTML & Slideshow Markdown)
// ============================================================
process GENERATE_REPORT {
    publishDir "${params.outdir}/05_final_report", mode: 'copy'

    input:
    path qc_html from qc_ch.qc_report
    path qc_json from qc_ch.qc_report_json
    path asm_info from asm_ch.assembly_info
    path checkv_dir from checkv_ch.checkv_results
    path verify_dir from verify_ch.viralverify_dir
    path pharokka_dir from pharokka_ch.pharokka_dir

    output:
    path "report.html", emit: html_report
    path "slideshow.md", emit: slide_deck
    path "assets/*", optional: true, emit: assets

    script:
    """
    # Create assets directory for any potential images/css if needed later
    mkdir -p assets

    # --- Parse Data (Simplified for stub/real compatibility) ---
    # In a real scenario, you would use python/R to parse TSV/JSON properly.
    # Here we use bash/awk for demonstration.
    
    QC_SUMMARY="QC Completed (See HTML details)"
    if [ -f ${qc_json} ]; then
        QC_SUMMARY=\$(cat ${qc_json} | head -c 200)
    fi

    ASM_SUMMARY="Assembly Completed"
    if [ -f ${asm_info} ]; then
        ASM_SUMMARY=\$(head -n 5 ${asm_info})
    fi

    CHECKV_SUMMARY="CheckV analysis available in folder"
    if [ -f ${checkv_dir}/quality_summary.tsv ]; then
        CHECKV_SUMMARY=\$(cat ${checkv_dir}/quality_summary.tsv)
    fi

    VERIFY_SUMMARY="ViralVerify analysis available"
    if [ -f ${verify_dir}/viralverify_report.txt ]; then
        VERIFY_SUMMARY=\$(cat ${verify_dir}/viralverify_report.txt)
    fi

    PHAROKKA_SUMMARY="Pharokka annotation completed"
    if [ -f ${pharokka_dir}/pharokka.gff ]; then
        CDS_COUNT=\$(grep -c "CDS" ${pharokka_dir}/pharokka.gff || echo "0")
        PHAROKKA_SUMMARY="Identified \$CDS_COUNT CDS features"
    fi

    # --- Generate HTML Report ---
    cat << EOF > report.html
    <!DOCTYPE html>
    <html>
    <head>
        <title>Phage Analysis Report: ${params.project_name}</title>
        <style>
            body { font-family: Helvetica, Arial, sans-serif; max-width: 900px; margin: 40px auto; line-height: 1.6; color: #333; }
            h1 { color: #2c3e50; border-bottom: 2px solid #eee; padding-bottom: 10px; }
            h2 { color: #2980b9; margin-top: 30px; }
            pre { background: #f4f4f4; padding: 15px; border-radius: 5px; overflow-x: auto; }
            table { width: 100%; border-collapse: collapse; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            .status { color: green; font-weight: bold; }
        </style>
    </head>
    <body>
        <h1>Phage Isolate Report: ${params.project_name}</h1>
        <p>Generated on: \$(date)</p>

        <h2>1. Quality Control (Fastplong)</h2>
        <p class="status">Status: Completed</p>
        <pre>$QC_SUMMARY</pre>

        <h2>2. Assembly (Flye & Medaka)</h2>
        <pre>$ASM_SUMMARY</pre>

        <h2>3. Quality Assessment (CheckV)</h2>
        <pre>$CHECKV_SUMMARY</pre>

        <h2>4. Taxonomic Verification (ViralVerify)</h2>
        <pre>$VERIFY_SUMMARY</pre>

        <h2>5. Annotation (Pharokka)</h2>
        <p class="status">$PHAROKKA_SUMMARY</p>
        
        <h2>6. Files Generated</h2>
        <ul>
            <li>Trimmed Reads: 01_qc/</li>
            <li>Assembly: 02_assembly/</li>
            <li>QC Metrics: 03_qc/</li>
            <li>Annotation: 04_annotation/</li>
        </ul>
        
        <p><i>To save as PDF: Use your browser's Print function (Ctrl+P) and select "Save as PDF".</i></p>
    </body>
    </html>
    EOF

    # --- Generate Slideshow Markdown (Marp/Reveal.js compatible) ---
    cat << EOF > slideshow.md
    ---
    marp: true
    theme: gaia
    class: lead
    ---

    # Phage Analysis Report
    ## ${params.project_name}
    
    Automated Pipeline Results
    
    ---

    # 1. Quality Control
    
    - Tool: Fastplong
    - Status: Completed
    - Summary: $QC_SUMMARY

    ---

    # 2. Assembly & Polishing
    
    - Assembler: Flye (Meta mode)
    - Polisher: Medaka
    - Result: $ASM_SUMMARY

    ---

    # 3. Quality & Verification
    
    ### CheckV
    $CHECKV_SUMMARY

    ### ViralVerify
    $VERIFY_SUMMARY

    ---

    # 4. Annotation
    
    - Tool: Pharokka
    - Result: $PHAROKKA_SUMMARY
    
    ---

    # Summary & Next Steps
    
    - All pipeline steps completed successfully.
    - Detailed files available in the 'results' directory.
    - Review HTML report for full tables.
    EOF

    echo "Report generation complete."
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
    qc_ch      = QC_TRIM(reads_ch)

    // 2. Assembly
    asm_ch      = ASSEMBLY(qc_ch.trimmed_reads)

    // 3. Polishing
    polish_ch   = POLISH(asm_ch.assembly, qc_ch.trimmed_reads)

    // 4. QC & Verification (Parallel)
    checkv_ch   = CHECKV(polish_ch.polished_assembly)
    verify_ch   = VIRAL_VERIFY(polish_ch.polished_assembly)

    // 5. Annotation
    pharokka_ch = PHAROKKA(polish_out.polished_assembly)

    // 6. Final Report Aggregation
    report_out   = GENERATE_REPORT(
        qc_ch.qc_report,
        qc_ch.qc_report_json,
        asm_ch.assembly_info,
        checkv_ch.checkv_results,
        verify_ch.viralverify_dir,
        pharokka_ch.pharokka_dir
    )

    report_out.html_report.view { path ->
        """
        =========================================
        Pipeline completed!
        HTML Report: ${path}
        Slideshow: ${report_out.slide_deck.first()}
        =========================================
        """
    }
}