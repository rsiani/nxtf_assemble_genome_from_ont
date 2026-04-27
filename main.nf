// main.nf
// Nextflow Pipeline for Phage Assembly
// Usage: 
//   Real run:  nextflow run main.nf --inputFile reads.fastq
//   Test run:  nextflow run main.nf --inputFile reads.fastq -stub

// --- Parameters ---
params.inputFile = null
params.outdir = './results'

// --- 1. QC & Trimming: fastplong ---
process QC_TRIM {
    input:
    path reads
    
    output:
    path "*.fastq.gz", emit: trimmed_reads
    path "fastplong_report.html", emit: qc_report
    path "fastplong_report.json", emit: qc_report_json

    script:
    """
    fastplong -i ${reads} -o ${reads.baseName}_trimmed.fastq.gz --n_base_limit 500 --thread ${task.cpus}
    """

    stub:
    """
    echo "[STUB] Running fastplong on ${reads}"
    echo "Dummy trimmed content" > ${reads.baseName}_trimmed.fastq.gz
    echo "<html><body>Stub QC Report</body></html>" > fastplong_report.html
    echo '{"summary": "stub"}' > fastplong_report.json
    """
}

// --- 2. Assembly: metaFlye ---
process ASSEMBLY {
    input:
    path trimmed_reads
    
    output:
    path "flye_out/assembly.fasta", emit: assembly
    path "flye_out/assembly_info.txt", emit: assembly_info

    script:
    """
    flye --nano-raw ${trimmed_reads} --meta --out-dir flye_out --threads ${task.cpus} --genome-size 50k
    """

    stub:
    """
    echo "[STUB] Running metaFlye assembly"
    mkdir -p flye_out
    echo ">contig_1 length=50000" > flye_out/assembly.fasta
    echo "ATCGN" >> flye_out/assembly.fasta
    echo "assembly_info: stubbed" > flye_out/assembly_info.txt
    """
}

// --- 3. Viral Refinement: viralFlye ---
process VIRAL_FLYE {
    input:
    path flye_dir
    path trimmed_reads
    
    output:
    path "viralflye_out/components_viralFlye.fasta", emit: final_assembly
    path "viralflye_out/host_prediction.txt", emit: host_pred
    path "viralflye_out/circulars_viralFlye.fasta", emit: circular_contigs

    script:
    """
    if [ ! -f Pfam-A.hmm.gz ]; then
        echo "Downloading Pfam-A.hmm.gz..."
        wget -q ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    fi
    viralFlye.py --dir ${flye_dir} --hmm Pfam-A.hmm.gz --reads ${trimmed_reads} --outdir viralflye_out
    """

    stub:
    """
    echo "[STUB] Running viralFlye refinement"
    mkdir -p viralflye_out
    echo ">viral_contig_1 circular=true" > viralflye_out/components_viralFlye.fasta
    echo "ATCGN" >> viralflye_out/components_viralFlye.fasta
    echo "Host: Escherichia coli (predicted)" > viralflye_out/host_prediction.txt
    cp viralflye_out/components_viralFlye.fasta viralflye_out/circulars_viralFlye.fasta
    """
}

// --- 4. Validation: CheckV ---
process CHECKV {
    input:
    path final_assembly
    
    output:
    path "checkv_out/*", emit: checkv_results

    script:
    """
    checkv end-to-end ${final_assembly} checkv_out --threads ${task.cpus}
    """

    stub:
    """
    echo "[STUB] Running CheckV validation"
    mkdir -p checkv_out
    echo "sample_id,completeness,contamination" > checkv_out/summary.csv
    echo "sample_1,95.5,1.2" >> checkv_out/summary.csv
    """
}

// --- 5. Validation: BUSCO ---
process BUSCO {
    input:
    path final_assembly
    
    output:
    path "busco_out/short_summary*.txt", emit: busco_report

    script:
    """
    busco -i ${final_assembly} -l viral_odb10 -m geno -o busco_out --cpu ${task.cpus}
    """

    stub:
    """
    echo "[STUB] Running BUSCO validation"
    mkdir -p busco_out
    echo "# BUSCO stub report" > busco_out/short_summary.txt
    echo "C:98.5%[S:97.0%,D:1.5%],F:0.5%,M:1.0%,n:100" >> busco_out/short_summary.txt
    """
}

// --- 6. Reporting: MultiQC ---
process REPORT {
    input:
    path qc_report
    path qc_report_json
    path checkv_results
    path busco_report
    
    output:
    path "multiqc_out/multiqc_report.html", emit: final_report

    script:
    """
    multiqc . -o multiqc_out
    """

    stub:
    """
    echo "[STUB] Generating MultiQC report"
    mkdir -p multiqc_out
    echo "<html><body><h1>Stub MultiQC Report</h1><p>All steps completed successfully.</p></body></html>" > multiqc_out/multiqc_report.html
    """
}

// --- Workflow Definition ---
workflow {
    // 1. Validate Input Parameter exists
    if (!params.inputFile) {
        error "Please provide an input file using --inputFile <path_to_fastq>"
    }

    // 2. Define Input Channel
    // 'checkIfExists: true' ensures we fail early if the path is wrong
    Channel.fromPath(params.inputFile, checkIfExists: true)
        .set { samples }

    // 3. Invoke Processes
    qc_results = QC_TRIM(samples)
    trimmed_reads_ch = qc_results.trimmed_reads
    qc_report_ch = qc_results.qc_report
    qc_json_ch = qc_results.qc_report_json

    asm_results = ASSEMBLY(trimmed_reads_ch)
    assembly_ch = asm_results.assembly

    viral_results = VIRAL_FLYE(assembly_ch, trimmed_reads_ch)
    final_assembly_ch = viral_results.final_assembly

    checkv_results = CHECKV(final_assembly_ch)
    checkv_out_ch = checkv_results.checkv_results

    busco_results = BUSCO(final_assembly_ch)
    busco_report_ch = busco_results.busco_report

    report_results = REPORT(qc_report_ch, qc_json_ch, checkv_out_ch, busco_report_ch)
    
    // 4. Final Output
    report_results.final_report.view { report_path -> 
        println "========================================="
        println "Pipeline completed successfully!"
        println "Final Report: $report_path"
        println "========================================="
    }
}