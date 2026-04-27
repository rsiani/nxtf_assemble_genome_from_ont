// main.nf
// Nextflow Pipeline for Phage Assembly
// Usage: 
//   Real run:  nextflow run main.nf --inputFile reads.fastq -profile conda
//   Test run:  nextflow run main.nf --inputFile reads.fastq -profile conda -stub

// --- Parameters ---
params.inputFile = null
params.outdir = './results'

// --- Process Configuration ---
process {
    // Enable conda
    conda = true
    condaEnvironment = "${projectDir}/environment.yml"
    
    // Default resources
    cpus = 4
    memory = '8 GB'
    
    // Error handling
    errorStrategy = 'finish'
    maxRetries = 1
}

// --- 1. QC & Trimming: fastplong ---
process QC_TRIM {
    tag "$sample"
    input:
    path reads
    
    output:
    path "*.fastq.gz", emit: trimmed_reads
    path "fastplong_report.html", emit: qc_report
    path "fastplong_report.json", emit: qc_report_json

    script:
    def sample = reads.simpleName
    """
    fastplong -i ${reads} -o ${sample}_trimmed.fastq -A --min_len 500 --threads ${task.cpus}
    gzip -c ${sample}_trimmed.fastq > ${sample}_trimmed.fastq.gz
    """

    stub:
    def sample = reads.simpleName
    """
    echo "[STUB] Running fastplong on ${reads}"
    echo "Dummy trimmed content" > ${sample}_trimmed.fastq
    gzip -c ${sample}_trimmed.fastq > ${sample}_trimmed.fastq.gz
    echo "<html><body>Stub QC Report</body></html>" > fastplong_report.html
    echo '{"summary": "stub"}' > fastplong_report.json
    """
}

// --- 2. Assembly: metaFlye ---
process ASSEMBLY {
    tag "$sample"
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
    tag "$sample"
    input:
    path flye_dir
    path reads
    
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
    viralFlye.py --dir ${flye_dir} --hmm Pfam-A.hmm.gz --reads ${reads} --outdir viralflye_out
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
    tag "$sample"
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
    tag "$sample"
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
    tag "$sample"
    input:
    path qc_report
    path qc_report_json
    path checkv_results
    path busco_report
    
    output:
    path "multiqc_report.html", emit: final_report

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
    if (!params.inputFile) {
        error "Please provide an input file using --inputFile <path_to_fastq>"
    }

    // Channel from input file
    Channel.fromPath(params.inputFile)
        .map { file -> [file.simpleName, file] }
        .set { samples }

    // Execute processes
    samples { |sample, reads|
        QC_TRIM(reads)
        ASSEMBLY(QC_TRIM.trimmed_reads)
        VIRAL_FLYE(ASSEMBLY.assembly, QC_TRIM.trimmed_reads)
        CHECKV(VIRAL_FLYE.final_assembly)
        BUSCO(VIRAL_FLYE.final_assembly)
        REPORT(QC_TRIM.qc_report, QC_TRIM.qc_report_json, CHECKV.checkv_results, BUSCO.busco_report)
        
        // Optional: Print final report path to console
        REPORT.final_report.view { "Pipeline finished! Report available at: $it" }
    }
}